classdef MIX_SP_BHTMM < handle
    %MIX_SP_BHTMM class to represent MIX_SP_BHTMM
    
    properties
        % number of mixture
        K;
        % number of hidden states in each SP-BHMT
        C;
        % number of max child nodes
        L;
        % number of emisssion labels
        M;
        
        % mixture distribution
        mix_distribution;
        
        % cell-array for mixture compoenent
        mix_components;
        
        % accumulate statistics to update mix_distribution
        new_mix_tab;
        
        logdir;
    end
    
    methods
        function obj = MIX_SP_BHTMM(C_,T_,L_,M_, par_logdir)
            %DP_MIX_SP_BHMT Construct an instance of this class
            obj.C = C_;
            obj.L = L_;
            obj.M = M_;
            obj.K = T_;
           
            % create all components
            for k=1:T_
                obj.mix_components{k} = SP_BHTMM(C_,L_,M_,1,1,1,1);
            end
            
            % create the distribution P(Tp=k)
            obj.mix_distribution.variables = 1;
            obj.mix_distribution.table = log(drchrnd(ones(1,T_),1)');
            
            obj.new_mix_tab = zeros(T_,1) + eps;
            
            % create log_dir
            obj.logdir = par_logdir;
            obj.init_log();
        end
        
        function [likelihood_vect, time_vect, K_vect, z] = train(obj,X_list,max_it)
            
            likelihood_vect = zeros(max_it,1);
            time_vect = zeros(max_it,1);
            K_vect = zeros(max_it,1);
            
            fprintf(1,'Train\n');
            for it=1:max_it
                tic;
                [loglike_trees,z] = obj.compute_posteriors(X_list, false);
                loglike = sum(loglike_trees);
                obj.update_params();
                t = toc;
                
                % print on screen
                fprintf(1,'It. %3d\tLoglike: %.10f\n',it,loglike);
                % save log
                obj.save_log_it(it,t,loglike);
                
                likelihood_vect(it) = loglike;
                time_vect(it) = t;
                K_vect(it) = obj.K;
            end
        end
        
        function [z, loglikeVect] = test(obj,X_list,max_it)
            loglikeVect = zeros(max_it,1);
            % initialise randomly z
            
            fprintf(1,'Test\n');
            for it=1:max_it
                [loglike_trees,z] = obj.compute_posteriors(X_list, true);
                loglike = sum(loglike_trees);
                
                % print on screen
                fprintf(1,'It. %3d\tLoglike: %.10f\n',it,sum(loglike));
                loglikeVect(it) = sum(loglike);
            end
        end
        
        function [loglike_trees,z] = compute_posteriors(obj,X_list,is_test)
            n = length(X_list);
            
            z = zeros(n,1);
            loglike_trees = zeros(n,1);
            % sample the assignment z
            textprogressbar('Data:');
            for i=1:n
                textprogressbar(100*(i/n));
                X_tree = X_list{i};
                
                % compute likelihood for each components (i.e. P(X | Tp=k))
                lk = zeros(obj.K,1);
                mix_comp_ = obj.mix_components;
                parfor j=1:obj.K
                    BU_model = mix_comp_{j};
                    % lk si P(X | TP=j)
                    lk(j) = BU_model.upward(X_tree);
                    mix_comp_{j} = BU_model;
                end
                obj.mix_components = mix_comp_;
                
                
                % create P(X | TP=k)
                lk_pot.variables = 1;
                lk_pot.table = lk;
                
                % compute P(TP=k, X)
                joint_tp_x = mulpot(lk_pot, obj.mix_distribution);
                % compute P(X)
                Z = sumpot(joint_tp_x);
                loglike_trees(i) = Z.table;
                % compute P(Tp=k | X)
                sampling_class_pot = condpot(joint_tp_x);
                % accumulate statistics
                obj.new_mix_tab = obj.new_mix_tab + exp(sampling_class_pot.table);
                
                % assign the tree to a mixture components
                s = samplepot(sampling_class_pot, 1);
                               
                z(i) = s;
                if(~is_test)
                    mix_comp_ = obj.mix_components;
                    parfor j=1:obj.K
                        BU_model = mix_comp_{j};
                         BU_model.downward(X_tree, sampling_class_pot.table(j));
                        mix_comp_{j} = BU_model;
                    end
                    obj.mix_components = mix_comp_;
                end
            end
            
            textprogressbar('Done');
        end
        
        function update_params(obj)            
            % do M-step
            for j=1:obj.K
                obj.mix_components{j}.M_step()
            end
            
            % update the mix_distribution
            obj.mix_distribution.table = log(obj.new_mix_tab);
            obj.mix_distribution = condpot(obj.mix_distribution);
            
            obj.new_mix_tab = zeros(obj.K,1) + eps;
        end
         
    end
    
    methods(Access=private)
        
        function logdir = init_log(obj)
            
            logdir = obj.logdir;
            
            mkdir([logdir '/iter']);
            
            fs = fopen([logdir '/results.out'],'w');
            
            fprintf(fs,'============== MODEL SUMMARY ==============\n');
            fprintf(fs,'Mixture of SP-BHMT\n');
            fprintf(fs,'M: %d\nL: %d\nC: %d\nT: %d\n',obj.M,obj.L,obj.C,obj.K);
            
            fclose(fs);
        end
        
        function save_log_it(obj, it, t, loglike_val)
            save([ obj.logdir '/iter/' num2str(it,'%03d') '.mat'], 'obj');
            
            fs=fopen([obj.logdir '/results.out'],'a');
            fprintf(fs,'============== ITERATION %03d ==============\n',it);
            fprintf(fs,'Elapsed time: %.6f s\n', t);
            fprintf(fs,'Loglike: %.10f\nK: %d\n',loglike_val,obj.K);
            fclose(fs);
        end
        
    end
end

