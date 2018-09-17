classdef INF_SP_BHTMM < handle
    %INF_SP_BHTMM class to represent INF-SP-BHTMM
    
    properties
        % number of mixture
        K;
        % number of hidden states in each SP-BHMT
        C;
        % number of max child nodes
        L;
        % number of emisssion labels
        M;
        
        % alpha dirichlet priori over tables
        alpha_A;
        alpha_pr;
        alpha_b;
        alpha_SP;
        
        eq_model;
        
        % hyper-parameters
        gamma;
        % cell-array for mixture compoenent
        mix_components;
        
        logdir;
    end
    
    methods
        function obj = INF_SP_BHTMM(gamma_,C_,L_,M_,alpha_A_,alpha_pr_,alpha_b_,alpha_SP_, par_logdir)
            %DP_MIX_SP_BHMT Construct an instance of this class
            obj.C = C_;
            obj.L = L_;
            obj.M = M_;
            obj.K = 1;
            
            obj.gamma = gamma_;
                     
            obj.alpha_A = alpha_A_;
            obj.alpha_pr = alpha_pr_;
            obj.alpha_b = alpha_b_;
            obj.alpha_SP = alpha_SP_;
            
            obj.mix_components{1} = SP_BHTMM(C_,L_,M_,alpha_A_,alpha_pr_,alpha_b_,alpha_SP_);
            
            obj.eq_model = SP_BHTMM_eq(obj.C,obj.L,obj.M);
            
            % create log_dir
            obj.logdir = par_logdir;
            obj.init_log();
        end
        
        function [likelihood_vect, time_vect, K_vect, z] = train(obj,X_list,max_it)
            
            likelihood_vect = zeros(max_it,1);
            time_vect = zeros(max_it,1);
            K_vect = zeros(max_it,1);
            
            n = length(X_list);
            % initialise randomly z
            z_in = randi(obj.K, n ,1);
            count_occ_in = zeros(obj.K, 1);
            for i=1:obj.K
                count_occ_in(i) = sum(z_in==i);
            end
            
            fprintf(1,'Train\n');
            for it=1:max_it
                tic;
                [loglike_trees,z,count_occ] = obj.compute_posteriors(X_list,z_in,count_occ_in, false);
                loglike = sum(loglike_trees);
                obj.update_params();
                % set the assignment obtained as initials state for the 
                % next iteration
                z_in = z;
                count_occ_in = count_occ;
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
            n = length(X_list);
            loglikeVect = zeros(max_it,1);
            % initialise randomly z
            z_in = randi(obj.K, n ,1);
            count_occ_in = zeros(obj.K, 1);
            for i=1:obj.K
                count_occ_in(i) = sum(z_in==i);
            end
            
            fprintf(1,'Test\n');
            for it=1:max_it
                [loglike_trees,z,count_occ] = obj.compute_posteriors(X_list,z_in,count_occ_in, true);
                loglike = sum(loglike_trees);
                % set the assignment obtained as initials state for the 
                % next iteration
                if(all(z_in==z))
                    % reached convergence
                    break;
                end
                z_in = z;
                count_occ_in = count_occ;
                
                % print on screen
                fprintf(1,'It. %3d\tLoglike: %.10f\n',it,sum(loglike));
                loglikeVect(it) = sum(loglike);
            end
        end
        
        function [loglike_trees,z,count_occ] = compute_posteriors(obj,X_list,z_in,count_occ_in,is_test)
            n = length(X_list);
            
            % assignment variables. z_i = j iff the i-th tree was
            % generated by the j-th mixture components
            z = z_in;
            count_occ = count_occ_in;
            loglike_trees = zeros(n,1);
            % sample the assignment z
            textprogressbar('Data:');
            for i=1:n
                textprogressbar(100*(i/n));
                X_tree = X_list{i};
                
                if(z(i)<=obj.K)
                    % remove my count
                    count_occ(z(i)) = count_occ(z(i)) - 1;
                end
                
                if (~is_test) && (count_occ(z(i))==0)
                    % remove mix_components no longer used
                    obj.mix_components(z(i)) = [];
                    count_occ(z(i)) = [];
                    % reduce K
                    obj.K = obj.K - 1;
                    % recompute assignments
                    idx = z >= z(i);
                    z(idx) = z(idx) - 1;
                end
                
                % +1 for the remaning stick
                lk = zeros(obj.K+1,1);
                mix_comp_ = obj.mix_components;
                parfor j=1:obj.K
                    BU_model = mix_comp_{j};
                    % lk si P(X | TP=j)
                    lk(j) = BU_model.upward(X_tree);
                    mix_comp_{j} = BU_model;
                end
                obj.mix_components = mix_comp_;
                
                % likelihood obtained integrating over all parameters
                % lk(obj.K+1) = X_tree.n * (log(1/M_));
                lk(obj.K+1) = obj.eq_model.upward(X_tree);
                
                
                % multuply lk to beat_stick in order to obtain P(TP=k | X)
                sampling_class_pot.variables = 1;
                sampling_class_pot.table = lk + log([count_occ; obj.gamma]);
                Z = sumpot(sampling_class_pot);
                loglike_trees(i) = Z.table;
                sampling_class_pot = condpot(sampling_class_pot);
                % assign the tree to a mixture components
                s = samplepot(sampling_class_pot, 1);
                               
                z(i) = s;
                
                if(~is_test)
                    if(s > obj.K)
                        % sample new mixture component
                        obj.mix_components{obj.K+1} = SP_BHTMM(obj.C,obj.L,obj.M,obj.alpha_A,obj.alpha_pr,obj.alpha_b,obj.alpha_SP);
                        obj.mix_components{obj.K+1}.upward(X_tree);
                        obj.K = obj.K + 1;
                        count_occ = [count_occ; 0];
                        
                    end                    
                    % do downward pass only on the assigned component
                    obj.mix_components{z(i)}.downward(X_tree, 0);
                end
                
                if(s<=obj.K)
                    count_occ(s) = count_occ(s) + 1;
                end
            end
            
            textprogressbar('Done');
        end
        
        function update_params(obj)            
            % do M-step
            for j=1:obj.K
                obj.mix_components{j}.M_step()
            end
        end
         
    end
    
    methods(Access=private)
        
        function logdir = init_log(obj)
            
            logdir = obj.logdir;
            
            mkdir([logdir '/iter']);
            
            fs = fopen([logdir '/results.out'],'w');
            
            fprintf(fs,'============== MODEL SUMMARY ==============\n');
            fprintf(fs,'Monte Carlo Markov Chain sampling\n');
            fprintf(fs,'M: %d\nL: %d\nC: %d\n',obj.M,obj.L,obj.C);
            fprintf(fs,'gamma: %.5f\nalpha_A: %.5f\nalpha_pr: %.5f\nalpha_SP: %.5f\nalpha_b: %.5f\n',obj.gamma, obj.alpha_A, obj.alpha_pr, obj.alpha_SP, obj.alpha_b);
            
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

