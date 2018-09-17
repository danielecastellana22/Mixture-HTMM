classdef SP_BHTMM_eq < handle
    %INPUTMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %number of maximum output degree
        L;
        %number of hidden values
        C;
        %number of visibile values
        M;
        
        %variable for the potential
        h_state_ch = 1;
        h_state_pa = 2;
        v_state = 3;
        ch_pos = 5;
        
        %the potentials
        prLeaf;
        b;
        A;
        SP;
        
        %variables to update potentials
        newPrTab;
        newSPTab;
        newATab;
        newBTab;
        
        bottomCh;
        
        %beta vals
        betaVals;
        etaChVals;
        etaVals;
    end
    
    methods
        
        function obj = SP_BHTMM_eq(C_,L_,M_)
            
            %set the bottom valeus
            bottomCh_ = C_ + 1;
            C_ = C_ + 1;
            
            obj.C = C_;
            obj.L = L_;
            obj.M = M_;
            obj.bottomCh = bottomCh_;
            
            %init the potentials
            obj.init_potentials_from_priori();

            %variables to update potentials
            obj.init_counting_table();
        end
        
        function [ loglike ] = upward(obj,X_tree)
            %UPWARDINPUT This function executes the upward pass on the input tree X.
                        
            %read variables from obj
            L_ = obj.L;
            C_ = obj.C;
            M_ = obj.M;
            
            bottomCh_ = obj.bottomCh;
            
            h_state_ch_ = obj.h_state_ch;
            h_state_pa_ = obj.h_state_pa;
            v_state_ = obj.v_state;
            ch_pos_ = obj.ch_pos;

            prLeaf_ = obj.prLeaf;
            b_ = obj.b;
            A_ = obj.A;
            SP_ = obj.SP;
            
            X = X_tree.adj;
            n = X_tree.n;
            nf = X_tree.nf;
            v = X_tree.v;
            
            %betaVals contains potentials over the h_state_pa variable
            betaVals_ = cell(n,1);
            
            loglikePot.variables = [];
            loglikePot.table = zeros(1,1);
            
            %compute betaVals for the leaf
            for u=n-nf+1:n
                [~,~,l] = find(X(:,u));
                evid = setpot(b_,v_state_,v(u));
                pr = setpot(prLeaf_,ch_pos_,l);
                
                betaPot = mulpot(evid, pr);
                
                Nu = sumpot(betaPot);
                loglikePot = mulpot(loglikePot, Nu);
                
                betaVals_{u} = condpot(betaPot);
            end
            
            ccc=0;
            %compute betaVals for the internal nodes
            for u=n-nf:-1:1
                [~,idCh,posCh] = find(X(u,:));
                allCh = zeros(L_,1);
                allCh(posCh) = idCh;
                %table for beta childs
                betaChTab = -inf(C_,L_);
                for l=1:L_
                    ch = allCh(l);
                    if(ch~=0)
                        betaChTab(:,l) = betaVals_{ch}.table;
                    else
                        ccc = ccc+1;
                        betaChTab(bottomCh_,l) = 0;
                        %update the loglike
                        Nu = setpot(prLeaf_,[h_state_pa_ ch_pos_],[bottomCh_ l]);
                        loglikePot = mulpot(loglikePot, Nu);
                    end
                end
                
                betaChPot.variables = [h_state_ch_ ch_pos_];
                betaChPot.table = betaChTab;
                
                evid = setpot(b_,v_state_,v(u));
                
                numPot = sumpot(mulpot(SP_,mulpot( A_, betaChPot)),[h_state_ch_ ch_pos_]);
                
                betaPot =  mulpot(evid, numPot);
                
                Nu = sumpot(betaPot);
                loglikePot = mulpot(loglikePot, Nu);
                
                betaVals_{u} = condpot(betaPot,h_state_pa_);
            end
            
            obj.betaVals = betaVals_;
            
            loglike = loglikePot.table;
        end
    
        function [all_posterior] = downward(obj,X_tree, p_to_be_selected)
         	%DOWNARDIN Summary of this function goes here
            %   Detailed explanation goes here
            
            %read variables from obj
            L_ = obj.L;
            C_ = obj.C;
            M_ = obj.M;
            
            bottomCh_ = obj.bottomCh;
            
            h_state_ch_ = obj.h_state_ch;
            h_state_pa_ = obj.h_state_pa;
            v_state_ = obj.v_state;
            ch_pos_ = obj.ch_pos;
            
            prLeaf_ = obj.prLeaf;
            b_ = obj.b;
            A_ = obj.A;
            SP_ = obj.SP;     
            
            X = X_tree.adj;
            n = X_tree.n;
            nf = X_tree.nf;
            v = X_tree.v;
            
            betaVals_ = obj.betaVals;
            
            newPrTab_ = obj.newPrTab;
            newATab_ = obj.newATab;
            newBTab_ = obj.newBTab;
            newSPTab_ = obj.newSPTab;
            
            pot_to_be_selected.variables = [];
            pot_to_be_selected.table = p_to_be_selected;
            
            etaVals_ = cell(n,1);
            
            etaVals_{1} = betaVals_{1};
            
            %update the b
            sumTab = -inf(size(b_.table));
            sumTab(:,v(1)) = etaVals_{1}.table;
            sumPot.variables = b_.variables;
            sumPot.table = sumTab;
                        
            newBTab_ = pluspot(newBTab_ , mulpot(sumPot,pot_to_be_selected));
            
            %for each internal node
            for u=1:n-nf
                [~,idCh,posCh] = find(X(u,:));
                allCh = zeros(L_,1);
                allCh(posCh) = idCh;
                betaChTab = -inf(C_,L_);
                
                for l=1:L_
                    ch = allCh(l);
                    if(ch~=0)
                        betaChTab(:,l) = betaVals_{ch}.table;
                    else
                        betaChTab(bottomCh_,l) = 0;
                    end
                end
                
                betaChPot.variables = [h_state_ch_ ch_pos_];
                betaChPot.table = betaChTab;
                
                etaPa = etaVals_{u};
                
                numPot = mulpot(betaChPot,mulpot( SP_ , A_));
                denPot = sumpot(numPot,[h_state_ch_ ch_pos_]);
               
                %P(Q_ch = j, Q_pa = i, Sp = l | X)
                etaChVals_u = divpot(mulpot(etaPa, numPot),denPot);
                
                %P(Q_ch = j, Q_pa = i | Sp = l  X)
                etaChValsCond_u = condpot(etaChVals_u,[h_state_ch_ h_state_pa_], ch_pos_);
                
                obj.etaChVals{u} = etaChVals_u;
               
                %update A_ 
                % WHICH IS THE CORRECT ONE?--------------------------------
                %newATab_ = newATab_ + condpot(etaChVals_u,[h_state_ch_ h_state_pa_],[tp_state_ ch_pos_]);
                %newATab_ = newATab_ + condpot(jointPot,[h_state_ch_ h_state_pa_],[tp_state_ ch_pos_]);
                %newATab_ = newATab_ + etaChVals_u;
                newATab_ = pluspot(newATab_, mulpot(etaChValsCond_u, pot_to_be_selected));
                %----------------------------------------------------------
                %update SP_
                newSPTab_ = pluspot(newSPTab_, mulpot(sumpot(etaChVals_u,ch_pos_,0), pot_to_be_selected));
                
                %set eta childs
                %P(Q_ch = j | Sp = l, X)
                etaChPot = sumpot(etaChValsCond_u, h_state_pa_);
                for l=1:L_
                    ch = allCh(l);
                    
                    %P(Q_ch = j | Sp = l, X)
                    etaPot = setpot(etaChPot, ch_pos_,l);
                    
                    if(ch~=0)
                        
                        etaVals_{ch} = changevar(etaPot, h_state_ch_, h_state_pa_); 
                        
                        %update b_
                        sumTab = -inf(size(b_.table));
                        sumTab(:,v(ch)) = etaPot.table;
                        sumPot.variables = b_.variables;
                        sumPot.table = sumTab;
                        
                        newBTab_ = pluspot(newBTab_, mulpot(sumPot, pot_to_be_selected));
                        
                        
                        if(ch>n-nf)
                            %is a leaf -> update priori
                            sumTab = -inf(size(prLeaf_.table));
                            sumTab(:,l) = etaPot.table;
                            sumPot.variables = prLeaf_.variables;
                            sumPot.table = sumTab;
                            
                            newPrTab_ = pluspot(newPrTab_, mulpot(sumPot, pot_to_be_selected));
                        end
                    else
                        %is a Bottom child! update the prioiri
                        sumTab = -inf(size(prLeaf_.table));
                        sumTab(:,l) = etaPot.table;
                        sumPot.variables = prLeaf_.variables;
                        sumPot.table = sumTab;
                            
                        newPrTab_ = pluspot(newPrTab_, mulpot(sumPot, pot_to_be_selected));
                    end
                end
            end
            
            obj.newATab = newATab_;
            obj.newBTab = newBTab_;
            obj.newSPTab = newSPTab_;
            obj.newPrTab = newPrTab_;
            
            obj.etaVals = etaVals_;
            
            all_posterior = etaVals_;
        end
        
        function M_step(obj)
            
            h_state_ch_ = obj.h_state_ch;
            h_state_pa_ = obj.h_state_pa;
            v_state_ = obj.v_state;
            ch_pos_ = obj.ch_pos;
            
            
            %M-STEP PR_LEAF OK
            obj.prLeaf = condpot(obj.newPrTab,h_state_pa_, ch_pos_);           

            %M-STEP SP OK
            obj.SP = condpot(obj.newSPTab);            

            %M-STEP A OK
            obj.A = condpot(obj.newATab,h_state_pa_,[ h_state_ch_ ch_pos_]);

            %M-STEP b OK
            obj.b = condpot(obj.newBTab,v_state_,h_state_pa_);
             
            %newObj = obj;
            obj.init_counting_table();
        end
             
        function cleanUp(obj)
            obj.betaVals = [];
            obj.newATab = [];
            obj.newBTab = [];
            obj.newPrTab = [];
            obj.newSPTab = [];
        end
        
        function init_counting_table(obj)
            %variables to update potentials
            L_ = obj.L;
            C_ = obj.C;
            M_ = obj.M;
            
            bottomCh_ = obj.bottomCh;
           
            %bigEps = 10^-8;
            
            tab = zeros(C_,L_) + eps;
            obj.newPrTab.variables = obj.prLeaf.variables;
            obj.newPrTab.table = log(tab);
            
            tab = zeros(L_,1) + eps;
            obj.newSPTab.variables = obj.SP.variables;
            obj.newSPTab.table = log(tab);
            
            tab = zeros(C_,C_,L_) + eps;
            tab(:,bottomCh_,:,:) = zeros(C_,1,L_);
            obj.newATab.variables = obj.A.variables;
            obj.newATab.table = log(tab);
            
            tab = zeros(C_,M_) + eps;
            tab(bottomCh_,:) = zeros(1,M_);           
            obj.newBTab.variables = obj.b.variables;
            obj.newBTab.table = log(tab);
            
        end   
        
        function init_potentials_from_priori(obj)
            %variables to update potentials
            L_ = obj.L;
            C_ = obj.C;
            M_ = obj.M;
            
            h_state_ch_ = obj.h_state_ch;
            h_state_pa_ = obj.h_state_pa;
            v_state_ = obj.v_state;
            ch_pos_ = obj.ch_pos;
            
            bottomCh_ = obj.bottomCh;
            
            tab = 1/C_ * ones(C_,L_);
            obj.prLeaf.variables = [h_state_pa_ ch_pos_];
            obj.prLeaf.table = log(tab);
            
            
            tab = 1/L_ * ones(1,L_);
            obj.SP.variables = ch_pos_;
            % tab must be a column vector
            obj.SP.table = log(tab');
            
            tab = 1/C_ * ones(C_,C_,L_);
            tab(:,bottomCh_,:) = zeros(C_,1,L_);
            obj.A.variables = [h_state_ch_ h_state_pa_ ch_pos_];
            obj.A.table = log(tab);
            
            
            tab = 1/M_ * ones(C_,M_);
            tab(bottomCh_,:) = zeros(1,M_);           
            obj.b.variables = [h_state_pa_ v_state_];
            obj.b.table = log(tab);
            
        end
    end
end

