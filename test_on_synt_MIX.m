function test_on_synt_MIX(max_it, C, T, n_worker, start_tr, end_tr)
    %LOAD THE LIBRARY------------------------------------------------------
    BRML_dir = './model/';
    p = genpath(BRML_dir);
    addpath(p); 
    
    % load dataset
    load('./data/synthetic/synt.mat');
    
    par_logdir = './logs/synthetic';
    
    parpool(n_worker);
    
    for i= start_tr:end_tr
        %create all log dir
        logdir = [par_logdir '/MIX_' num2str(T) '_' num2str(C) '/trial_' num2str(i) ];
        mkdir(logdir)
        
        m = MIX_SP_BHTMM(C,T,L,M, logdir);

        [likelihood_vect, time_vect, K_vect, z_tr] = m.train(Xtr,max_it);

        [z_test,log_test] = m.test(Xtest,10);

        save([logdir '/trial_' num2str(i) '.mat'])
    end
    
end