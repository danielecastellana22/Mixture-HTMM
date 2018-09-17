function test_on_inex05_INF(max_it, C, alpha_A, alpha_pr, alpha_SP, alpha_b, gamma, n_worker, start_tr, end_tr)
    %LOAD THE LIBRARY------------------------------------------------------
    BRML_dir = './model/';
    p = genpath(BRML_dir);
    addpath(p); 

    % load dataset
    load('./data/inex05/inex05.mat');
    
    par_logdir = './logs/inex05';
    
    parpool(n_worker);
    
    for i= start_tr:end_tr
        %create all log dir
        logdir = [par_logdir '/INF_' num2str(alpha_A) '_' num2str(C) '_' num2str(gamma) '/trial_' num2str(i) ];
        mkdir(logdir)
        
        m = INF_SP_BHTMM(gamma,C,L,M,alpha_A,alpha_pr,alpha_b,alpha_SP, logdir);

        [likelihood_vect, time_vect, K_vect, z_tr] = m.train(Xtr,max_it);

        [z_test,log_test] = m.test(Xtest,10);

        save([logdir '/trial_' num2str(i) '.mat'])
    end
    
end