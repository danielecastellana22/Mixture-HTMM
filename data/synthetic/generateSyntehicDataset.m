function generateSyntehicDataset()
    BRML_dir = '../../model/';
    p = genpath(BRML_dir);
    addpath(p); 
    
    % to create the trees
    L = 3;
    M = 5;
    nMin = [20 50 80];
    nMax = [25 60 100];
    labelType = 1;

    % data set size
    Ntr = 600;
    Ntest = 180;

    %% create the training set
    Xtr = cell(Ntr,1);
    dd = Ntr/3;
    for i=1:3
        Xtr((i-1)*dd+1:i*dd) = treeGen(dd,M,L,i-2,nMin(i),nMax(i),labelType);
    end

    %the Y list represnt the type of trees:
    % - 1 left asymmetric
    % - 2 symmetric
    % - 3 right asymmetric
    Ytr = [ ones(dd,1); 2*ones(dd,1); 3*ones(dd,1)];

    %% create the test set

    Xtest = cell(Ntest,1);
    dd = Ntest/3;
    for i=1:3
        Xtest((i-1)*dd+1:i*dd) = treeGen(dd,M,L,i-2,nMin(i),nMax(i),labelType);
    end

    %the Y list represnt the type of trees:
    % - 1 left asymmetric
    % - 2 symmetric
    % - 3 right asymmetric
    Ytest = [ ones(dd,1); 2*ones(dd,1); 3*ones(dd,1)];
    
    save('synt.mat','L','M','Ntest','Ntr','Xtest','Xtr','Ytest','Ytr');
end