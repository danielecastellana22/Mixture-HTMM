function generateInex05Dataset()

    BRML_dir = '../../model/';
    p = genpath(BRML_dir);
    addpath(p); 
    
    filenameTest = 'inex05_raw_data/inex05.test.elastic.tree';
    filenameTrain = 'inex05_raw_data/inex05.train.elastic.tree';


    [Xtr,Ytr,Ltr,Mtr] = datasetParser(filenameTrain);
    [Xtest,Ytest,Ltest,Mtest] = datasetParser(filenameTest);

    M = max(Mtr,Mtest);
    L = max(Ltr,Ltest);
    Ntr = length(Xtr);
    Ntest = length(Xtest);
    
    save('inex05.mat','L','M','Ntest','Ntr','Xtest','Xtr','Ytest','Ytr');
end