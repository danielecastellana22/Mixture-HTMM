function [X,Y,L,M] = datasetParser(fileName)
%PARSEDATASET Summary of this function goes here
%   Detailed explanation goes here
    [Y,X_str] = fileParser(fileName);
    X = cell(length(X_str),1);
    L=0;
    M=0;
    for i=1:length(X_str)
        [X{i},l,m] = treeParser(X_str{i});
        L = max(L,l);
        M = max(M,m);
    end
end

