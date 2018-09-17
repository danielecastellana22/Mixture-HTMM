function [ treeList ] = treeGen(n,M,L,treeType,nMin,nMax,labelType,distrToLabel)
%TREEGEN Generate tree with specific structural characteristics
%   Detailed explanation goes here
    %n:             number of tree
    %[1,M]:         range for the label
    %L:             number max of children
    %nMin:          min number of nodes in the tree
    %nMax:          max number of nodes in the tree
    %type:          -1 asymetric left tree
    %               0 symmetric tree
    %               1 asymetric right tree
    %labelType:     specify how to generate the label
    %               0 random generation
    %               1 the number of child +1 to avoid 0
    %               2 fixed, always the same label
    %               3 the label are generated accoridng to a distribution
    %               random if nothing is specified
    %distrToLabel:  distribution to generate label. It is used only id
    %               labelType = 3
    
    if(nargin < 7)
        labelType = 0;
    end
    
    
    treeList = cell(1,n);
    
    a = log10(1/10);
    b = log10(9/10);
    
    if(treeType == 0)
        prCh = (4/5)*ones(1,L);
    end
    
    if(treeType==-1)
        prCh = logspace(b,a,L);
    end
    
    if(treeType==1)
        prCh = logspace(a,b,L);
    end
    
    
    for i=1:n
        
        %choose the number od nodes in the tree
        nNodes = round((nMax-nMin)*rand(1) + nMin);
        isBigEnough = false;
        while(~isBigEnough)
            %initialise
            adjMat = zeros(nNodes,nNodes);
            newId = 2;
            idPa = 1;
            %number of leaf
            nf=0;
            v = zeros(nNodes);
            v(1) = 1;
            endId = false;

            while(~endId && idPa < newId)
                %ensure property is true near root
                if(idPa <=2)
                    isLeaf = false;
                    if(treeType==0)
                        for l=1:L
                            adjMat(idPa,newId)=l;
                            
                            if(labelType == 3)
                                v(newId) = getLabel(v(idPa),l);
                            end
                            newId = newId+1;
                            if(newId>nNodes)
                                endId = true;
                                break
                            end
                        end
                    else
                        [~,l] = max(prCh);
                        adjMat(idPa,newId)=l;
                        
                        if(labelType == 3)
                            v(newId) = getLabel(v(idPa),l);
                        end
                        
                        newId = newId+1;
                        if(newId>nNodes)
                            endId = true;
                        end
                    end
                else
                    %genrate the children
                    isLeaf = true;
                    for l=1:L
                        p = rand(1);
                        if(p<prCh(l))
                            isLeaf = false;
                            %add the l-th child of p
                            adjMat(idPa,newId)=l;
                            
                            if(labelType == 3)
                                v(newId) = getLabel(v(idPa),l);
                            end
                            
                            newId = newId+1;
                            if(newId>nNodes)
                                endId = true;
                                break
                            end
                        end
                    end
                end
                if(isLeaf)
                    %no generation
                    nf = nf+1;
                end
                idPa = idPa + 1;
            end
            
            n = newId - 1;
            if(n==nNodes)
                isBigEnough = true;
            end
        end    
        n = newId - 1;
        adjMat = adjMat(1:n,1:n);      
        v = v(1:n);
        nf = nf + (newId-idPa);
        
        if(labelType == 0)
            v=ceil(rand(1,n)*M);
        end
        
        if(labelType == 1)
            v=sum(adjMat~=0,2)+1;
        end
        
        if(labelType == 2)
            v=ones(1,M);
        end
        
        tr = treeRep(adjMat,n,nf,v,L);
        treeList{i} = tr.format();
        
    end
    
    function lab = getLabel(vPa,posCh)
        distr = distrToLabel(:,vPa,posCh);
        cumDistr = cumsum(distr);
        
        tc = rand(1);
        lab=1;
        while(tc >= cumDistr(lab))
            lab = lab+1;
        end
    end
end

