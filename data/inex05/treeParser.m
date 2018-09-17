function [t,L,M] = treeParser(s)
   %INEX05TreeParser parse a string which represents a INEX05 tree to data
   %structure
    
    currIdx = 1;
    nextId = 1;
    adj = [];
    v = [];
    nf = 0;
    L=0;
    parser();
    n = nextId -1;
    newAdj = zeros(n,n);
    szAdj = size(adj);
    M=max(v);
    newAdj(1:szAdj(1),1:szAdj(2)) = adj;
    t = treeRep(newAdj,n,nf,v,L);
    t = t.format();

    function parser()
        
        %read label
        idx=currIdx;
        while(s(idx) ~= '(')
            idx = idx+1;
        end
        
        label = str2double(s(currIdx:idx-1));
        v(nextId) = label;
        paId = nextId;
        nextId = nextId +1;
        idx = idx + 1;
        currIdx = idx;
        
        if(s(currIdx) == '$')
            nf = nf +1;
            currIdx = currIdx + 2;
        else
            countL = 0;
            while(s(currIdx)~=')')
                nextChId = nextId;
                countL = countL + 1;
                parser();
                adj(paId,nextChId) = countL;
            end
            
            L = max(L,countL);
            currIdx = currIdx + 1;
        end
    end
end
