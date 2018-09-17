function [ diffA ] = fastSetdiff(A,B)
%FASTSETDIFF Compute A/B efficently. A and B must be sorted and with no
%repeated elements.
    i=1; j=1; 
    nA = length(A); nB = length(B);
    C = true(1,nA);
    while(i<=nA && j<=nB)
        if(A(i) == B(j))
            C(i) = false;
            i = i+1;
            j = j+1;
        elseif(A(i) < B(j))
            i = i+1;
        else
            j = j+1;
        end
    end
    diffA = A(C);
end

