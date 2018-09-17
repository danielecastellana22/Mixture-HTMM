function [ C ] = fastUnion( A,B )
%FASTUNION Compute C=A U B efficently. A and B must be sorted and with no
    %repeated elements.
    nA = length(A);
    nB = length(B);
    nVars = length(A) + length(B);
    C = zeros(nVars,1);
    idxC = 1;
    idxA = 1;
    idxB = 1;
    while(idxA <= nA && idxB <= nB)
        if(A(idxA) < B(idxB))
            C(idxC) = A(idxA);
            idxA = idxA + 1;
        elseif(B(idxB) < A(idxA))
            C(idxC) = B(idxB);
            idxB = idxB + 1;
        else
            C(idxC) = A(idxA);
            idxA = idxA + 1;
            idxB = idxB + 1;
        end
        idxC = idxC+1;
    end

    %copy element still in p1
    while(idxA <= length(A))
        C(idxC) = A(idxA);
        idxA = idxA + 1;
        idxC = idxC+1;
    end

    %copy element still in p2
    while(idxB <= length(B))
        C(idxC) = B(idxB);
        idxB = idxB + 1;
        idxC = idxC+1;
    end

    C = C(1:idxC-1);
end

