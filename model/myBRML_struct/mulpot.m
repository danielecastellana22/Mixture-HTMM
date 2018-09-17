function outPot = mulpot(pot1,pot2)
%MULPOT Multiply the two logpot pot1 and pot2, i.e. outPot = pot1*pot2

    [newVars,sizeVars,posVar1,posVar2] = mergeVariables(pot1, pot2);

    if(length(newVars) <= 1)
        %there is ONE variable in each pot (and it is the same)
        newTab = pot1.table + pot2.table;
    else
        toReshape1 = ones(1,length(newVars));
        toReshape2 = ones(1,length(newVars));
        toRepmat1 = sizeVars;
        toRepmat2 = sizeVars;

        toReshape1(posVar1) = sizeVars(posVar1); 
        toReshape2(posVar2) = sizeVars(posVar2);

        toRepmat1(posVar1) = ones(length(pot1.variables),1);
        toRepmat2(posVar2) = ones(length(pot2.variables),1);

        newTab1 = reshape(pot1.table,toReshape1);
        newTab1 = repmat(newTab1,toRepmat1);
        newTab2 = reshape(pot2.table,toReshape2);
        newTab2 = repmat(newTab2,toRepmat2);

        newTab = newTab1 + newTab2;
    end
    if(isscalar(newTab))
        newVars = [];
    elseif(isrow(newTab))
        newTab = newTab';
    end
    outPot.variables = newVars;
    outPot.table = newTab;
end

