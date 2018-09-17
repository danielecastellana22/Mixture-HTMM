function outPot = setpot(pot,varsToSet,valuesToSet)
%SETPOT Set the value of variables in a logpot.
%   In particular:
%   - varsToSet is the list of variables that will be setted
%   - valueToSet is the value to be setted

    vars = pot.variables;
    newVarsIdx = true(1,length(vars));
    [varsToSet,idxSorted] = sort(varsToSet);
    valuesToSet = valuesToSet(idxSorted);
    idxCell = cell(1,length(vars));
    idxCell(:) = {':'};
    idx = 1;
    i = 1;
    while(i<=length(vars) && idx<=length(varsToSet))
        if(vars(i) == varsToSet(idx))
            idxCell(i) = {valuesToSet(idx)};
            newVarsIdx(i) = false;
            idx = idx +1;
        elseif(vars(i) > varsToSet(idx))
            idx = idx +1;
        else
            i=i+1;
        end
    end
    newTab = pot.table(idxCell{:});
    newVars = vars(newVarsIdx);
    newTab = squeeze(newTab);
    if(isscalar(newTab))
        newVars = [];
    elseif(isrow(newTab))
        newTab = newTab';
    end
    outPot.variables = newVars;
    outPot.table = newTab;
end

