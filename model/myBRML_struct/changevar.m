function outPot = changevar(pot,oldVar,newVar)
%CHANGEVAR Create new logpot replacing the variables
%   Detailed explanation goes here

    outPot.variables = newVar;
    outPot.table = pot.table;
end

