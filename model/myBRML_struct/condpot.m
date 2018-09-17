function outPot = condpot( pot,leftVars,rightVars)
%CONDPOT Condition the logpot pot s.t. it represents 
%   P(leftVars | rightVars)

    if(nargin < 2)
        leftVars = pot.variables;
    end
    
    if(nargin < 3)
        rightVars = [];
    end
    
    leftVars = sort(leftVars);
    rightVars = sort(rightVars);
    
    numVars = fastUnion(leftVars,rightVars);
    denVars = rightVars;
    
    numPot = sumpot(pot,numVars,0);
    denPot = sumpot(pot,denVars,0);
    
    outPot = divpot(numPot,denPot);
end

