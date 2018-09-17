function outPot = pluspot(pot1,pot2)
%PLUSPOT Return the sum of the logpot, i.e. outPot = pot1 + pot2

    if(haveSameVariables(pot1,pot2))
        %find max for each variable (i.e. dimension)
        maxP1 = max(pot1.table);
        maxP2 = max(pot2.table);
        logprefactor = max(maxP1(:),maxP2(:));
        infPos = isinf(logprefactor);
        logprefactor(infPos) = 10^8 .* sign(logprefactor(infPos));
        szTab = size(pot1.table);
        nVars = ndims(pot1.table);
        logprefactor = reshape(logprefactor,[1 szTab(2:end)]);
        logprefactorREP = repmat(logprefactor,[szTab(1) ones(1,nVars-1)]);
        newTab = log(exp(pot1.table-logprefactorREP) + exp(pot2.table-logprefactorREP))+logprefactorREP;

        outPot.variables = pot1.variables;
        outPot.table = newTab;
    else
        error('The sum must be on potentials wiht the same variables');
    end
end

