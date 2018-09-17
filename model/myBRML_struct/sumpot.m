function outPot = sumpot(pot,vars,sumOver)
%SUMPOT Marginalise the logpot. 
%   In particular:
%   - vars is a list of variables
%   - sumOver if 1, we marginalise over vars; otherwise, we marginalise on
%   all the others variables
    if(nargin<2)
        vars = pot.variables;
    else
        % sort the list given by the user
        vars = sort(vars);
    end

    if(nargin<3)
        sumOver = 1;
    end
    
    if(sumOver)
        varsToSumOver = vars;
    else
        varsToSumOver = fastSetdiff(pot.variables,vars);
    end
    % varsToSumOver is always sorted
    
    % compute the marginalisation
    if(isempty(varsToSumOver))
        outPot = pot;
    else
        vars = pot.variables;
        newVarsIdx = true(1,length(vars));
        newTab = pot.table;
        idx = 1;
        i = 1;
        while(i<=length(vars) && idx<=length(varsToSumOver))
            if(vars(i) == varsToSumOver(idx))
                %avoid underflow
                logprefactor=max(newTab,[],i);
                %remove infinity
                infPos = isinf(logprefactor);
                logprefactor(infPos) = 10^8 .* sign(logprefactor(infPos));
                %reshape the logprefactor
                if(isscalar(logprefactor))
                    szRep = size(newTab);
                else
                    szTab = size(newTab);
                    szPre = ones(1,ndims(newTab));
                    szPre(1:ndims(logprefactor)) = size(logprefactor);
                    szRep = szTab - szPre + 1;
                end
                logprefactorREP = repmat(logprefactor,szRep);

                newTab=log(sum(exp(newTab-logprefactorREP),i))+logprefactor;

                newVarsIdx(i) = false;
                idx = idx + 1;
            elseif(vars(i) > varsToSumOver(idx))
                    idx = idx + 1;
            else
                i = i+1;    
            end
        end
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
end

