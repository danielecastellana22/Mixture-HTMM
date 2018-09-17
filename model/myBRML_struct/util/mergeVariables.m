function [mergedVars,varsSize,posVar1,posVar2] = mergeVariables(pot1,pot2)
%PLUSPOT Merge the variables contained in pot1 and pot2.
%   In particular:
%   - mergedVars contains the pot1.variables U pot2.variables
%   - varsSize contains the dimension for each variable
%   - posVar1 indacates positions occupied by pot1.variables in mergeVars
%   - posVar2 indacates positions occupied by pot2.variables in mergeVars

        var1 = pot1.variables;
        szTab1 = size(pot1.table);
        var2 = pot2.variables;
        szTab2 = size(pot2.table);
        nVars = length(var1) + length(var2);
        posVar1 = zeros(1,length(var1));
        posVar2 = zeros(1,length(var2));
        mergedVars = zeros(1,nVars);
        varsSize = zeros(1,nVars);
        idxW = 1;
        idxP1 = 1;
        idxP2 = 1;
        while(idxP1 <= length(var1) && idxP2 <= length(var2))
            if(var1(idxP1) < var2(idxP2))
                mergedVars(idxW) = var1(idxP1);
                varsSize(idxW) = szTab1(idxP1);
                posVar1(idxP1) = idxW;
                idxP1 = idxP1 + 1;
            else

                if(var2(idxP2) < var1(idxP1))
                    mergedVars(idxW) = var2(idxP2);
                    varsSize(idxW) = szTab2(idxP2);
                    posVar2(idxP2) = idxW;
                    idxP2 = idxP2 + 1;
                else
                    if(szTab2(idxP2) ~= szTab1(idxP1))
                        error('Size of vriable mismatch.');
                    end
                    mergedVars(idxW) = var1(idxP1);
                    varsSize(idxW) = szTab1(idxP1);
                    posVar1(idxP1) = idxW;
                    posVar2(idxP2) = idxW;
                    idxP1 = idxP1 + 1;
                    idxP2 = idxP2 + 1;
                end
            end
            idxW = idxW+1;
        end

        %copy element still in p1
        while(idxP1 <= length(var1))
            mergedVars(idxW) = var1(idxP1);
            varsSize(idxW) = szTab1(idxP1);
            posVar1(idxP1) = idxW;
            idxP1 = idxP1 + 1;
            idxW = idxW+1;
        end

        %copy element still in p2
        while(idxP2 <= length(var2))
            mergedVars(idxW) = var2(idxP2);
            varsSize(idxW) = szTab2(idxP2);
            posVar2(idxP2) = idxW;
            idxP2 = idxP2 + 1;
            idxW = idxW+1;
        end

        mergedVars = mergedVars(1:idxW-1);
        varsSize = varsSize(1:idxW-1);        
end

