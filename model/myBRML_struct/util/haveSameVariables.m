function ris = haveSameVariables(pot1,pot2)
%HAVESAMEVARIABLES Return true if pot1 and pot2 have the same variables.
    ris = all(pot1.variables == pot2.variables);
end

