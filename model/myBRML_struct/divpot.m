function outPot = divpot(pot1,pot2)
%DIVPOT Divide the two logpot pot1 and pot2, i.e. outPot = pot1/pot2

    %avoid the inf in the denominator (i.e. zero)
    newDenPotTab = pot2.table;
    infPos = isinf(newDenPotTab);
    newDenPotTab(infPos) = 10^8 .* sign(newDenPotTab(infPos));
    pot2.table = -newDenPotTab;
    outPot = mulpot(pot1,pot2);
end

