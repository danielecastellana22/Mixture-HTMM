function outSamples = samplepot(pot,nSamples)
%SAMPLEPOT Sample the logpot pot nSamples times.
    
    %remove the log
    tab = exp(pot.table);
    nStates = size(tab);
    nVars = length(pot.variables);
    outSamples = zeros(nVars,nSamples);
    idx=cell(nVars,1);
    for samp=1:nSamples
        randindex = randgen(tab(:));
        [idx{:}]= ind2sub(nStates,randindex);
        outSamples(:,samp) = cell2mat(idx);
    end
end

