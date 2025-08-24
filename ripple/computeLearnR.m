function LearnR = computeLearnR(trial,rippleRate)

if numel(trial)<10 || mean(rippleRate>0)<0.1
    % Discard calculations when there are few observations (i.e., limited correct trials) 
    % or sparse ripple events (i.e., no ripples detected in most trials) 
    % to avoid unstable results.
    LearnR = nan;
    return

else  % calculate correlation between trial number and ripple rate
        r = corr(trial,rippleRate,"Type","Spearman");
end


% Fishers' Z transform, for further statistics
LearnR = 0.5 * log((1 + r) ./ (1 - r));


end
 