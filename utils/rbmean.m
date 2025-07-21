function varagout = rbmean(data,dim,k,fit)
if ~exist('dim','var'),dim = 1;end
if ~exist('k','var'),k   = 3;end
if ~exist('fit','var'),fit   = 0;end
if size(data,dim)<4
    varagout = mean(data,1);
else
    varagout = robustMean(data,dim,k,fit);
end