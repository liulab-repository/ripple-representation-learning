function [b,DEV,STATS] = f_pinv_GLM(X,y)
% for mex compability 
% author, zhibing xiao
int = ones(numel(y),1);
if ~isequal(X(:,1),int)
    X = [int, X];
end

b = pinv(X)*y;

e = y-X*b;
SSE = sum(e.^2);

n = size(X, 1);  %  
p = size(X, 2);  %  
X_corr = X' * X;  %  
cov_mat = SSE / (n - p) * pinv(X_corr);  % 
se = sqrt(diag(cov_mat));  %  
 
t = b ./ se;
df = n - p;
p = 2 * (1 - tcdf(abs(t), df));

STATS.se= se;
STATS.t = t;
STATS.p = p;

DEV = nan;
