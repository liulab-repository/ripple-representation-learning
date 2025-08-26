function [R,P,N,D]=corrWithNan(x,y,varargin)
% compute correlation with nan data, nan data will be exclud by pairs
if ~exist('y','var')
    y = x;
end
R = nan(size(x,2),size(y,2));
P = R;
N = R;
D = R;
for i = 1:size(x,2)
    for j = 1:size(y,2)
        d = [x(:,i), y(:,j)];
        nanind = isnan(d);
        d(sum(nanind,2)>0,:)=[];
        [r,p] = corr(d(:,1),d(:,2),varargin{:});
        n = size(d,1);
        R(i,j) = r;
        P(i,j) = p;
        N(i,j) = n;

        if nargout<1
            fprintf('n=% 5d  ||  r = %.3f, p = %.3f\n',n,r,p)
        end        
    end
end