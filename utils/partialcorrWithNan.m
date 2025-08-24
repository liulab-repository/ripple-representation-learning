function [R,P,N,D]=partialcorrWithNan(x,y,z,rmExtrem,coreFun,varargin)
% compute correlation with nan data, nan data will be exclud by pairs
if ~exist('y','var'),y = x;end
if ~exist('rmExtrem','var'),rmExtrem = false;end
if ~exist('coreFun','var'),coreFun = 1;end
R = nan(size(x,2),size(y,2));
P = R;
N = R;
if size(z,2)<size(x,2), z = repmat({z},1,size(x,2));end
D={};
for i = 1:size(x,2)
    for j = 1:size(y,2)
        if iscell(z)
            d = [x(:,i), y(:,j), z{i}];
        else
            d = [x(:,i), y(:,j), z];
        end
        nanind = isnan(d);
        d(sum(nanind,2)>0,:)=[];
        if rmExtrem
            [~,rmind] = rmoutliers(d(:,1:2),'median');
            d(rmind,:)=[];
        end
        if coreFun==1
            [r,p] = partialcorr(d(:,1),d(:,2),d(:,3:end), varargin{:});
            pname = 'r';
        elseif coreFun==2
     
        elseif coreFun ==3
            
        end
        n = size(d,1);
        R(i,j) = r;
        P(i,j) = p;
        N(i,j) = n;
        D{i,j} = d;
        if nargout<1
            [flag]=f_pValue2flag(p);
            if p>=0.05,flag='';end
            fprintf('n=% 5d  ||  %s = %7.3f, p = %.3f  %s\n',n,pname,r,p,flag)
        end
        
    end
end


function [flag]=f_pValue2flag(p)
threshs=[0.05 0.01 0.001 ];
flags  = {'*' '**' '***' };

flag = cell(size(p));
for i = 1:numel(p)
    ip = p(i);
    ind = find(ip<threshs,1,'last');
    if isempty(ind)
        flag{i}= 'n.s.';
    else
%         flag{i} = ['$' flags{ind} '$'];
        flag{i} = [flags{ind}];
    end
end

if numel(p)==1
    flag = flag{1};
end
