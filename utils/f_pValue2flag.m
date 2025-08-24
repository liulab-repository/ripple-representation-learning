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
