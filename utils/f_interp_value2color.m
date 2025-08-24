function [mappedColors, vals_col] = f_interp_value2color(values,cmap,cbar_lim)
% author, ZHibing Xiao

if ~exist('cbar_lim','var')
    cbar_lim = [min(values) max(values)];
end

%% compute colors:
n_colors=size(cmap,1);
cbar_min=cbar_lim(1);
cbar_max=cbar_lim(2);
cmapTH = cmap;
if length(cbar_lim)>2
    thr = cbar_lim(3);
    thrColorInd1 = floor(interp1([cbar_min cbar_max],[1 n_colors],max(-thr,cbar_min),'linear',n_colors));
    thrColorInd2 = floor(interp1([cbar_min cbar_max],[1 n_colors],min(thr,cbar_max),'linear',n_colors));
    tmp = ones(size(cmap(thrColorInd1:thrColorInd2,:)));
    cmapTH(thrColorInd1:thrColorInd2,:)=[tmp(:,1)*belowThrColor(1),tmp(:,2)*belowThrColor(2),tmp(:,3)*belowThrColor(3),]; % below thr "null" color
    clear tmp
end

% compute electrode colors:
vals_col=ceil(interp1([cbar_min cbar_max],[1 n_colors],values,'linear',n_colors));
mappedColors = cmapTH(vals_col,:);


% vals_col=ceil(interp1([cbar_min cbar_max],[1 n_colors],values,'linear',n_colors));
