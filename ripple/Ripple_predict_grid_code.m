clear,clc
restoredefaultpath
addpath('../utils/')
load('../grid-like code/Results/grid_code.mat')
dtacc = readtable('../data/dd_subject_performance.csv');
channel = readtable('../data/Channel_ROI.csv');
channel.numid = arrayfun(@(x) {num2str(x)},1:size(channel,1))';

dtripple = readtable('../data/dd_subject_rest_and_overall_ripple_rate.csv');
%% Test for selective six-fold modulation in entorhinal
ec_ind = channel.Entorhinal>0;
subject = channel.subject(ec_ind);
beta_ec = Mean_BETAs(ec_ind ,[4 5 6 7 8]);
beta_ec_noext = beta_ec;
beta_mpfc = Mean_BETAs(channel.DMNmPFC>0,[4 5 6 7 8]);
beta_mpfc_noext = beta_mpfc;

% remove outliers
beta_ec_noext(abs(zscore(beta_ec,1,1))>3)=nan;
beta_mpfc_noext(abs(zscore(beta_mpfc,1,1))>3)=nan;

% Test if six-fold modulation is greater than zero
[~,ecP,  ~, stat] = ttest(beta_ec_noext)
[~,mpfcP,~, stat] = ttest(beta_mpfc_noext)

% collect data for surface visualization 
ec_sixfold = table(beta_ec(:,3),channel.numid(ec_ind)); 
ec_sixfold.Properties.VariableNames = {'beta','channelid'};

% Average six-fold modulation at all time points.
beta_ec_t = squeeze(mean(Point_BETAs(ec_ind,time>0,[4 5 6 7 8]),2));
[beta_ec_s,gname]  = groupsummary(beta_ec_t,subject,@nanmean);
gc_ec_s = table(beta_ec_s(:,3),gname,'VariableNames',{'beta_ec','subject'});

% combine grid-code with subject performance
dtGCandACC = innerjoin(dtacc,gc_ec_s,"Keys",'subject');
 
% combine grid-code with subject ripple rate
dtGCandRipple = innerjoin(dtripple, gc_ec_s,'Keys','subject');
dtGCandRipple = innerjoin(dtGCandRipple,dtacc(:,{'subject','testAfter'}),'Keys','subject');

%% correlation between grid-like code and ripple/performance

[r_behav,p_behav]   = partialcorr(dtGCandACC.testAfterInfer, dtGCandACC.beta_ec, dtGCandACC.testAfterMemory)
[r_ripple,p_ripple] = partialcorr(dtGCandRipple.RR_Rest,dtGCandRipple.beta_ec,dtGCandRipple.RR_Overall)

%% visualization, main figure 4

H = figure(4);clf
plt = myFigure([],[],H);
w = 0.15;

ax = axes('Position',[0.05 0.7 w w]);
plt.DATA = [dtacc.testAfter,dtacc.compAll];
plt.Colors = plt.CM.gray;
[r,p] = corr(plt.DATA(:,1),plt.DATA(:,2));
plt.addScatter([],[])
plt.addShadeErrorForGLM
plt.addPstar(sprintf('r = %.2f',r),p)
xticks(0.7:0.1:1);
xlabel('Feature test'),ylabel('Compound test')


ax = axes('Position',[0.35 0.7 w w]);
plt.DATA=zscore(beta_ec,1,'all');  
plt.P  = ecP;
plt.X = [4 5 6 7 8];
plt.Colors = repmat(plt.CM.gray,5,1);
plt.Colors(3,:) = plt.CM.sixfold;
plt.addBarWithSpread('PstarOffset',0.1)
ylabel(sprintf('Modulation \\beta (z-scored)'))
plt.ylim([-1 1]*0.5)
title('Entorhinal','FontWeight','normal')


ax = axes('Position',[0.05 0.4 w*1.8 w]);
act = squeeze(Mean_SixFold.binActivation(ec_ind,:));
plt.DATA = zscore(act,1,'all');  
plt.Colors = repmat([plt.CM.sixfold;plt.CM.gray],6,1);
plt.addBarWithSpread()
plt.ylim([-1 1])
title('Entorhinal','FontWeight','normal')
ylabel('Oscillatory power (z-score)')


plt.DATA = [dtGCandACC.testAfterInfer,dtGCandACC.beta_ec,dtGCandACC.testAfterMemory];
ax = axes('Position',[0.5   0.4 w  w]);
plt.Colors = plt.CM.sixfold;
plt.addPartialRegressionScatter;
plt.addPstar(sprintf('r = %.3f',r_behav),p_behav)
ylim([-2.5 2.5]),yticks(-2.5:1:2.5)
xlabel({'Inference' '(adjusted)'})
ylabel({'Hexadirectional modulation' '(adjusted)'})
title('Entorhinal','FontWeight','normal')

plt.DATA = [dtGCandRipple.RR_Rest,dtGCandRipple.beta_ec,dtGCandRipple.RR_Overall];
ax = axes('Position',[0.79  0.4 w  w]);
plt.Colors = plt.CM.sixfold;
plt.addPartialRegressionScatter;
ylim([-2.5 2.5]),yticks(-2.5:1:2.5)
plt.addPstar(sprintf('r = %.3f',r_ripple),p_ripple)
xlabel({'Ripple rate during rest' '(adjusted)'})
ylabel({'Hexadirectional modulation' '(adjusted)'})
title('Entorhinal','FontWeight','normal')
plt.setFont

plt.savefig('Figures/ripple_grid_code.pdf')
 
%% visualization, supplementary figure 3, sub-regions of DMN
H = figure(5);clf
plt.H = H;
roiShorts = {'DMNmPFC','DMNTPJ','DMNLOFC','DMNPCC','DMNTL','DMNDLPFC'};
roiNames = plt.f_DMNsubNames(roiShorts);

ax = axes('Position',[0.45 0.7 w w]);
for iroi = 1:numel(roiShorts)
    roi_ind = channel.(roiShorts{iroi})>0;
    beta_roi = Mean_BETAs(roi_ind ,[4 5 6 7 8]);
    beta_roi_noext = beta_roi;
    beta_roi_noext(abs(zscore(beta_roi,1,1))>3)=nan;
    [h,p,ci,stat] = ttest(beta_roi_noext);
    
    b = beta_roi_noext(:,3);
    plt.DATA = b;plt.X =iroi;
    if iroi ==1, plt.Colors=plt.CM.sixfold;
    else,plt.Colors = plt.CM.gray;end
    plt.P = p(3);
    hs =plt.addSpreadWithBar();
    hold on
end
xticks(1:6)
xticklabels(roiNames)
ylabel('Hexadirectional modulation')
plt.ylim([-1 1.5])
plt.updatePstartLoc()
xlim([0 7])
%
ax = axes('Position',[0.78 0.7 w w]);
plt.DATA = Mean_BETAs(channel.DMNmPFC>0,[4 5 6 7 8]);
plt.DATA(abs(zscore(plt.DATA,1,1))>3)=nan;
plt.X = [ 4 5 6 7 8];
plt.Colors = repmat(plt.CM.gray,5,1);
plt.Colors(3,:)  = plt.CM.sixfold;
[h,p,ci,stat]=ttest(plt.DATA);
plt.addSpreadWithBar('P',p,'Pstaroffset',1)
plt.ylim([-1 1.5])
ylabel('Modulation coefficients')
s = numel(unique(channel.subject(channel.DMNmPFC>0)));
title(sprintf('mPFC (n=%d, N=%d)',size(plt.DATA,1), s),'FontWeight','normal')
xlabel('Symmetry')
plt.setFont
plt.savefig('Figures/supp_grid_code_DMN.pdf')

%% vsualization, main figure 4, grid code on EC surface
 
CH = load('../data/Contact_distribution_DMN.mat');
axs = findobj(CH.H,'type','axes');
delete(axs([1 3]))
axs = axs([4 2]);

elecH  = findobj(CH.H,'tag','Electrode');
elecID = get(elecH ,'UserData'); 

[~,ia,ib] = intersect(ec_sixfold.channelid,elecID,'stable');
bd = ec_sixfold.beta(ia);
bH = elecH(ib);
bids = elecID(ib);

mycmap =  brewermap(32,'RdYlBu');
mycmap = flip(mycmap);

lm =[-2 2];
bd(bd<-2)=-2;
bd(bd>2)=2;
cols=f_interp_value2color(bd,mycmap(:,:) );

for ic = 1:numel(bH)
    id = bids{ic};
    ind = ismember(elecID,id);    
    set(elecH(ind),'MarkerFaceColor',cols(ic,:) , ...
        'MarkerEdgeColor',cols(ic,:)*1,"MarkerSize",2, ...
        'DisplayName',num2str(bd(ic)))
    
end
disp done!
hb = colorbar("east",'Position',[0.93 0.3 0.02 0.4]);
hb.Ticks=lm;hb.Label.String = 'Hexadirectional modulation';
colormap(mycmap)
clim(lm)
 
[rm,ia]=setdiff(elecID,ec_sixfold.channelid);
delete(elecH(ismember(elecID,rm)))

%== pull out electrods for visualization
for ia = 1:2
    allElec = findobj(axs(ia),'tag','Electrode');
    elecValue = get(allElec ,'DisplayName');
    elecValue = cellfun(@str2num,elecValue );
    nonGrayElec = allElec;
    [~,ord]=sort(elecValue,1,'ascend');
    nonGrayElec = nonGrayElec(ord);
    get(nonGrayElec ,'DisplayName');

    if ~isempty(nonGrayElec)
        X = [nonGrayElec.XData];
        Y = [nonGrayElec.YData];
        Z = [nonGrayElec.ZData];
        S = [nonGrayElec.MarkerSize];
        pullingExtent = 5;
        maxEnlargingExtent = 1;
        cmap = mycmap;
        ncolors = size(cmap,1)-1;
        pulledCrd = [];
        for i = 1:length(nonGrayElec)
            if ia==1,ifc = i; else,ifc = -i;end
            pulledCrd=[X(i)+ifc,Y(i),Z(i)];
            set(nonGrayElec(i),'XData',pulledCrd(1),'YData',pulledCrd(2),'ZData',pulledCrd(3));%,'markersize',S(i)*(1+enlarging_factor*maxEnlargingExtent));
        end

    end
end
exportgraphics(CH.H,'Figures/grid_code_distribute_EC.png','Resolution',600)
 
%% visualization, supplementary fig 3, grid-code on DMN surface
roi_ind = channel.Y7_Default>0;
beta_DMN = Mean_BETAs(roi_ind ,[4 5 6 7 8]);
subject = channel.subject(roi_ind);
 
DMN_sixfold = table(beta_DMN(:,3),channel.numid(roi_ind));
DMN_sixfold.Properties.VariableNames = {'beta','channelid'};

 
CH = load('../data/Contact_distribution_DMN.mat');
axs = findobj(CH.H,'type','axes');
axs = axs([4 2 3 1]);


elecH  = findobj(CH.H,'tag','Electrode');
elecID = get(elecH ,'UserData'); 

[~,ia,ib] = intersect(DMN_sixfold.channelid,elecID,'stable');
bd = DMN_sixfold.beta(ia);
bH = elecH(ib);
bids = elecID(ib);

mycmap =  brewermap(32,'RdYlBu');
mycmap = flip(mycmap);

lm =[-2 2];
bd(bd<-2)=-2;
bd(bd>2)=2;
cols=f_interp_value2color(bd,mycmap(:,:) );

for ic = 1:numel(bH)
    id = bids{ic};
    ind = ismember(elecID,id);    
    set(elecH(ind),'MarkerFaceColor',cols(ic,:) , ...
        'MarkerEdgeColor',cols(ic,:)*1,"MarkerSize",2, ...
        'DisplayName',num2str(bd(ic)))
    
end
disp done!
 
hb = colorbar("east",'Position',[0.93 0.3 0.02 0.4]);
hb.Ticks=lm;hb.Label.String = 'Hexadirectional modulation';
colormap(mycmap)
clim(lm)

[rm,ia]=setdiff(elecID,DMN_sixfold.channelid);
delete(elecH(ismember(elecID,rm)))
 
exportgraphics(CH.H,'Figures/grid_code_distribute_DMN.png','Resolution',600)
