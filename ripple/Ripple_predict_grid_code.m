clear,clc


restoredefaultpath
addpath('../utils/')
GC = load('../data/mat/grid_code.mat');
dtacc   = readtable('../data/table/dd_subject_performance.csv');
channel = readtable('../data/table/Channel_ROI.csv');
channel.numid = arrayfun(@(x) {num2str(x)},1:size(channel,1))';
 
dtripple = readtable('../data/table/dd_subject_rest_and_overall_ripple_rate.csv');


%% Test for selective six-fold modulation in entorhinal
ec_ind = channel.Entorhinal>0;
subject = channel.subject(ec_ind);
beta_ec = GC.Mean_BETAs(ec_ind ,:);
beta_mpfc = GC.Mean_BETAs(channel.DMNmPFC>0,:);

% remove outliers
beta_ec_noext = beta_ec;
beta_ec_noext(abs(zscore(beta_ec,1,1))>3)=nan;
beta_mpfc_noext = beta_mpfc;
beta_mpfc_noext(abs(zscore(beta_mpfc,1,1))>3)=nan;

[beta_ec_subj, g_subj] = groupsummary(beta_ec_noext,subject,@nanmean);
beta_ec_subj = table(beta_ec_subj,g_subj,'VariableNames',{'beta','subj'});

% Test if six-fold modulation is greater than zero
[~,ecP,  ~, stat] = ttest(beta_ec_noext)
[~,mpfcP,~, stat] = ttest(beta_mpfc_noext)

% collect data for surface visualization 
ec_sixfold = table(beta_ec(:,3),channel.numid(ec_ind)); 
ec_sixfold.Properties.VariableNames = {'beta','channelid'};

% Average six-fold modulation at all time points.
beta_ec_t = squeeze(mean(GC.Point_BETAs(ec_ind,GC.time>0,:),2));
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

[r_acc,p_acc]=corr(dtacc.testAfter,dtacc.compAll) 
%% visualization, main figure 4

H = figure(4);clf
plt = myFigure([],[],H);
w = 0.13;

ax = axes('Position',[0.05 0.7 w w]);
plt.DATA = [dtacc.testAfter,dtacc.compAll];
plt.Colors = plt.CM.gray;
plt.addScatter([],[])
plt.addShadeErrorForGLM
plt.addPstar(sprintf('r = %.2f',r_acc),p_acc)
xticks(0.7:0.1:1);
xlabel('Feature test'),ylabel('Compound test')


ax = axes('Position',[0.35 0.7 w w]);
dz = beta_ec_noext;
plt.DATA = dz;
plt.P  = ecP;
plt.X = [4 5 6 7 8];
plt.Colors = repmat(plt.CM.gray,5,1);
plt.Colors(3,:) = plt.CM.sixfold;
plt.addBar('PstarOffset',0.255)
ylabel(sprintf('Modulation coefficients'))
plt.ylim([-1 1])
title('Entorhinal','FontWeight','normal')
xlabel('Symmetry')


ax = axes('Position',[0.05 0.4 w*2.4 w]);
act = squeeze(GC.bin_Activation(ec_ind,:));
dz = zscore(act,1,'all');  
plt.DATA = dz;

plt.Colors = repmat([plt.CM.sixfold;plt.CM.gray],6,1);
plt.addBar('ErrorbarStyle','T')
plt.ylim([-1 1]*0.5)
title('Entorhinal','FontWeight','normal')
ylabel('Oscillatory power (z-score)')
xlabel('Direction relative to \phi (°)')
xticks([1 12])
xticklabels([0 330])
plt.setFont
 
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


%% vsualization, main figure 4, grid code on EC surface
 
CH = load('../data/mat/Contact_distribution_DMN.mat');
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
 
%% visualization, supplementary figure 3, sub-regions of DMN
H = figure(5);clf
plt.H = H;
roiShorts = {'DMNmPFC','DMNTPJ','DMNLOFC','DMNPCC','DMNTL','DMNDLPFC'};
roiNames = plt.f_DMNsubNames(roiShorts);
w = 0.15
ax = axes('Position',[0.45 0.7 w w]);
for iroi = 1:numel(roiShorts)
    roi_ind = channel.(roiShorts{iroi})>0;
    beta_roi = GC.Mean_BETAs(roi_ind ,:);
    beta_roi_noext = beta_roi;
    beta_roi_noext(abs(zscore(beta_roi,1,1))>3)=nan;
    [h,p,ci,stat] = ttest(beta_roi_noext);
    
    b = beta_roi_noext(:,3);
    plt.DATA = b;plt.X =iroi;
    if iroi ==1, plt.Colors=plt.CM.sixfold;
    else,plt.Colors = plt.CM.gray;end
    plt.P = p(3);
    hs =plt.addBar('pstaroffset',0.1);
    hold on
end
xticks(1:6)
xticklabels(roiNames)
ylabel('Hexadirectional modulation')
plt.ylim([0 0.5])
yticks([0 0.1 0.2 0.3 0.4 0.5])
plt.updatePstartLoc()
xlim([0 7])
%
ax = axes('Position',[0.78 0.7 w w]);
plt.DATA = GC.Mean_BETAs(channel.DMNmPFC>0,:);
plt.DATA(abs(zscore(plt.DATA,1,1))>3)=nan;
plt.X = [ 4 5 6 7 8];
plt.Colors = repmat(plt.CM.gray,5,1);
plt.Colors(3,:)  = plt.CM.sixfold;
[h,p,ci,stat]=ttest(plt.DATA);
plt.addBar('P',p,'Pstaroffset',0.2)
plt.ylim([-0.2 0.6])
ylabel('Modulation coefficients')
s = numel(unique(channel.subject(channel.DMNmPFC>0)));
title(sprintf('mPFC (n=%d, N=%d)',size(plt.DATA,1), s),'FontWeight','normal')
xlabel('Symmetry')
plt.setFont


%%

cpath = path;
restoredefaultpath
%%

addpath('/HDD/lbcode-dev-xzb/externalPackages/fieldtrip-20220714')
ft_defaults
addpath(genpath('/HDD/lbcode-dev-xzb/externalPackages/eeglab2022.1/functions/'))
addpath('../utils/')
time =GC.time/1000;
trng = [min(time),max(time)];
 
ylims = [-2 3; -0.5 1];


ROIs  =  {'Entorhinal' 'DMNmPFC'};
for ii = 1:2
    ROI = ROIs(ii);
    subplot(4,2,4+ii)

    braw = GC.Point_BETAs(channel.(ROI{1})>0,:,3);
    b = braw-mean(braw(:,time<0),2);
    se = std(b)/sqrt(size(b,1));
    m = mean(b,1);
    shadedErrorBar(time,m,se)
    xlim(trng)
    ylim(ylims(ii,:))
    hold on
    plot(trng,[0,0],'k--')
    ylm = get(gca,'YLim');
    plot([0,0],ylm,'k--')
    title(ROI{1})

    bs = b;
    bs = movmean(bs,10,2);
 
    erps = {bs'; bs'*0};
    M_group = {'groupstats', 'off'};
    M_cond = {'condstats', 'on', 'paired',{'on'}};
    M_state = [M_cond,M_group];

    [pcond, pgroup, pinter, statscond, statsgroup, statsinter] = std_stat(erps, ...
        M_state{:} , ...
        'fieldtripmcorrect', 'cluster', ...
        'fieldtripmethod', 'montecarlo', ...
        'mode', 'fieldtrip');
    adj_p = pcond{1};


    T_mask = adj_p<0.05;
    mt =time;   mt(~T_mask)=nan;
    plot(mt,ylm(2)*ones(size(time)),'-','LineWidth',5,'Color',[0.7 0 0])
    box off
end 

restoredefaultpath
addpath(cpath)

plt.savefig(sprintf('Figures/supp_grid_code_DMN_%s.pdf','1'))

%% visualization, supplementary fig 3, grid-code on DMN surface
roi_ind = channel.Y7_Default>0;
beta_DMN = GC.Mean_BETAs(roi_ind ,:);
subject = channel.subject(roi_ind);
 
DMN_sixfold = table(beta_DMN(:,3),channel.numid(roi_ind));
DMN_sixfold.Properties.VariableNames = {'beta','channelid'};

 
CH = load('../data/mat/Contact_distribution_DMN.mat');
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
