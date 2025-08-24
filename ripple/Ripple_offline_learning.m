clear,clc
restoredefaultpath
addpath('../utils/')
mkdir('Results')
%% preparing data
% read data
dtSubjACC  = readtable('../data/table/dd_subject_performance.csv');
dtRipple = readtable("../data/table/dd_ripple_during_learning_raw.csv");
  
% # exclude extreme values
threshold_prc  = 99; 
dtRipple.exclude = false(size(dtRipple,1),1);
for phase  = {'stimulus' 'feedback' 'interval'}
    ind = ismember(dtRipple.phase,phase);
    dd = dtRipple.ripple(ind);
    thresh = prctile(dd(dd>0),threshold_prc);
    dtRipple.exclude(ind) = dtRipple.ripple(ind)>thresh;
end

dtRipple = dtRipple(~dtRipple.exclude,:); %
dtRipple2 = dtRipple(~ismember(dtRipple.phase,'stimulus'),:);

writetable(dtRipple,'Results/dt_ripple_stimulus_feedback_interval.csv')
writetable(dtRipple2,'Results/dt_ripple_feedback_interval.csv')

dtRippleM = unstack(dtRipple,'ripple','phase');  
dtRippleM.ACC = ismember(dtRippleM.ACC,'correct'); 
writetable(dtRippleM,'Results/dt_trial_level_Ripple.csv')

%% run models
% model 1, Compare ripple rates between feedback and interval
Rfile = fullfile(pwd,"R_scripts/Learning_RippleRate_ITIvsFeedback.R");
RUN_R_Script('',Rfile)

% model 2, Interaction between trials and ACC on ripple rate during each phase.
% This analysis was based on data from the first learning session to avoid 
% confounding effects from relearning or testing.
Rfile = fullfile(pwd,"R_scripts/Learning_RippleRate_ACCtrial.R");
RUN_R_Script('',Rfile)

% model 3, Compare the ripple rates in correct and incorrect response trials within each phase.
Rfile = fullfile(pwd,"R_scripts/Learning_RippleRate_ACCphases.R");
RUN_R_Script('',Rfile)


warning off
lmeResults = table;
n = 1;
for phase = {'stimulus' 'feedback' 'interval'}
    rp1 = readtable(sprintf(['Results/Learning_Ripple_Trial/interaction_' ...
        '%s_ACCtrial_slope.csv'],phase{:}));
    rp2 = readtable(sprintf(['Results/Learning_Ripple_Trial/mixed_model_' ...
        '%s.csv'],phase{:}));
    rp3 = readtable(sprintf(['Results/Learning_ripple_ACC_phase/mixed_effects_model_' ...
        '%s.csv'],phase{:}));    
    lmeResults.('phase')(n)=phase;
    lmeResults.interaction_p(n)   = rp2(end,"Pr___t__").Variables;
    lmeResults.simple_correct_p(n)= rp1.p(2);
    lmeResults.simple_wrong_p(n)  = rp1.p(1);
    lmeResults.M_ACC0vs1_p(n) = rp3.Pr__F_;
    n = n+1;
end
clear rp1 rp2 rp3 n
lmeResults(:,2:end).Variables = round(lmeResults(:,2:end).Variables,5);
warning on
lmeResults
 
%% Calculate learn-R in correct trials
% select correct trials
dt_correct = dtRipple(ismember(dtRipple.ACC,'correct'),:)  ; 
dt_correct.ACC = [];

% calculate the learn-R
T3     = groupsummary(dt_correct,{'channel' 'phase' }, ...
    @(x,y)computeLearnR(x,y),{"trial" "ripple" });

% remove outliers
phases = {'stimulus' 'feedback' 'interval'};
for i = 1:3
    indpha = find(ismember(T3.phase,phases{i}));
    [~,ind]=rmoutliers(T3.fun1_trial_ripple(indpha),"mean","ThresholdFactor",3);
    ind = find(ind);
    T3.fun1_trial_ripple(indpha(ind))=nan;
end
 
% reshape  data
T3.GroupCount=[];
T3m = unstack(T3,"fun1_trial_ripple",'phase');

% average to subject level
T3m.subject = regexprep(T3m.channel,'.*_','');  T3m.channel=[];
subjLearnR = grpstats(T3m,{'subject'},@nanmean);
subjLearnR.Properties.VariableNames = strrep(subjLearnR.Properties.VariableNames,'nanmean_','');
 
% combine with subject level performance
[~,ia,ib]= intersect(dtSubjACC.subject,subjLearnR.subject,'stable');
subjLearnR = subjLearnR(ib,:);
assert(isequal(ia,(1:height(dtSubjACC))'))
 
subjLearnR.eleInfer_perf = dtSubjACC.testBeforeInfer;
subjLearnR.eleMemory_perf= dtSubjACC.testBeforeMemory;
writetable(subjLearnR,'Results/dt_sublevel_LearnR_behav.csv')


%  compute correlation between Learn Increased Ripple and test performance
Rfile = fullfile(pwd,"R_scripts/LearnR_predict_performance.R");
RUN_R_Script('',Rfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off
learnR_results = readtable('Results/LearnR_predict_infer/learnR_predict_feature_performance.csv');
learnR_results.depend = repmat({'memory';'inference'},3,1);
warning on
  
learnR_results
 
  
%% visualization for main figure 2
plt  = myFigure(); clf

% plot 1, ripple rate: phase x acc
contrast_01 = readtable('Results/Learning_ripple_iti_vs_feedback/contrast_ACC_by_phase_Learning_ripple_iti_vs_feedback.csv');
contrast_fi = readtable('Results/Learning_ripple_iti_vs_feedback/contrast_phase_by_ACC_Learning_ripple_iti_vs_feedback.csv');
eemean  = readtable('Results/Learning_ripple_iti_vs_feedback/contrast_emmeans_ACC_by_phase_Learning_ripple_iti_vs_feedback.csv');
eemean  = sortrows(eemean ,'phase')
P  = nan(4);
P(1,2)  = contrast_01.p_value(1);  
P(3,4)  = contrast_01.p_value(2); 
P(1,3)  = contrast_fi.p_value(1); 
P(2,4)  = contrast_fi.p_value(2); 
P(P>0.05)=nan;
P = plt.f_symmetricMatrix(P);


axes(plt.H,"Position",[0.05 0.7 0.2 0.2])
plt  = myFigure(eemean.emmean',[],plt.H); 
plt.X = [1 2 4 5]; plt.Colors = [plt.CM.correct;plt.CM.white];
plt.addBar('E',eemean.SE','BarEdgeColor',plt.Colors*0.7,'P',P,'PLineOffset',0.01);
xticks([1.5 4.5])
plt.ylim([0.25 0.35])
xticklabels({'feedback' 'interval'})
ylabel('Ripple rate (events/sec)')
plt.setFont


% plot 2, ripple rate change with trials
interAct = readtable(fullfile(pwd,'Results/Learning_Ripple_Trial', ...
    'interaction_interval_ACCtrial_fit.csv'));
interAct  =sortrows(interAct,2);
interP = readtable(fullfile(pwd,'Results/Learning_Ripple_Trial', ...
    'interaction_interval_ACCtrial_slope.csv'));

d1=dtRippleM(dtRippleM.ACC==1,{'interval' 'trial' 'subject'});
d1.ripple = d1.interval;  d1.interval=[];
d1=d1(~isnan(d1.ripple),:);
d1m=groupsummary(d1,{'trial' 'subject'},'mean'); d1m.subject=[];
d1m=groupsummary(d1m,{'trial' },'mean');
 
ax = axes('Position',[0.4  0.7 0.2 0.2]);
hold on
plt.DATA = [d1m.trial,d1m.mean_mean_ripple];
plt.Colors = plt.CM.Flearning;
d =interAct(interAct.ACC==1,:);  
plt.addScatter([],[])
hs=plt.addShadeErrorBar(d.trial,d.fit,d.se);
plot(d.trial,d.fit,'color',plt.Colors,'LineWidth',2);
hs.patch.FaceColor = plt.Colors;
hs.patch.FaceAlpha=0.4;

hf=hs.patch; 
hs.mainLine.Color = plt.CM.Flearning; 

xlim([0,45]);ylim([0.1 0.5]);xticks([0,45]),yticks(0.1:0.1:0.5)
box off,xlabel('Trial'),ylabel('Ripple rate (Hz)'), title('Interval','FontWeight','normal')
plt.addPstar(['\beta = ' sprintf('%.4f',interP.Est_(2))],interP.p(2))
plt.setFont()


% plot 3, correlation between learn-R and feature inference performance
slope = learnR_results.Estimate(6); int = learnR_results.Estimate(5);
p = learnR_results.Pr___t__(end);
plt = myFigure([subjLearnR.interval, subjLearnR.eleInfer_perf],[],plt.H);
ax = axes('Position',[0.785 0.7 0.2 0.2]);
plt.Colors=plt.CM.FeleInfer;
hl=plt.addScatter([],[]);hold on
plt.ylim([0.6 1]);  
plt.addShadeErrorForGLM()
ylabel('Feature inference'); xlabel({'Learning-Increased Ripple'})
plt.addPstar(['\beta ' sprintf('= %.2f',slope)],p  )

title('Interval','FontWeight','normal')
plt.setFont
plt.savefig('Figures/Figure 1.pdf')

%% ###################################
%% visualization for supplementary figure 2A-J
H = figure(2); clf
 
plt = myFigure([],[],H);
load('../data/mat/raster_during_learning.mat')
ax = axes('Position',[0.05 0.7 0.4  0.2]);
plot([-3 5.5],repmat(raster.compare1vbsl.avgbsl,2),':', ...
    'Color',plt.CM.correct*0.5,'LineWidth',1 )
hold on
shadedErrorBar(raster.cond1.binscenters,raster.cond1.avg,raster.cond1.stdev,...
        {'color',plt.CM.white*0.6,'linewidth',1},1)

shadedErrorBar(raster.cond2.binscenters,raster.cond2.avg,raster.cond2.stdev,...
        {'color',plt.CM.correct,'linewidth',1},1)
xlim([-3 5.5]),xticks([-3 0 3 5])

box off
ylim([0 0.5])
plot([0 0],get(gca,'YLim'),':k','LineWidth',1)
plot([3 3],get(gca,'YLim'),':k','LineWidth',1)
plot([5 5],get(gca,'YLim'),'-k','LineWidth',0.2)

cms  = [plt.CM.white*0.6;plt.CM.correct];
for i = 1:2
    sig = raster.compare1vbsl.sig(i,:);
    min_cluster_pval = raster.compare1vbsl.min_cluster_pval(i,:);
    if min_cluster_pval>0.05,continue;end
    binscenters = raster.(sprintf('cond%d',i)).binscenters;
    idx=sprintf('%d',sig);
    onsets = binscenters(regexp(idx, '1{2,}', 'start'));
    offsets = binscenters(regexp(idx, '1{2,}', 'end'));
    tmp = get(gca,'Ylim');
    X=[onsets;offsets];

    Y=repmat(tmp(1)+range(tmp)*0.86 + range(tmp)*0.06*i,size(X));    
    line(X,Y,'LineWidth',4,'Color',cms(i,:))
    text(max(X)+0.01*range(get(gca,'xlim')),mean(Y), ...
        sprintf('p<%.4f',min_cluster_pval),'fontsize',4,'horizontalalignment','left'); 
end
ylabel('Ripple rate (Hz)')
xlabel('Time from response (s)')
%
ax = axes('Position',[0.55  0.7 0.35 0.2]);
xx = raster.xx; yy = raster.yy; nn = raster.nn;
binsizeX = xx(2)-xx(1);  binsizeY = yy(2)-yy(1);
CM = flipud(brewermap(32,'RdYlBu'));
tmph=imagesc(xx+binsizeX/2, yy+binsizeY, nn');
colormap(CM)
set(gca,'YDir','normal')
hold on
plot(get(gca,'XLim'),raster.split(2,:),'--k')
plot([0 0],get(gca,'YLim'),':k','LineWidth',1)
plot([3 3],get(gca,'YLim'),':k','LineWidth',1)
plot([5 5],get(gca,'YLim'),'-k','LineWidth',0.2)
set(gca,"CLim",[0 7])
xlim([-3 5.5]),xticks([-3 0 3 5])
box off
text(5.7,raster.split(2,1)/2,'error','Rotation',90,'horizontalalignment','center')
text(5.7,(max(get(gca,'ylim'))+raster.split(2,1))/2,'correct','Rotation',90,'horizontalalignment','center')
ytck = get(gca,"YTick"); ytcklbl = ytck./1000 + "K";
yticklabels(ytcklbl)
ylabel(sprintf('Single trials\n(pooled across elec.)'))
cb=colorbar(ax,'eastoutside','Position',[0.93 0.75 0.015 0.1]);
xlabel(cb,'Ripple count')
xlabel(ax,'Time from response (s)')
plt.setFont
% 
dt_chan = groupsummary(dtRipple,{'channel','phase'},@nanmean,'ripple');
dt_chan = unstack(dt_chan,"fun1_ripple","phase");

w=0.16;
Phases = {'stimulus','feedback'};
yloc   = [0.4 0.1];
for ii = 1:2
    phase = Phases{ii};
    % histogram
    ax = axes('Position',[0.05 yloc(ii) w w]);
    histogram(dt_chan.(phase) ,0.0:0.05:0.6,'FaceColor',plt.CM.gray),box off
    ylim([0 50])
    ylabel('Contact')
    xlabel('Ripple rate (Hz)')
    title([upper(phase(1)),phase(2:end)],FontWeight="normal")
    plt.setFont

    % bar plot
    t = readtable(sprintf('Results/Learning_ripple_ACC_phase/emmeans_ACC_%s.csv',phase));
    RRacc = t.emmean([2 1]);e = t.SE([2,1]);
    t = readtable(sprintf('Results/Learning_ripple_ACC_phase/contrast_ACC_%s.csv',phase));
    p = t.p_value; pp = nan(2); pp(eye(2)==0)=p;
    ax = axes('Position',[0.3 yloc(ii) w w]);
    plt.DATA=RRacc'; plt.E = e'; plt.X = [1 1.5]; plt.Colors=[plt.CM.correct;plt.CM.white];
    plt.X = [1 2 ];
    bh=plt.addBar('E',e,'P',pp,'PStarShowNS',true,'PStarOffset',0.02,'BarWidth',0.5, ...
        'BarEdgeColor',plt.Colors*0.7);
    plt.ylim([0.2 0.4])
    xticks(plt.X),xticklabels({'Correct','Wrong'})
    ylabel('Ripple rate (Hz)'),title([upper(phase(1)),phase(2:end)],FontWeight="normal")


    % scatter
    ax = axes('Position',[0.55 yloc(ii) w w]);
    interAct = readtable(fullfile(pwd,'Results/Learning_Ripple_Trial', ...
        sprintf('interaction_%s_ACCtrial_fit.csv',phase)));
    interAct  =sortrows(interAct,2);
    interP = readtable(fullfile(pwd,'Results/Learning_Ripple_Trial', ...
        sprintf('interaction_%s_ACCtrial_slope.csv',phase)));

    d1=dtRippleM(dtRippleM.ACC==1,{phase 'trial' 'subject'});
    d1.ripple = d1.(phase);  d1.(phase)=[];
    d1=d1(~isnan(d1.ripple),:);
    d1m=groupsummary(d1,{'trial' 'subject'},'mean'); d1m.subject=[];
    d1m=groupsummary(d1m,{'trial' },'mean');
    hold on
    plt.DATA = [d1m.trial,d1m.mean_mean_ripple];
    plt.Colors = plt.CM.Flearning;
    d =interAct(interAct.ACC==1,:); hold on
    plt.addScatter([],[])
    hs=plt.addShadeErrorBar(d.trial,d.fit,d.se);
    
    plot(d.trial,d.fit,'color',plt.Colors,'LineWidth',2);
    hs.patch.FaceColor = plt.Colors;
    hs.patch.FaceAlpha=0.4;

    hf=hs.patch;
 
    
    xlim([0,45]);ylim([0.1 0.45]);xticks([0,45]),yticks(0.15:0.1:0.45)
    box off,xlabel('Trial'),ylabel('Ripple rate (Hz)'),
    plt.addPstar(sprintf('\\beta = %.2e ',interP.Est_(2)),interP.p(2))
    title([upper(phase(1)),phase(2:end)],FontWeight="normal")


    % bar plot
    t = readtable('Results/LearnR_predict_infer/learnR_predict_feature_performance.csv');
    t = t(ismember(t.Var1,phase),:);
    ax = axes('Position',[0.8 yloc(ii) w w]);
    plt.X = [1 2 ]; plt.DATA = t.Estimate'; plt.E = t.Std_Error'; pp = t.Pr___t__';
    plt.Colors = [plt.CM.correct;plt.CM.correct];
    bh=plt.addBar('P',pp,'PStarShowNS',true,'PStarOffset',0.2,'BarWidth',0.5, ...
        'BarEdgeColor',plt.Colors*0.7);
    plt.ylim([-1 1])
    ylabel('Predicted performance (\beta)')
    xticklabels({'memory','feature inference'})
    title([upper(phase(1)),phase(2:end)],FontWeight="normal")

end
plt.setFont
plt.savefig('Figures/Figure S2.pdf')

%% ###################################

%% For Figure S2K-M

clear

load('../data/mat/HFB_activity_during_learning.mat')

% 
restoredefaultpath
addpath('../utils/')
addpath('../toolboxs/fieldtrip-20220714/')
addpath(genpath('../toolboxs/eeglab2022.1/functions/'))
ft_defaults

plt = myFigure;
AllStateMasks = struct;

% 
%%

clf
drawnow

rois = ["Y7_DorsalAttention"  "Y7_Frontoparietal"  "Hippocampus"];

isubplot = 1: numel(rois);
isubplot = isubplot-2;

for i = 1:3
    roi = rois(i);
    droi = threeROI(threeROI.roi==roi ,:) ;
    epoch_01 = groupsummary(droi,["eventACC" "roi" "subject"],'mean','data');
    epoch_01.GroupCount = [];
    epoch_01  = unstack(epoch_01, 'mean_data','eventACC');
    erps_cell = {epoch_01.correct';epoch_01.wrong'};


    for c= 1:2
        for ich = 1:size(erps_cell{c},2)
            erps_cell{c}(:,ich) = movmean(erps_cell{c}(:,ich),50,1);
        end
    end

    col = [    0.2627    0.5765    0.7647; 0.8392    0.3765    0.3020];
    ax = subplot(8,2,i);
    [hs, rg1 ]= plt.shadeErrorBar_matlab(Times,erps_cell{1}',[],col(1,:),ax);
    hs.mainLine.LineWidth  = 1; hs.shade.FaceAlpha=0.25;
    [hs, rg2 ] = plt.shadeErrorBar_matlab(Times,erps_cell{2}',[],col(2,:),ax);
    hs.mainLine.LineWidth  = 1; hs.shade.FaceAlpha=0.25;
    ylim([-0.04 0.1])

 
    if  isfield(AllStateMasks,roi)
        condmask = AllStateMasks.(roi) ;
    else
        [pcond, pgroup, pinter, statscond, statsgroup, statsinter] = std_stat(erps_cell, ...
            'condstats', 'on', 'paired',{'on'} , 'groupstats', 'off', ...
            'fieldtripmcorrect', 'cluster', ...
            'fieldtripmethod', 'montecarlo', ...
            'mode', 'fieldtrip' );
        
        condmask = pcond{1}<0.05;
        AllStateMasks.(roi) = condmask;
    end
    

    tnan = Times; tnan(~condmask)=nan;
    plot(ax,tnan,ones(size(tnan))*ax.YLim(2),'LineWidth',5,'Color',[0.5 0 0])
 

    title(ax,strrep(roi,'_',' '))
    xlim(ax,Times([1 end]))
    
    hold on
    plot(ax,[0 0],ax.YLim,'--k')
    plot(ax,[3 3],ax.YLim,'--k')
    plot(ax,[5 5],ax.YLim,'--k')
    plot(ax,ax.XLim,[0 0],'--k')
    xticks(ax,[0 3 5])
    box(ax,'off')
    plt.setFont
    
    ax.Position(3:4) = ax.Position(3:4).*[0.9  0.8];
    xlabel('Time from response onset (s)')
    ylabel('HFB power (dB)')
    drawnow
end


axs = findall(plt.H,'type','axes');
axs(1).Position(2) = axs(1).Position(2)-0.03;
%%

plt.savefig('Figures/Figure S2K-M.png')


 