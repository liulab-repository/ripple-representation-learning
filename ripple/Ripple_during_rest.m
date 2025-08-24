% For Figure 3

clear,clc
restoredefaultpath                
addpath('../utils/')

% read data
dtRipple = readtable('../data/table/dd_subject_rest_and_overall_ripple_rate.csv');
dtACC    = readtable('../data/table/dd_subject_performance.csv');
%% Enhanced ripple rate during rest predict enhancement of feature inference
% prepare data for R
dt = dtACC(:,{'testAfter','testAfterInfer','testAfterMemory', ...
             'testBefore','testBeforeInfer','testBeforeMemory', ...
             'subject'});

dt = join(dt,dtRipple,"Keys",'subject');
writetable(dt,'Results/dt_rest_ripple_feature_behave.csv')


% test whether performance changed after rest
% (This effect is more pronounced at the trial level (LME), and here we
% statistically analyze at the subject level to align with subject-level correlation analysis.)
[p0,~,stat0]=signrank(dt.testAfterMemory, dt.testBeforeMemory ,'Tail','right') ;
[p1,~,stat1]=signrank(dt.testAfterInfer, dt.testBeforeInfer ,'Tail','right')   ;

% test whether ripple rate enhanced during rest
[h,p,ci,stat] = ttest(dt.RR_Rest,dt.RR_Overall);


% run robust fit in R
Rfile = fullfile(pwd,"R_scripts/Increased_ripple_during_rest_facilitate_feature_inference.R");
RUN_R_Script('',Rfile)

% read results from R
report1 = readtable('Results/Rest_Ripple_feature/lmrob_testMemory.csv');
report2 = readtable('Results/Rest_Ripple_feature/lmrob_testInfer.csv') ;

%% visualization, main figure 3
dtchan = readtable('../data/table/dd_channel_rest_and_overall_RR.csv');
H = figure(3);clf
plt = myFigure([],[],H);  w= 0.22;

% plot 1, ripple rate distribution during rest and the entire task duration
ax = axes('Position',[0.05 0.7 w w]);
x = sort(dtchan.Overall,'ascend');xi = linspace(min(x),max(x),100);
hh1=histogram(x,0:0.03:0.6,'FaceColor',plt.CM.gray,'EdgeAlpha',0, ...
    'FaceAlpha',0.9,'Normalization','pdf');box off, hold on
[f,xi] = ksdensity(x,xi,'Bandwidth',0.03);
plot(xi,f,"Color",plt.CM.gray*0.9,'LineWidth',1)

x = sort(dtchan.Rest,'ascend'); xi = linspace(min(x),max(x),100);
hh2=histogram(x,0:0.03:0.6,'FaceColor',plt.CM.rest,'EdgeAlpha',0, ...
    'FaceAlpha',0.9,'Normalization','pdf');box off, hold on
[f,xi] = ksdensity(x,xi,'Bandwidth',0.03);
plot(xi,f,"Color",plt.CM.gray*0.6,'LineWidth',1)

xlim([0 0.7])
ylabel('Contacts (probability density)')
xlabel('Ripple rate (Hz)')
hl=legend([hh1,hh2],{'overall','rest'},'location',[0.2 0.65+w 0.1 0.05],'Box','off');
hl.ItemTokenSize = [9 18]'; 
 
% plot 2, bar plot for ripple rate during rest and entire task duration
ax = axes('Position',[0.4 0.7 w w]);
plt.DATA = [dt.RR_Overall,dt.RR_Rest];
ind = abs(zscore(plt.DATA,1,1))>3;
plt.DATA(sum(ind,2)>0,:)=[];
[h,p,ci,stat]=ttest(plt.DATA(:,1),plt.DATA(:,2)) ;
pp = nan(2); pp(1,2)=p; pp(2,1)=p;
plt.P = pp;
plt.X  = [1 2];
plt.Colors = [plt.CM.gray;plt.CM.rest;];
hb=plt.addBar('barwidth',0.35,'Plineoffset',0.17,'PstarOffset',0.03);
hold on
plot(plt.X,plt.DATA,'-','Color',[1 1 1]*0.75)

s1=scatter(repmat(plt.X(1),size(plt.DATA,1)),plt.DATA(:,1),25,plt.CM.gray, ...
    'filled','CData',plt.CM.gray_light,'MarkerFaceAlpha',1,'MarkerEdgeColor',plt.CM.gray*0.9);
s1=scatter(repmat(plt.X(2),size(plt.DATA,1)),plt.DATA(:,2),25,plt.CM.gray, ...
    'filled','CData',plt.CM.gray_light,'MarkerFaceAlpha',1,'MarkerEdgeColor',plt.CM.gray*0.9);
ylim([0.1 0.6])
xticklabels({'overall','rest'})


% plot 3, bar plot for behavioural performance before and after rest
ax = axes('Position',[0.05 0.35 w w]);
pp = nan(4); pp(1,3)=p0;pp(2,4)= p1; plt.P = pp;
plt.DATA = dt(:,{'testBeforeMemory','testBeforeInfer','testAfterMemory','testAfterInfer'}).Variables;
plt.X = [1 2 4 5]; plt.E = plt.f_sem(); plt.P = pp;
plt.Colors = repmat([plt.CM.FeleMemory;plt.CM.FeleInfer],2,1);
hb=plt.addBar('barWidth',0.7,'PStarShowNS',1,'Plineoffset',0.02,'Pstaroffset',0.01);
plt.legend(hb(1:2),{' memory',' feature inference'},'location',[0.1 0.33+w 0.15 0.05])
ylim([0.8 1])
yticks(0.8:0.1: 1)
ylabel('Accuracy')
xticks([1.5 4.5]),xticklabels({'before','after'})

% plot 4, scatter plot for the correlation between ripple enhancement and
% behavioural enhancement
ax = axes('Position',[0.4 0.35 w w]);
plt.DATA = [dt.RR_Rest_enhance,dt.testAfterInfer-dt.testBeforeInfer];
plt.Colors = plt.CM.FeleInfer;
plt.addScatter(report2.Estimate(2),report2.Estimate(1))
plt.addShadeErrorForGLM
xlabel('Increased ripple rate during rest')
ylabel('Enhanced feature inference')
plt.ylim([-0.4 0.4]),yticks(-0.4:0.2:0.4)
plt.addPstar(sprintf('\\beta=%.2f',report2.Estimate(2)),report2.Pr___t__(2))

plt.setFont
plt.savefig('Figures/Figure 3.pdf')
 
%% visualization for supplimentary figure 2

phases = {'Rest','Overall'};
histlim = [30 50];
for iph = 1:numel(phases)
    phase = phases{iph}; 
    rr = dtchan.(phase); 
    hh = load('../data/mat/hippocampal_channels_view.mat');
    climv=[0.0 0.6];
    CM = flipud(brewermap(64,'Reds'));
    CM = flip(CM);
    DD = rr;
    DD(DD>climv(2))=climv(2); DD(DD<climv(1))=climv(1);

    ncolors = size(CM,1);
    for ii = 1:size(DD,1)
        colind = floor(interp1(climv,[1 ncolors],DD(ii),'linear',ncolors));
        set(hh.Helec(ii),'facecolor',CM(colind,:));
    end
    hc = colorbar(gca,'Position',[0.85 0.2 0.025 0.5 ],'Ticks',climv);
    clim(climv)
    colormap(hc,CM)
    ylabel(hc,'Ripple rate (Hz)')
    fontname(hc,'Arial'),fontsize(hc,8,'points')

    ax = axes('Position',[0.25 0.55 0.35 0.3]);
    histogram(rr,0:0.03:0.6,'FaceColor',plt.CM.gray*0.1)
    box off,ylim([0 histlim(iph)]),yticks([0 histlim(iph)])
    ylabel('Contact'),xlabel('Ripple rate')
    title(sprintf('n = %d contacts (%d subjects)',numel(rr),size(dt,1)),'FontWeight','normal')
    plt.setFont(findobj(hh.Hes,"type",'axes'))
    plt.savefig(sprintf('Figures/hipp_ripple_view_%s.pdf',phase),hh.Hes)
end
%%


