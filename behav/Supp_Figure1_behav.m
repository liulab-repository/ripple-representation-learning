% For Figure S1

clear all
addpath('utils')
%%
sBehavTable = readtable('../data/table/Supp_behav_feature_dimension_performance.csv')
plt = myFigure;


%%
clf
% Figure S1A
ax=subplot(4,5,1);
histogram(sBehavTable.nLearningSession,15,'FaceColor',plt.CM.gray)
ylim([0    40])
xticks(1:5:11)
box off
ylabel('Number of participants')
xlabel('Number of learning sessions')


% Figure S1B
subplot(4,5,2)
plt.DATA = sBehavTable(:,1:3);
plt.Colors = repmat([plt.CM.white;],3,1);
plt.MarkerSize = 10;
hb=plt.addSpreadWithBar();
plt.addLineForPairedSpread('color',[plt.CM.gray 0.5])
ylim([0 1])
ylabel('Feature test accuracy')
xticklabels( regexprep(plt.VarName(1:end),{'feature_ACC_'  },{'' ''}) )
xlabel('Feature dimension')


% Figure S1C
plt.DATA = sBehavTable(:,4:6) ;
subplot(4,5, 3 )
plt.Colors = repmat([plt.CM.white;  ],3,1);
plt.MarkerSize = 10;
hb=plt.addSpreadWithBar();
plt.addLineForPairedSpread('-','color',plt.CM.gray)
ylim([0 1])
ylabel('2D inference accuracy')
xticklabels(strrep(plt.VarName(1: end),{'irr_'},{''}))
xlabel('Irrelevant dimension')

%%  Figure S1D

formula = 'ACC~trialType*session+ (1|subject)' ;
d1 = readtable('../data/table/learning_behav_without_boundary.csv');
[LME2] = fitLME_byR(d1,formula)

ord = [4 3 2 1];
ax = subplot(4,5,4); cla; w = ax.Position(3);
assert(isequal(LME2.interaction.eemeans.trialType(ord(1:2)),{'memory';'inference'}))
plt.DATA =  LME2.interaction.eemeans.emmean(ord )';
plt.X = [1 2 4 5];
plt.E =  LME2.interaction.eemeans.SE(ord )';
plt.P =  LME2.interaction.pMat(ord,ord);
plt.P(plt.P>0.05)=nan;
plt.Colors = repmat([plt.CM.FeleMemory;plt.CM.FeleInfer],2,1);
hb=plt.addBar('pLineOffset',0.05);
plt.ylim([0.5 1 ])
hl = legend(hb(1:2),{'memory' 'feature inference'},'location','eastoutside','box','off');
hl.ItemTokenSize = [9 18]';
ax.Position(3) = w;
xticks([1.5 4.5])
xticklabels({'before' 'after'})
ylabel('Feature test performance')
xlabel('Session')
ylim([0.5 1])
plt.setFont
plt.subplot_compress_margins(0.1,0.1)

drawnow
%%

clearvars -except plt
 
%% For Figure S1E-F
% First, run Ripple_offline_learning to create ../ripple/Results/dt_trial_level_Ripple.csv.

dtLearning1 = readtable('../ripple/Results/dt_trial_level_Ripple.csv')
plt = myFigure;

LME_real = fitLME_byR(dtLearning1,'interval~trial*ACC+(1|subject/channel)','trial',"real");

real_P = LME_real.interaction.eemeans.p(3:4)'
real_T = LME_real.interaction.eemeans.tVal_(3:4)'
real_beta = LME_real.interaction.eemeans.Est_(3:4)'
%%

if exist(fullfile('Results','Outlier_ALLLME.mat'),'file')
    load(fullfile('Results','Outlier_ALLLME.mat'))
else
    dtL1_Error = dtLearning1(dtLearning1.ACC==0,:) ;
    dtL1_Correct = dtLearning1(dtLearning1.ACC==1,:) ;

    N_Correct = height(dtL1_Correct)
    N_Error   = height(dtL1_Error)
    N_Test = N_Error

    idx_Correct = find(dtLearning1.ACC==1);
    ALLLME = {}; Tvalues = [];
    parfor i = 1:10
        idx_rng =  sort(randperm(N_Correct,N_Test)');
        d = cat(1,dtL1_Error, dtL1_Correct(idx_rng,:));
        tic
        LME = fitLME_byR(d,'interval~trial*ACC+(1|subject/channel)','trial',"tmp"+i,"tmp/tmp"+i);
        toc

        ALLLME{i} = LME;
    end
    mkdir("Results")
    save(fullfile('Results','Outlier_ALLLME.mat'),'ALLLME')

end


Tvalues = arrayfun(@(x) ALLLME{x}.interaction.eemeans.tVal_(3:4)',1:numel(ALLLME),'un',0);
Tvalues  = cat(1,Tvalues{:});

beta = arrayfun(@(x) ALLLME{x}.interaction.eemeans.Est_(3:4)',1:numel(ALLLME),'un',0);
beta = cat(1,beta{:});

pvals= arrayfun(@(x) ALLLME{x}.interaction.eemeans.p(3:4)',1:numel(ALLLME),'un',0);
pvals = cat(1,pvals{:});

%%

[~,minidx]  = min(abs(pvals(:,2)-0.05))
[~,minidx0] = min(abs(pvals(:,1)-0.05))
 
ax = subplot(4,2,3);
ax.Position(2) = ax.Position(2)-0.05;
h1=histogram(Tvalues(:,1),10,'FaceColor',plt.CM.wrong,'EdgeAlpha',0.001)
hold on
h2=histogram(Tvalues(:,2),10,'FaceColor',plt.CM.correct,'EdgeAlpha',0.001)
xlim([-5 5])
box off
xlabel('T value'), ylabel('Number of permutations')
mx = min(Tvalues(minidx,2));
plot(mx*[1 1],ax.YLim,'--k')
text(mx,ax.YLim(end),0,sprintf('p = %.2f',max(pvals(minidx,2))),'HorizontalAlignment','center','VerticalAlignment','bottom')
plot(repmat(real_T,2,1),ax.YLim,'-k')
plt.setFont
hl = legend([h1,h2],{'Error' 'Correct'},'Location','northwest','Box','off');
hl.ItemTokenSize = [9 18]';

ax = subplot(4,2,4);
ax.Position(2) = ax.Position(2)-0.05;
h1=histogram(beta(:,1),10,'FaceColor',plt.CM.wrong,'EdgeAlpha',0.001)
hold on
h2=histogram(beta(:,2),10,'FaceColor',plt.CM.correct,'EdgeAlpha',0.001)
xlim([-3 3]*0.001)
box off
xlabel('\beta value'), ylabel('Number of permutations')
mx = min(beta(minidx,2));
plot(mx*[1 1],ax.YLim,'--k')
text(mx,ax.YLim(end),0,sprintf('p = %.2f',max(pvals(minidx,2))),'HorizontalAlignment','center','VerticalAlignment','bottom')
plot(repmat(real_beta,2,1),ax.YLim,'-k')
plt.setFont
hl = legend([h1,h2],{'Error' 'Correct'},'Location','northwest','Box','off');
hl.ItemTokenSize = [9 18]';
 
plt.setFont
plt.savefig('Results/Figure S1')


mean(pvals<0.05 & beta>0)
mean(Tvalues)
 

 
