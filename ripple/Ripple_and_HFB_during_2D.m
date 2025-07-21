

clear,clc
addpath('../utils/')

dtripple = readtable('../data/dd_subject_rest_and_overall_ripple_rate.csv');
dtacc    = readtable('../data/dd_subject_performance.csv');
 
%% compare performance of adjacent inference and non-adjacent inference
plt = myFigure; H = plt.H;
%
clf
dd = dtacc{:,{'compAdjacent','compNonAdjacent'}};

ax = axes('Position',[0.1 0.7 0.2 0.2]);
plt.DATA = dd; plt.X = [1 3];
plt.Colors=[plt.CM.FadjacentInfer;plt.CM.FnonadjacentInfer];
plt.addBar('P',[0.5 0.5; 0.5 0.5],'PStarOffset',0.2,'PLineOffset',0.05, ...
    'Barwidth',0.7,'BarEdgeColor',[.5 .5 .5],'PstarshowNS',true),hold on
plot(plt.X,plt.DATA,'Color',[1 1 1]*0.7)
scatter(plt.X,plt.DATA,30,[1 1 1]*0.6,'filled','MarkerEdgeColor',[1 1 1])
ylim([0 1]),ylabel('Accuracy')
xticklabels({'adjacent','non-adjacent' })

%% correlation between ripple rate during rest and 2d behavior

dd = [dtripple.RR_Rest,dtacc.compAdjacent];
[r,p]=partialcorr(dd(:,1),dd(:,2),dtripple.RR_Overall);
ax = axes('Position',[0.45 0.7 0.2 0.2]);
plt.DATA = dd;
plt.Colors = plt.CM.FadjacentInfer;
hold on
plt.addPartialRegressionScatter()

xlabel({'Ripple rate during rest (Hz)','(adjusted)'}),ylabel({'Adjacent inference' '(adjusted)'})
ylim([-0.3 0.3]);yticks(-0.3:0.1:0.3)
text(0.1,-0.2,0,sprintf('r=%.2f n.s.',r))


%
dd = [dtripple.RR_Rest,dtacc.compNonAdjacent];
[r,p]=partialcorr(dd(:,1),dd(:,2),dtripple.RR_Overall);

ax = axes('Position',[0.78 0.7 0.2 0.2]);
plt.DATA = dd;
plt.Colors = plt.CM.FnonadjacentInfer;
hold on
plt.addPartialRegressionScatter()

xlabel({'Ripple rate during rest (Hz)','(adjusted)'}),ylabel({'Non-adjacent inference' '(adjusted)'})
ylim([-0.3 0.3]);yticks(-0.3:0.1:0.3)
text(0.1,-0.2,0,sprintf('r=%.2f n.s.',r))
 
%% raster
 
figure(plt.H)
load('../data/OnTask2D_Stim2Resp_raster.mat')
 
%
ax = axes('Position',[0.07  0.35 0.22 0.2]);
xx = raster.xx; yy = raster.yy; nn = raster.nn;
binsizeX = xx(2)-xx(1);  binsizeY = yy(2)-yy(1);
CM = flipud(brewermap(32,'RdYlBu'));
tmph=imagesc(xx+binsizeX/2, yy+binsizeY, nn');
colormap(CM)
set(gca,'YDir','normal')
hold on
for i = 2:4
    plot(get(gca,'XLim'),raster.split(i,:),'--k','LineWidth',1)
end
plot([0 0],get(gca,'YLim'),'-k','LineWidth',1)
for i = 1:size(raster.marker,1)
    scatter(-raster.marker{i,1},raster.marker{i,2},1,'o', ...
    'markerFaceColor',[0 0 0],'markerEdgeColor','none'); hold on
end
set(gca,"CLim",[0 8])
xlim([-8 0.5]),xticks([-8 0 0.5])
box off

ytck = get(gca,"YTick"); ytcklbl = ytck./1000 + "K";
yticklabels(ytcklbl)
ylabel(sprintf('Single trials\n(pooled across elec.)'))
cb=colorbar(ax,'northoutside','Position',[0.1 0.56 0.06 0.01 ]);
cb.set("Ticks",[0 8])
cb.set("Ticks")
xlabel(cb,'Ripple count')
xlabel(ax,'Time from response (s)')
plt.setFont
%% Difference in ripple rate during the 2D task. corresponding to raster plot
cbt = readtable('../data/dd_RR_during_2D_withinResp.csv')
 
plt.DATA = cbt{:,2:end}; 
ax = axes('Position',[0.45 0.35 0.2 0.2]); 
plt.Colors = [plt.CM.FadjacentInfer;plt.CM.white; 
              plt.CM.FnonadjacentInfer;plt.CM.white; ];

plt.X = [1 2 4 5];
cla
plt.addBar('barwidth',0.8,'BarEdgeColor',[.5 .5 .5])
% plt.addSpreadWithBar('barwidth',0.8,'BarEdgeColor',[.5 .5 .5])
xx = repmat(plt.X,size(plt.DATA,1),1);
scatter(xx(:),plt.DATA(:),30,[1 1 1]*0.6,'filled','MarkerEdgeColor',[1 1 1])
hold on
hl=plot([1 2],plt.DATA(:,1:2),'-','Color',[0.6 0.6 0.6 0.3]);
plot([4 5],plt.DATA(:,3:4),'-','Color',[0.6 0.6 0.6 0.3])
plt.ylim([0 1.2])
ylabel('Ripple rate (Hz)'),xticks([1.5,4.5]),xticklabels({'memory' 'inference'})
xlabel('')
title('Within response','FontWeight','normal')
%% relationship between hfb and performance during 2d inference task.
b = readtable('../data/peri-rippleHFB_during_compound_inference/tvalues_peri-rippleHFB_during_compound_inference.csv');
ax = axes('Position',[0.78 0.35 0.2 0.2]); 
plt.DATA = b.Estimate(2);
plt.Colors = plt.CM.FnonadjacentInfer;
plt.addBar('E',b.Std_Error(2),'P',b.Pr___t__(2),'barwidth',0.4, ...
    'BarEdgeColor',[.5 .5 .5],'PstarShowNS',true,'PstarOffset',0.001)
ylim([0 0.01]),yticks([0 0.01])
ylabel('Predicted non-adjacent inference (\beta)'),
title('Within response','FontWeight','normal')
xticklabels('Peri-ripple HFB within DMN')
% text(0.1,0.0095,'DMN','FontWeight','normal')
%%
plt.setFont
sfilename = fullfile('Figures','supp_fig_ontask_2D.pdf');
plt.savefig(sfilename)