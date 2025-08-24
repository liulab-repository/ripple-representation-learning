clear,clc,close all
restoredefaultpath
addpath('../utils/')

%% load data  

HFB = load('../data/mat/peri_ripple_HFB_dB.mat');
channel = readtable('../data/table/Channel_ROI.csv');
dtACC   = readtable('../data/table/dd_subject_performance.csv');
Ripple  = readtable('../data/table/dd_channel_rest_and_overall_RR.csv');
Ripple  = groupsummary(Ripple,'subject','mean');Ripple.GroupCount=[];

%% Select all channels that have a hippocampal channel in the same subject.

[~,~,periRRchan] = intersect(HFB.EEG_ALLchan.channel,channel.channel,'stable');
EEG_channel  = channel(periRRchan,:);
ROI_y7 = {'Y7_Default'	'Y7_DorsalAttention'	'Y7_Frontoparietal'	'Y7_Limbic'	'Y7_Somatomotor'	'Y7_VentralAttention'	'Y7_Visual'};
ROI_y7_dmn_sub = {'DMNTL'	'DMNTPJ'	'DMNDLPFC'	'DMNLOFC'	'DMNmPFC'	'DMNPCC'};
%% Set the ROI order based on the averaged peri-ripple HFB for visualization.

% Compute the average peri-ripple HFB within each ROI.
roiMeanHFB = table;
ChannelLevel =[];
for roi = [ROI_y7,ROI_y7_dmn_sub]  
    roi_ind = EEG_channel.(roi{1})>0;
    hfb  = mean(HFB.EEG_ALLEpoch(roi_ind,abs(HFB.times)<250),2);
    hfbg = HFB.EEG_ALL_generalAvg(roi_ind) - HFB.EEG_ALL_wholeAvg(roi_ind);
    roiname = repmat(roi,sum(roi_ind),1);
    subject = EEG_channel.subject(roi_ind);
    ChannelLevel = cat(1,ChannelLevel,table(hfb,hfbg,roiname,subject));
    roiMeanHFB.(roi{1}) = mean(hfb,'all') ;    
end
roiMeanHFB{2,:} = [ones(1,numel(ROI_y7)),ones(1,numel(ROI_y7_dmn_sub))*2];
[~,order_real]  = sortrows(roiMeanHFB{:,:}',[2 1],{'ascend' 'descend' });
roiMeanHFB_real = roiMeanHFB(:,order_real);

% highligh Default model network
roiMeanHFB.Y7_Default(1)=999;

% sort ROI order for visualization
[~,order]=sortrows(roiMeanHFB{:,:}',[2 1],{'ascend' 'descend' });

roiMeanHFB = roiMeanHFB(:,order); 
%% visualiztion, main figure 5,  time course of peri-ripple HFB
H = figure(5);clf
plt =myFigure([],[],H);

% attributes of subplots 
locs=[[0.05 0.3 0.4 0.5 0.3 0.4 0.5;
     [0.7 0.79 0.79 0.79 0.7 0.7 0.7]+0.1;],...
    [[0.3 0.4 0.5 0.3 0.4 0.5]+0.4;
     [0.79 0.79 0.79 0.7 0.7 0.7]+0.1]];
locs = array2table(locs,'VariableNames',roiMeanHFB.Properties.VariableNames);
wh = [0.08 0.06];
ylms = [repmat([-1 2],13,1) ] ;
subcolors = [nan(7,3);flip( brewermap(numel(ROI_y7_dmn_sub),'OrRd') )];
 
 
SubjectLevel =[];
ROILevel = table;
for n = 1:numel(roiMeanHFB.Properties.VariableNames)
    roi = roiMeanHFB.Properties.VariableNames(n);
    roi_ind = EEG_channel.(roi{1})>0;
    roi_chanTbl = EEG_channel(roi_ind,:);
    
    hfb = HFB.EEG_ALLEpoch(roi_ind,:);
    hfbg = HFB.EEG_ALL_generalAvg(roi_ind) - HFB.EEG_ALL_wholeAvg(roi_ind);
    [hfb,~ ] = groupsummary(hfb,roi_chanTbl.subject,@temporalRobmean);
    [hfbg,subject ] = groupsummary(hfbg,roi_chanTbl.subject,'mean');
    roiname = repmat(roi,numel(subject),1);
    SubjectLevel = cat(1,SubjectLevel,table(hfb,hfbg,roiname,subject));
    
    dd = zscore(hfb,1,2);  
    m = mean(dd,1);
    sem = std(dd)/sqrt(size(dd,1));
 
    if contains(roi,'Y7_')       
        col = plt.CM.(roi{1});
    else
        col = subcolors(n,:);
    end
    pos = [locs.(roi{1})',wh];
    if isequal(roi,{'Y7_Default'}),pos= [locs.(roi{1})',0.15 0.15];end
    axes('Position',pos,'Tag','timecourse'); 
    plt.addShadeErrorBar(HFB.times/1000,m,sem,{'color',col})
    xticks([]),ylim(ylms(n,:)),yticks([])
    hold on
    plot(HFB.times/1000,0*HFB.times,'--k')
    plot([0,0],get(gca,'ylim'),'--k')    
    roiname = plt.f_DMNsubNames(roi{1},true,true);
    title(sprintf('%s (n=%d)',roiname{1},sum(roi_ind)),'FontWeight','normal')    
    if isequal(roi,{'Y7_Default'})
        ylabel('HFB power (z-score)');     
        xticks([-1 0 1]),
        yticks([-1 0 1 2])
        title(sprintf('%s (n=%d, N=%d)',roiname{1},sum(roi_ind),numel(subject)),'FontWeight','normal')    
    end
    
    box off
    drawnow
end

findall(plt.H,'Tag','timecourse')
tcaxs = findall(plt.H,'type','axes');
fontsize(tcaxs,7,"points")
 
%% Prepare subject-level data 
% for correlation analysis between HFB activity and behavioral performance.

% Average activation within 250 ms of the ripple peak.
SubjectLevel.hfb_500 = mean(SubjectLevel.hfb(:,HFB.times>-250&HFB.times<250),2);
Sub_roi_hfb = unstack(SubjectLevel,'hfb_500','roiname','GroupingVariables','subject');

% Genearal HFB activation during rest
Sub_roi_hfb_g = unstack(SubjectLevel,'hfbg','roiname','GroupingVariables','subject');
Sub_roi_hfb_g.Properties.VariableNames(2:end) = strcat('g_',Sub_roi_hfb_g.Properties.VariableNames(2:end));


subjHfbBehav = innerjoin(Sub_roi_hfb,dtACC,"Keys","subject");
subjHfbBehav = innerjoin(subjHfbBehav,Sub_roi_hfb_g,"Keys",'subject');
subjHfbBehav = innerjoin(subjHfbBehav,Ripple,"Keys","subject"); 
 
mkdir('Results/Peri_Ripple_HFB_behav')
writetable(subjHfbBehav,'Results/Peri_Ripple_HFB_behav/dt_periRippleHFB_behave.csv')
 
 
%% Correlation analysis between HFB and behavior.
 
sum(channel(:,3:end).Variables)
sum(EEG_channel(:,3:end).Variables)
 
% use robust fit to calculate partial correlation
Rfile = fullfile(pwd,"R_scripts/robust_Corr_rest_periRippleHFB_behav.R");
RUN_R_Script('',Rfile)
 
% read results
resultsMain = readtable('Results/Peri_Ripple_HFB_behav/R_robCorr_periHFB_CompInfer.csv');
resultsMain.t_value = [];
resultsMain.z_value = norminv(1-resultsMain.p_value/2);
resultsMain.z_value(resultsMain.cor<0) = -resultsMain.z_value(resultsMain.cor<0);
resultsMain.SE = 1./sqrt(resultsMain.df+2-3);
% multiple comparison correction
ind_y7 = contains(resultsMain.x_var,'Y7');
ind_sig = resultsMain.p_value<0.05; 
resultsMain.p_value(ind_y7&ind_sig)  = resultsMain.p_value(ind_y7&ind_sig)  * 7;
resultsMain.p_value(~ind_y7&ind_sig) = resultsMain.p_value(~ind_y7&ind_sig) * 6;
resultsMain.p_flag = categorical(resultsMain.p_value<0.05,[1 0],["*" "-"])


resultsAdditional = readtable('Results/Peri_Ripple_HFB_behav/R_robCorr_periHFB_CompInfer_additionalControll.csv');
resultsAdditional.t_value=[];
resultsAdditional.z_value = norminv(1-resultsAdditional.p_value/2);
resultsAdditional.z_value(resultsAdditional.cor<0) = -resultsAdditional.z_value(resultsAdditional.cor<0);
resultsAdditional.SE = 1./sqrt(resultsAdditional.df+2-3);
Corr_xy_2_plot = readtable('Results/Peri_Ripple_HFB_behav/R_robCorr_periHFB_CompInfer_xy.csv');

%%
delete(setdiff(findall(plt.H,'type','axes'),tcaxs))
rd= resultsMain([order;order+13],:);
Y7names = rd.x_var(1:7);
Y7Cols = [];
for i = 1:7
    Y7Cols(i,:) =  plt.CM.(Y7names{i});
end
w = 0.2;
figure(plt.H)
ax = axes('Position',[0.05 0.5 w w]);
id =(1:7) ;
p = rd.p_value(id)';
plt.DATA = rd.z_value(id)';
plt.E = rd.SE(id)';
plt.Colors = Y7Cols;

bh=plt.addBar('P',p,'Pstaroffset',0.05);
plt.ylim([-3 4]),
xticklabels ''
plt.legend(bh,strrep(Y7names,'Y7_',''));ax.Position(3) = w;
title('Non-adjacent inference','FontWeight','normal')
ylabel('Correlation (z-score)')

ax = axes('Position',[0.45  0.5 w w]);
id =(1:7)+13;
p = rd.p_value(id)';
plt.DATA = rd.z_value(id)';
plt.E = rd.SE(id)';
plt.Colors = Y7Cols;
plt.addBar('P',p,'Pstaroffset',0.05)
plt.ylim([-3 4]),
xticklabels ''
title('Adjacent inference','FontWeight','normal')
ylabel('Correlation (z-score)')

ax = axes('Position',[0.78 0.5 w w]);
ind = ismember(Corr_xy_2_plot.x_var,'Y7_Default') & ismember(Corr_xy_2_plot.y_var,'compNonAdjacent');
plt.DATA = Corr_xy_2_plot{ind,{'res_x','res_y'}};
plt.Colors = plt.CM.Y7_Default;
plt.addScatter([],[])
plt.addShadeErrorForGLM()
xlabel({'Peri-ripple HFB power (dB) (adjusted)'})
ylabel({'Non-adjacent inference (adjusted)'})
title('DMN','FontWeight','normal')
plt.addPstar(sprintf('r = %.2f',rd.cor(1)),rd.p_value(1))


ax = axes('Position',[0.05 0.15 w w]);
id =(1:6)+7 ;
subnames = rd.x_var(id);
plt.DATA = rd.z_value(id)';
plt.E = rd.SE(id)';
plt.Colors = subcolors(end-5:end,:)    ;
p = rd.p_value(id)';

bh=plt.addBar('P',p,'Pstaroffset',0.05);
plt.ylim([-1 4]),
xticklabels ''
plt.legend(bh,plt.f_DMNsubNames(subnames));ax.Position(3) = w;
title('Non-adjacent inference','FontWeight','normal')
ylabel('Correlation (z-score)')

ax = axes('Position',[0.45  0.15 w w]);
id =(1:6)+7+13;
subnames = rd.x_var(id);
plt.DATA = rd.z_value(id)';
plt.E = rd.SE(id)';
p = rd.p_value(id)';
plt.addBar('P',p,'Pstaroffset',0.05)
plt.ylim([-1 4]),
xticklabels ''
title('Adjacent inference','FontWeight','normal')
ylabel('Correlation (z-score)')

ax = axes('Position',[0.78 0.15 w w]);
ind = ismember(Corr_xy_2_plot.x_var,'DMNmPFC') & ismember(Corr_xy_2_plot.y_var,'compNonAdjacent');
plt.DATA = Corr_xy_2_plot{ind,{'res_x','res_y'}};
plt.Colors = subcolors(end-2,:);
plt.addScatter([],[])
plt.addShadeErrorForGLM()
xlabel({'Peri-ripple HFB power (dB) (adjusted)'})
ylabel({'Non-adjacent inference (adjusted)'})
title('mPFC','FontWeight','normal')
iv = find(contains(rd.x_var,'DMNmPFC'),1);
plt.addPstar(sprintf('r = %.2f',rd.cor(iv)),rd.p_value(iv))
plt.FontSize = 7;
plt.setFont
plt.savefig('Figures/Fig5_peri_ripple_HFB.pdf')

%% visualization, peri-ripple HFB on surface

CH=load('../data/mat/Contact_distribution_DMN.mat')

EEG_ALLEpoch_500= mean(HFB.EEG_ALLEpoch(:,HFB.times>-250&HFB.times<250),2);
channelid = arrayfun(@(x) num2str(x),1:size(channel,1),'un',0)';
channelid = channelid(periRRchan);
[~,ia,ib] = intersect(channelid,CH.elecID,'stable');
bd = EEG_ALLEpoch_500(ia);
bH = CH.elecH(ib);
bids = CH.elecID(ib);

mycmap =  brewermap(32,'RdYlBu');
mycmap = flip(mycmap);

lm =[-0.1 0.1];
bd(bd<-0.1)=-0.1;
bd(bd>0.1)=0.1;
cols=f_interp_value2color(bd,mycmap(:,:) );

for ic = 1:numel(bH)
    id = bids{ic};
    ind = ismember(CH.elecID,id);
    set(CH.elecH(ind),'MarkerFaceColor',cols(ic,:) , ...
        'MarkerEdgeColor',cols(ic,:)*1,"MarkerSize",2)
    
end
disp done!
 
hb = colorbar("east",'Position',[0.93 0.3 0.02 0.4]);
hb.Ticks=lm;hb.Label.String = 'Peri-ripple HFB power (dB)';
colormap(mycmap)
clim(lm)
 
exportgraphics(CH.H,'Figures/Peri-ripple-HFB-Y7.png','Resolution',300)


%% visualization, supplementary figure 4, add additional covar
Corr_xy_2_plot_Ad = readtable('Results/Peri_Ripple_HFB_behav/R_robCorr_periHFB_CompInfer_additionalControll_xy.csv');

w=0.18;
H=figure(16);clf
plt = myFigure([],[],H);

hfbbin = unstack(ChannelLevel,"hfb",'roiname');
hfbbin = hfbbin(:,roiMeanHFB_real.Properties.VariableNames);

ax = axes('Position',[0.08 0.7 w w]);
roinames = hfbbin(:,contains(hfbbin.Properties.VariableNames,'Y7_')).Properties.VariableNames;
plt.DATA = hfbbin{:,roinames};
plt.Colors = reshape(plt.CM{:,roinames},3,[])';
[h,p,ci,stat] = ttest(plt.DATA);
plt.addBar('P',p,'ErrorbarStyle' ,'I' )
xticklabels(plt.f_DMNsubNames(roinames,false))
ylabel('Peri-ripple HFB (dB)')

ax = axes('Position',[0.3  0.7 w w]);
roinames = hfbbin(:,contains(hfbbin.Properties.VariableNames,'DMN')).Properties.VariableNames;
plt.DATA = hfbbin{:,roinames};
plt.Colors = subcolors(8:end,:) ;
[h,p,ci,stat] = ttest(plt.DATA);
plt.addBar('P',p,'ErrorbarStyle' ,'I' )
xticklabels(plt.f_DMNsubNames(roinames,false))
ylabel('Peri-ripple HFB (dB)')


ax = axes('Position',[0.55 0.7 w w]);
ind = ismember(Corr_xy_2_plot_Ad.x_var,'Y7_Default') & ismember(Corr_xy_2_plot_Ad.y_var,'compNonAdjacent');
plt.DATA = Corr_xy_2_plot_Ad{ind,{'res_x','res_y'}};

plt.Colors = plt.CM.Y7_Default;
plt.addScatter([],[])
plt.addShadeErrorForGLM()
xlabel({'Peri-ripple HFB power (dB)';'(adjusted)'})
ylabel({'Non-adjacent inference (adjusted)'})
title('DMN','FontWeight','normal')
iv = find(contains(resultsAdditional.x_var,'Default'),1);
plt.addPstar(sprintf('r = %.2f',resultsAdditional.cor(iv)),resultsAdditional.p_value(iv))
plt.setFont



ax = axes('Position',[0.8 0.7 w w]);
ind = ismember(Corr_xy_2_plot_Ad.x_var,'DMNmPFC') & ismember(Corr_xy_2_plot_Ad.y_var,'compNonAdjacent');
plt.DATA = Corr_xy_2_plot_Ad{ind,{'res_x','res_y'}};

plt.Colors = subcolors(end-2,:);
plt.addScatter([],[])
plt.addShadeErrorForGLM()
xlabel({'Peri-ripple HFB power (dB)';'(adjusted)'})
ylabel({'Non-adjacent inference (adjusted)'})
title('mPFC','FontWeight','normal')
iv = find(contains(resultsAdditional.x_var,'DMNmPFC'),1);
plt.addPstar(sprintf('r = %.2f',resultsAdditional.cor(iv)),resultsAdditional.p_value(iv))
plt.setFont

plt.savefig('Figures/supp_bin_HFB_significance.pdf')


%% visualization, supplementary figure 4, 
% scatter plots for correlations in main figure 5
clf
w = 0.18;

posx = [0 linspace(0.05,0.78,3),linspace(0.05,0.78,3) linspace(0.05,0.78,3),linspace(0.05,0.78,3)]
posx(end) = posx(end-1);
posy = [0 ones(1,3)*0.8  ones(1,3)*0.55  ones(1,3)*0.3  ones(1,3)*0.05]
for i = [2:11 13]
    xvar = rd.x_var{i};
    ax = axes('Position',[posx(i) posy(i) w w]);
    ind = ismember(Corr_xy_2_plot.x_var,xvar) & ismember(Corr_xy_2_plot.y_var,'compNonAdjacent');
    plt.DATA = Corr_xy_2_plot{ind,{'res_x','res_y'}};
    if i>7
        plt.Colors = subcolors(i,:);
    else
    plt.Colors = plt.CM.(xvar);
    end

    plt.addScatter([],[]);
    plt.addShadeErrorForGLM();

    xlabel({'Peri-ripple HFB power (dB) (adjusted)'})
    ylabel({'Non-adjacent inference (adjusted)'})
   
    title(plt.f_DMNsubNames(xvar)); 
    plt.addPstar(sprintf('r = %.2f ',rd.cor(i)),rd.p_value(i)*7)
end
plt.FontSize = 6.9;
plt.setFont 
plt.savefig('Figures/supp_Fig_peri_ripple_HFB_scatters.pdf')
%% visualization, supplementary figure 5, 
Ripple_aligned_activity_during_rest_whole_brain_S7

 

 
