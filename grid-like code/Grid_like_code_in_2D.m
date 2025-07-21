clear,clc
addpath('../utils/')


% load data
D = load('../data/Data_power (3-7Hz).mat');
channels = readtable('../data/Channel_ROI.csv');
load('../data/Allangles.mat')
ROIs =channels.Properties.VariableNames(3:end);
 
%% compute grid like code

% compute grid-like code for each time point during mental simulation
[Point_BETAs,Point_BETAs_p, Point_OMEGA,Point_binAvgActivation, Point_SixFold]= ...
    f_estimate_gridCode_byChannelAndTime(D.ALLSignals,Allangles,1);

% compute grid-like code for averaged activity during mental simulation
AllSignalMean = cellfun(@(x) mean(x(1,D.time>0,:)),D.ALLSignals,'UniformOutput',0);

[Mean_BETAs,Mean_BETAs_p, Mean_OMEGA,Mean_binAvgActivation, Mean_SixFold]= ...
    f_estimate_gridCode_byChannelAndTime(AllSignalMean,Allangles,1);

time = D.time;
save('Results/grid_code.mat','Point_SixFold','Mean_SixFold','Mean_BETAs','Point_BETAs','time','-nocompression')
%%
% run ../ripple/  for statistic analysis and visualization
%%
% report = table('size',[0,4],'VariableTypes',{'string','double','double','double'}, ...
%     'VariableNames',{'ROI','df','t','p'});
% for ROI = channels.Properties.VariableNames(3:end)
%     % ROI = {'Entorhinal'};
%     b = Mean_SixFold.beta(channels.(ROI{1})>0);
%     [~,rmind1]=rmoutliers(b,"mean","ThresholdFactor",3); % 3*std
%     b(rmind1,:)=[];
%     [h,p,ci,st]=ttest(b);
%     t = st.tstat; df = st.df;
% 
%     report = cat(1,report,table(ROI,df,t,p));
% end
% report.h = report.p<0.05;
% report
% 
% cpath = path;
% restoredefaultpath
% %%
% 
% addpath('/HDD/lbcode-dev-xzb/externalPackages/fieldtrip-20220714')
% ft_defaults
% addpath(genpath('/HDD/lbcode-dev-xzb/externalPackages/eeglab2022.1/functions/'))
% addpath('../utils/')
% time =D.time;
% trng = [min(time),max(time)];
% figure(1),clf
% t =tiledlayout(8,2);
% for ROI = ROIs
%     nexttile
%     braw = Point_SixFold.beta(channels.(ROI{1})>0,:);
%     b = braw-mean(braw(:,time<0),2);
%     se = std(b)/sqrt(size(b,1));
%     m = mean(b,1);
%     shadedErrorBar(time,m,se)
%     xlim(trng)
%     hold on
%     plot(trng,[0,0],'k--')
%     ylm = get(gca,'YLim');
%     plot([0,0],ylm,'k--')
%     title(ROI{1})
% 
%     bs = b;
%     bs = movmean(bs,10,2);
% 
%     if isequal(ROI{1},'Entorhinal')
%         erps = {bs'; bs'*0};
%         M_group = {'groupstats', 'off'};
%         M_cond = {'condstats', 'on', 'paired',{'on'}};
%         M_state = [M_cond,M_group];
% 
%         [pcond, pgroup, pinter, statscond, statsgroup, statsinter] = std_stat(erps, ...
%             M_state{:} , ...
%             'fieldtripmcorrect', 'cluster', ...
%             'fieldtripmethod', 'montecarlo', ...
%             'mode', 'fieldtrip');
%         adj_p = pcond{1};
%     else
%         [~,p]=ttest(bs);
%         [~,~,~,adj_p]=fdr_bh(p,0.05,'pdep');
%     end
% 
%     T_mask = adj_p<0.05;
%     mt =time;   mt(~T_mask)=nan;
%     plot(mt,ylm(2)*ones(size(time)),'-','LineWidth',5)
% 
% end 
% 
% % restoredefaultpath
% % addpath(cpath)
% 
% %%
% 
% 
% GCmeanByT = table;
% GCmeanByT.subject = unique(channels.subject);
% GG_t_sub = GCmeanByT;
% for ROI = ROIs
% 
%     chanind = logical(channels.(ROI{1}));
%     subject = channels.subject(chanind);
%     beta    = Mean_SixFold.beta(chanind);
%     betaByT  = Point_SixFold.beta(chanind,:);
% 
% 
%     betaByTmean = mean(betaByT(:,time>0),2); 
% 
%     [betaByTmean, gname ]= groupsummary(betaByTmean,subject,@(x)robustMean(x,1));
%     [betamean, ~ ]= groupsummary(beta,subject,@(x)robustMean(x,1));
% 
%     t = table;
%     t.subject = gname;
%     t.(['gcbyTmean_' ROI{1}]) = betaByTmean;
%     t.(['gcmean_' ROI{1}]) = betamean;
% 
%     GCmeanByT = outerjoin(GCmeanByT,t,"Keys","subject");
%     GCmeanByT.subject = GCmeanByT.subject_GCmeanByT;
%     GCmeanByT(:,contains(GCmeanByT.Properties.VariableNames,'subject_'))=[];
% 
%     [betam, gname2 ] = groupsummary(betaByT,subject,@(x)mean(x,1));
%     tt = table; tt.subject = gname2;
%     tt.(['gcbyT_' ROI{1}]) = betam;
% 
%     GG_t_sub = outerjoin(GG_t_sub,tt,"Keys","subject");
%     GG_t_sub.subject = GG_t_sub.subject_GG_t_sub;
%     GG_t_sub(:,contains(GG_t_sub.Properties.VariableNames,'subject_'))=[];
% 
% 
% 
% end 
% 
% mkdir('Results')
% writetable(GCmeanByT,'Results/dt_subject_grid_code.csv')
% %%
% dtripple = readtable('../data/dd_subject_rest_and_overall_ripple_rate.csv');
% dtacc    = readtable('../data/dd_subject_performance.csv');
% dtsubinfo = readtable('../data/dd_subject_info.csv');
% dtacc.Properties.VariableNames = strcat('ACC_',dtacc.Properties.VariableNames);
% dtacc.subject = dtacc.ACC_subject;
% dt = innerjoin(GCmeanByT,dtripple,"Keys","subject");
% dt = innerjoin(dt,dtacc,'Keys','subject');
% dt = innerjoin(dt,GG_t_sub,'Keys','subject');
% % dt = innerjoin(dt,dtsubinfo,'Keys', 'subject');
% clc
% for ROI = ROIs
%     fprintf('ROI: %s ',ROI{1})
%     partialcorrWithNan(dt.RR_Rest, dt.(['gcbyTmean_' ROI{1}]),dt.RR_Overall);
%     % gcmean_ ,gcbyTmean_
% end
% % [r,p]=
% partialcorrWithNan(dt.ACC_testAfterInfer,dt.("gcbyTmean_Entorhinal"),dt.ACC_testAfterMemory)
% 
% % [r,p]=partialcorrWithNan(dt.ACC_compNonAdjacent,dt.("gcbyTmean_Entorhinal"),dt.RR_Overall)
% 
% % [r,p]=corrWithNan(dt.ACC_compNonAdjacent,dt.("gcbyTmean_Entorhinal"))
% %%
% for ROI = ROIs
%     z = [dt.gender, dt.age];
%     fprintf('% 20s: ',ROI{1})
%     partialcorrWithNan(dt.ACC_compAll, dt.(['gcbyTmean_' ROI{1}]),z);
% %     robustPartialCorrelation(dt.ACC_compNonAdjacent,dt.(['gcmean_' ROI{1}]),z);
% 
% end
% % for ROI = ROIs
% %     corrWithNan(dt.ACC_compNonAdjacent, dt.(['gcmean_' ROI{1}]));
% %     corrWithNan(dt.ACC_compAdjacent, dt.(['gcmean_' ROI{1}]));
% % end
% 
% %%
% clc
% 
% figure(3),clf
% t =tiledlayout(8,2);
% for ROI = ROIs
%     nexttile
%     braw =  dt.(['gcbyT_' ROI{1}]);
%     braw = braw - mean(braw(:,time<0),2);
% 
%     [r,p]=corrWithNan(dt.ACC_compNonAdjacent, braw);
% 
% 
%     plot(time,r)
% 
% 
%     [~,~,~,adj_p]=fdr_bh(p,0.05,'pdep');
% 
%     ylim([-1 1]); xlim(time([1 end]))
%     hold on
%     T_mask = adj_p<0.05;
%     mt =time;   mt(~T_mask)=nan;
%     plot(mt,ylm(2)*ones(size(time)),'-','LineWidth',5)
%     title(ROI{1})
% end
