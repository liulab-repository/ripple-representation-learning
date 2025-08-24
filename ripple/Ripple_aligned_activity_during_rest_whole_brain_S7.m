
%% prepare data, using S7 atlas
H = figure(17); clf
plt.H = H;
locy = 0.25;

OverallRR = readtable('../data/table/dd_subject_rest_and_overall_ripple_rate.csv');
OverallRR = OverallRR(:,{'subject','RR_Rest'});

atlas = readtable('../data/table/S7_atlas_labels.csv');


atlas = atlas(periRRchan,:);
atlas.hfb = HFB.EEG_ALLEpoch;
atlas.hfbg   = HFB.EEG_ALL_generalAvg -HFB.EEG_ALL_wholeAvg;
atlas(EEG_channel.Hippocampus>0,:)=[];
atlas.hfbbin = mean(atlas.hfb(:,abs(HFB.times)<250),2);

atlas.atlas_S7r = atlas.atlas_S7;
atlas.Y7_Default = contains(atlas.atlas_S7r,'Default')>0;
MPFCsub = atlas.atlas_S7r(atlas.ismPFC>0&(atlas.Y7_Default>0));
 
atlas.atlas_S7r(atlas.ismPFC&(atlas.Y7_Default>0)) = repmat({'Default_mPFC'},numel(MPFCsub),1);
atlas.atlas_S7r(contains(atlas.atlas_S7r,'Background'))={'Backgroud+Medial Wall'};
atlas.atlas_S7r = regexprep(atlas.atlas_S7r,{'7Networks_LH_' '7Networks_RH_'}, ...
    {'' ''});
atlas.atlas_S7r = regexprep(atlas.atlas_S7r,'_\d+','');

atlas.roi = atlas.atlas_S7r;
 
atlas(contains(atlas.channel,'CREF'),:)=[];

NetOrder = {'Backgroud+Medial Wall' 'Vis' 'SomMot' 'DorsAttn' 'SalVentAttn' 'Limbic' 'Cont' 'Default'};
Ncolors  = plt.CM(1,{'gray','Y7_Visual','Y7_Somatomotor','Y7_DorsalAttention','Y7_VentralAttention','Y7_Limbic','Y7_Frontoparietal','Y7_Default'});
NetNames = cell2table({'' '' '' 'DorsAttn' 'SalVentAttn' 'Limbic' 'Cont' 'Default'},'VariableNames',NetOrder);
NetNames.("Backgroud+Medial Wall")={''};
NetNames.("Vis")={''};
NetNames.("SomMot")={''};
Ncolors.Properties.VariableNames = NetOrder;

atlas_sub = groupsummary(atlas,{'roi' 'subject'},'mean',{'hfbbin' 'hfb' 'hfbg'});
atlas_sub = innerjoin(atlas_sub,subjHfbBehav(:,{'subject','compAdjacent','compNonAdjacent','testAfter'}),"Keys",'subject');
atlas_sum = groupsummary(atlas,'roi','mean',{'hfbbin' 'hfb' 'hfbg'});

  
mhfb = groupsummary(atlas_sum,'roi','mean',{'mean_hfbbin' 'mean_hfb' 'mean_hfbg'});
mhfb.net = regexprep(mhfb.roi,'_.*','');
ord = [];
for ii =1 :numel(NetOrder)
    ind = ismember(mhfb.net,NetOrder{ii});
    ord(ind)=ii;
end
mhfb.order = ord';
mhfb = sortrows(mhfb,'order','ascend');

mhfb.ylabel = mhfb.roi;
mhfb.hfbbin = mhfb.mean_mean_hfbbin;
mhfb.hfb    = mhfb.mean_mean_hfb;
mhfb.hfbg    = mhfb.mean_mean_hfbg;
mhfb = sortrows(mhfb,{'order','net','hfbbin'},{'ascend' 'descend' 'ascend'});


ax1 = axes('Position',[0.05 locy  0.3 0.7]);
dd = mhfb.hfb;
ylb = mhfb.ylabel;
yy = 1:size(dd,1);
imagesc(dd,XData=HFB.times)
yticks(yy),
yticklabels ''
xlim([-1000 1000])
hold on,plot([0 0],get(gca,'YLim'))
xlabel('Time to ripple peak (ms)')
set(gca,'YDir','normal')
cm = flipud(brewermap(64,'RdYlBu'));
colormap(cm)
clim([-0.1  0.15 ])
cbar = colorbar('Location','north','position',[0.1 locy+0.703  0.2 0.013], ...
    'AxisLocation','out');

title(cbar,'Ripple aligned 60-160 Hz power')
set(ax1,'YAxisLocation','right','FontSize',6)
  
ax = axes('Position',[0.45 locy  0.1 0.7]);
for ii = 1:numel(ylb)
text(0,yy(ii),1,strrep(ylb{ii},'_',' '),"HorizontalAlignment","center", ...
        'Rotation',0,'fontsize',6,'FontName','arial')
end
yticks(yy),ylim(yy([1 end])+[-0.5 0.5]);
axis off


% add herizontal bar for peri-ripple HFB
ax2 = axes('Position',[0.585 locy  0.1 0.7]);
dd = mhfb.hfbbin;% mean(dd(:,abs(HFB.times)<250),2);
yy = 1:size(dd,1);
barh(yy,dd,'grouped','facecolor',[1 1 1])
yticks(yy)
yticklabels('')
ylim(yy([1 end])+[-0.5 0.5])
xlabel('HFB power')

% add lobe label
ax3 = axes('Position',[0.52  locy  0.1 0.7]);
hold on
nt = regexprep(mhfb.roi,'_.*','');
[g,gind]=findgroups(nt);

for ii =1:numel(gind)
    yrg = [find(g==ii,1,"first")-0.5,find(g==ii,1,"last")+0.5];
    text(0,mean(yrg),1,strrep(NetNames{1,gind{ii}},'Y7_','') ,"HorizontalAlignment","center", ...
        'Rotation',-90,'fontsize',6)
    plot([0 0],yrg,'-','LineWidth',10,'Color',Ncolors.(gind{ii}))
end
ylim(yy([1 end])+[-0.5 0.5])
xlim([-0.5 0.5])
axis off

% compute correlation
rr=[];pp=[];nn=[];NN=[]; nthresh = 5;
for iroi = 1:numel(ylb)
    
    ind = ismember(atlas_sub.roi,ylb{iroi});
    tmp = atlas_sub(ind,:);
    nn(iroi,:) = size(tmp,1);
    tmp = groupsummary(tmp,'subject',@rbmean,["mean_hfbbin" "mean_hfbg" "compAdjacent" "compNonAdjacent" "testAfter"]);
    tmp.Properties.VariableNames = strrep(tmp.Properties.VariableNames,'fun1_','');
    tmp = innerjoin(tmp,OverallRR,"Keys",'subject');
    xyz = tmp{:,["mean_hfbbin" "compNonAdjacent" "testAfter" ]};% "mean_hfbg" "RR_Rest"
    xyz(sum(isnan(xyz),2)>0,:)=[];
    if size(xyz,1)>nthresh
        [r1,p1]=robustPartialCorrelation(xyz(:,1),xyz(:,2),xyz(:,3:end));   
    else
        r1 = nan;p1=nan;
    end
    rr(iroi,:)=r1;
    pp(iroi,:)=p1;
    NN(iroi,:)=size(xyz,1);
end
zz = norminv(1-pp/2);
zz(rr<0)=-zz(rr<0);
rrpp = table(rr,pp,ylb);
pp(pp>0.05) =nan;
rrppM =table(rr,pp,ylb,NN);
writetable(rrppM,'Results/RRPPM_s7.csv')


% add counts
ax2 = axes('Position',[0.7 locy  0.1 0.7]);
dd = nn; 
yy = 1:size(dd,1);
barh(yy,dd,'grouped','facecolor',[1 1 1]*0.95)
yticks(yy)
yticklabels('')
ylim(yy([1 end])+[-0.5 0.5])
hold on
tmpy = yy; tmpy(1)=tmpy(1)-0.5;tmpy(end)=tmpy(end)+0.5;
plot(nthresh*ones(size(yy)),tmpy,':k')
xlabel('Number of subjects')


% add correlation bar
ax2 = axes('Position',[0.815  locy  0.15 0.7]);
dd = zz(:,1);
yy = 1:size(dd,1);
hb=barh(yy,dd,'facecolor','flat');
yticks(yy),
yticklabels('')
ylim(yy([1 end])+[-0.5 0.5])
xlim([-2 4])
hb.CData  = repmat(plt.CM.FnonadjacentInfer,numel(yy),1);
hb.CData(~isnan(pp(:,1)),:) = repmat(plt.CM.FnonadjacentInfer*0.6,sum(~isnan(pp(:,1))),1);
xlabel('Correlation (z-score)')
pstar = f_pValue2flag(pp(~isnan(pp(:,1)),1));
text(3.5,find(~isnan(pp(:,1)))-0.45,pstar,'Rotation',90 )
 

ind = ismember(atlas_sub.roi,'Default_mPFC' );
tmp = atlas_sub(ind,:);
tmp = groupsummary(tmp,'subject',@rbmean,["mean_hfbbin" "mean_hfbg" "compAdjacent" "compNonAdjacent" "testAfter" "mean_hfb"]);
tmp.Properties.VariableNames = strrep(tmp.Properties.VariableNames,'fun1_','');

xyz = tmp{:,["mean_hfbbin" "compNonAdjacent" "testAfter" "mean_hfbg" ]};
xyz(sum(isnan(xyz),2)>0,:)=[];

ax = axes('Position',[0.82  0.05 0.15 0.15]);
plt.DATA = xyz; plt.Colors = plt.CM.FnonadjacentInfer;
plt.addPartialRegressionScatter
plt.addPstar(sprintf('r = %.3f',rrppM.rr(27)),rrppM.pp(27))
xlabel({'Peri-ripple 60~160 Hz power (dB)';'(adjusted)'})
ylabel({'Non-adjacent inference';'(adjusted)'})

ax = axes('Position',[0.585  0.05 0.15 0.15]);
d = tmp{:,'mean_hfb'}; d = zscore(d,1,2);
m = mean(d,1); se = plt.f_sem(d);
plt.DATA =  d; plt.Colors = [190 161 165] ;
plt.addShadeErrorBar(HFB.times,m,se,{'color',plt.Colors/255})
hold on
plot([0 0],get(gca,'ylim'),':k')
xlabel('Time to ripple peak (s)')
ylabel('60~160 Hz power (z-score) ')

plt.FontSize=7;
plt.setFont
plt.savefig('Figures/supp_peri-ripple-HFB-S5.pdf')
%%
 
%%
CH = load('../data/mat/Contact_distribution_S7_DMNmPFC.mat')
%%
figure(CH.H)
EEG_ALLEpoch_500 = mean(HFB.EEG_ALLEpoch(:,HFB.times>-250&HFB.times<250),2)
channelid = arrayfun(@(x) num2str(x),1:size(channel,1),'un',0)';
channelid  = channelid(periRRchan);
mpfcchannelid = atlas.channel(contains(atlas.atlas_S7r,'Default_mPFC'));
mpfcchannelid = arrayfun(@num2str,find(ismember(channel.channel,mpfcchannelid )),'un',0)

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
        'MarkerEdgeColor',cols(ic,:)*1,"MarkerSize",3)
    
end
disp done!

[~,ia]=setdiff(CH.elecID,mpfcchannelid)
delete(CH.elecH(ia))
exportgraphics(CH.H,'Figures/Peri-ripple-HFB-S7-mPFC.png','Resolution',600)