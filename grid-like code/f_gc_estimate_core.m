function [paramEsti]= f_gc_estimate_core(neuralSignal,angles,covars,noTest,Nrep)
 
if ~exist('covars','var')
    covars = [];
end
if ~exist('Nrep','var') 
    Nrep = 1; % repetation for stability. Only needed when random split train/test
end

DATA = table(neuralSignal,angles);
%% estmate omega by crossvalidation
for irep = 1:Nrep 
     
    N_trial = numel(angles);
    ind = ones(N_trial ,1);
    ind(2:2:end)=2;


    trainData = DATA(ind==1,:);
    testData  = DATA(ind==2,:);
    COV = struct;
    if ~isempty(covars)
        COV.train = covars(ind==1,:);
        COV.test  = covars(ind==2,:);
    else
        COV = [];
    end
    
    paramEsti  = estmate(trainData ,testData,COV);
 
    if ~isempty(covars)
        cov_new = struct; cov_new.train = COV.test; cov_new.test = COV.train;
    else
        cov_new =[];
    end
    paramEsti2 = estmate(testData ,trainData,cov_new); % switch train test
    paramEsti.beta_estimate  = mean([paramEsti.beta_estimate,paramEsti2.beta_estimate],2);
    paramEsti.beta_estimate_p  = nan(size(paramEsti.beta_estimate_p));
    paramEsti.omega_estimate  = mean([paramEsti.omega_estimate,paramEsti2.omega_estimate],2);
    paramEsti.beta_cov  = mean([paramEsti.beta_cov,paramEsti2.beta_cov],2);

    tmp.beta_estimate(:,:,irep) =   paramEsti.beta_estimate  ;
    tmp.beta_estimate_p(:,:,irep) =   paramEsti.beta_estimate_p  ;
    tmp.omega_estimate(:,:,irep) =   paramEsti.omega_estimate  ;
    tmp.beta_cov(:,:,irep) =   paramEsti.beta_cov  ;

end
paramEsti.beta_estimate  = mean(tmp.beta_estimate,3);
paramEsti.beta_estimate_p  = mean(tmp.beta_estimate_p,3);
paramEsti.omega_estimate  = mean(tmp.omega_estimate,3);
paramEsti.beta_cov    = mean(tmp.beta_cov,3);

omega_estimate = paramEsti.omega_estimate(6);
 
% use all data to compute bin activation, just for visualization
testData = DATA; %[trainData;testData]; 
 
%% align vs misalign
if exist('noTest','var') 
    if logical(noTest)
        return
    end
end

alignedD_360 = mod(testData.angles - omega_estimate,360);
anglebinNum = round(alignedD_360/30)+1;
anglebinNum(anglebinNum==13)=1;
if numel(unique(anglebinNum))<4,error('angles are too sparse'),end
 
labels=cell(size(anglebinNum));
binind = mod(anglebinNum,2)~=0;
labels(binind) ={'align'};
labels(~binind)={'misalign'};
 
labels = categorical(labels) ;
y_test = testData.neuralSignal;
dplt = [table(y_test ),table(labels),table(anglebinNum)];
 
% yg1 = dplt.y_test(dplt.labels=='align');
% yg2 = dplt.y_test(dplt.labels=='misalign');
% [~,p,~,stat] = ttest2(yg1, yg2 );
% T = stat.tstat;
% n1 = numel(yg1);
% n2 = numel(yg2);
% s1 = var(yg1);
% s2 = var(yg2);
% m1 = mean(yg1);
% m2 = mean(yg2);

% 
% sP=((n1-1).*s1 + (n2-1).*s2)./(n1+n2-2); % pooled (within-groups) variance
% es=(m1-m2)./sqrt(sP);
binAvg = grpstats(dplt(:,[1 3]),'anglebinNum','mean');
binAvg.GroupCount = [];
binAvgActivation = nan(12,1);
binAvgActivation(binAvg.anglebinNum)=binAvg.mean_y_test;
paramEsti.binAvgActivation= binAvgActivation;

function paramEsti = estmate(trainData ,testData, covars)
if isempty(covars)
    covars_train = [];
    covars_test  = [];
else
    covars_train = covars.train;
    covars_test  = covars.test;    
end
nantmp = nan(12,1);

paramEsti= struct;
paramEsti.beta_estimate = nantmp;
paramEsti.beta_estimate_p = nantmp;
paramEsti.omega_estimate = nantmp;
paramEsti.beta_cov = nantmp;

for ifold = 4:8
    x0 = ones(size(trainData,1),1);
    x1 = sin(deg2rad(ifold*trainData.angles));
    x2 = cos(deg2rad(ifold*trainData.angles));
    X =  [x0 x1 x2 covars_train]; % add constant and covariable
    betaTrain = pinv(X)*trainData.neuralSignal;
 
    betaTrain = betaTrain(2:3);
    omega_estimate = rad2deg(atan2(betaTrain(1),betaTrain(2)))/ifold;

    % test the omega effect
    x = cos(deg2rad(ifold*(testData.angles - omega_estimate)));
    xx = [x covars_test];
    [bata , ~, stats ]=  f_pinv_GLM(xx,testData.neuralSignal);    
    beta_p = stats.p(2);
    beta_estimate = bata(2);
    paramEsti.beta_estimate(ifold,1) = beta_estimate;
    paramEsti.beta_estimate_p(ifold,1)= beta_p;
    paramEsti.omega_estimate(ifold,1) = omega_estimate;
    if ~isempty(covars_test)
        paramEsti.beta_cov(ifold,1) = bata(3);
    end
 
end

 