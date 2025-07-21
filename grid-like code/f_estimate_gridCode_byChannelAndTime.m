function [BETAs,BETAs_p, OMEGA,binAvgActivation, SixFold] = ...
    f_estimate_gridCode_byChannelAndTime(...
    Signals,Allangles,Nrep)
 

if ~exist('Nrep','var'),Nrep=1;end
 
[nbchan, ~, ~] = size(Signals);
point = size(Signals{1},2);

nanResults = nan([nbchan point 12]);

BETAs = nanResults;
OMEGA = nanResults;
BETAs_p = nanResults;
BETAs_cov = nanResults;

binAvgActivation = BETAs;
 
tic
TT  = nan(nbchan,point);ES  = TT;T_p = TT;
warning  off
parfor ichan = 1:nbchan

    fprintf('.')
    angles = Allangles{ichan};
    angles = mod(angles,360); % scale angles to [0 360]
    csignal = Signals{ichan};

    for p = 1 :point
        neuralSignal = squeeze(csignal(1,p,:));
        [pEst ]= f_gc_estimate_core(neuralSignal,angles,[],false,Nrep);
 
        BETAs_p(ichan,p,:) = pEst.beta_estimate_p;
        BETAs(ichan,p,:) = pEst.beta_estimate;
        OMEGA(ichan,p,:) = pEst.omega_estimate;
        binAvgActivation(ichan,p,:) = pEst.binAvgActivation;
        % T(ichan,p) = t;
        % T_p(ichan,p) = tp;
        % ES(ichan,p) = es;
        BETAs_cov(ichan,p,:) = pEst.beta_cov;
    end
end



SixFold = struct;
SixFold.beta = squeeze(BETAs(:,:,6));
SixFold.beta_p = squeeze(BETAs_p(:,:,6));
SixFold.omaga  = squeeze(OMEGA(:,:,6));
SixFold.beta_cov = squeeze(BETAs_cov(:,:,6));
% SixFold.T = T;
% SixFold.T_p = T_p;
% SixFold.es = ES;
SixFold.binActivation  = binAvgActivation;

 
toc
fprintf('\n')
warning on

 

