function xm = temporalRobmean(x,avgDIM,varDim)
    if nargin<2
        avgDIM = 1;
        varDim = 2;
    end
    if nargin<3
        varDim = 2;
    end
    sz = size(x);
    szm = ones(size(sz));
    szm(varDim)=sz(varDim);

    % variation
    sd   =  std(x,1,varDim); % temporal variation
    sdsd =  std(sd,1,avgDIM);  % cross epochs
    sdmean = mean(sd,avgDIM); 
    extremeind = sd> sdmean+3*sdsd;    
    exind1 = repmat(extremeind,szm);
    

    % amplitude level
    mm   =  mean(x,varDim); % temporal 
    sdmm =  std(mm,1,avgDIM);  % cross epochs
    mmmean = mean(mm,avgDIM); 
    extremeind = (mm > mmmean+3*sdmm)| (mm<mmmean-3*sdmm);        
    exind2 = repmat(extremeind,szm);

    x(exind1|exind2) = NaN;
    
   
    xm = nanmean(x,avgDIM);

    


