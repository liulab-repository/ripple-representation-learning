
function Results = fitLME_byR(dataForLME,formula,detailVar, TAG, outputpath)
% Author, Zhibing Xiao @bnu
% automatically report contrast, pMat, lsmeans, and interactions(if need)
% Simple version.
% detailVar = 'TrialType'
if ~exist('TAG','var'),TAG = 'temp';end
if ~exist('outputpath','var'),outputpath = pwd;end
TAG = char(TAG); outputpath = char(outputpath);
outputpath = fullfile(outputpath,TAG,filesep);

if exist(outputpath,'dir'),rmdir(outputpath,"s"),end
mkdir(outputpath);
if ~isequal(outputpath(1),filesep),outputpath = fullfile(pwd,outputpath);end
if ~isequal(outputpath(end),filesep),outputpath = fullfile(outputpath,filesep);end

dfile = sprintf('dataForLME_%s.csv',TAG);
if exist(dfile,'file'),delete(dfile);end

Results = struct;


vars = regexp(formula,'\w+','match');
indUseVar = intersect(dataForLME.Properties.VariableNames,vars);
dataForLME = dataForLME(:,indUseVar);
nanIdx = false(size(dataForLME));
for ivar = 1:width(dataForLME)
    if isnumeric(dataForLME.(ivar))
        nanIdx(:,ivar) = isnan(dataForLME.(ivar));
    end
end
dataForLME(any(nanIdx,2),:) =[];

writetable(dataForLME,fullfile(outputpath,dfile))


if ~exist('detailVar','var')
    detailVar = vars{2};
end



Rscript = R_s_tmplate;
 
if isempty(detailVar)
    ind=find(contains(Rscript,'DETAILVAR'));
    Rscript(ind(end):end,:)=[];
end



if exist('outputpath','var')
    if ~exist('outputpath','dir'),mkdir(outputpath);end
    Rscript = strrep(Rscript,'getwd()',['''' outputpath '''']);
end

scriptFile = fullfile(outputpath,sprintf('%s_mixed_model_script.R',TAG) );

intr_str = regexp(formula,'(?<=[\-\+]?)\w+*\w+(?=[\-\+]?)','match');
if ~isempty(intr_str)
    vs = split(intr_str,'*');
    vsIsFactor = [false false];
    for ivs = 1:2
        if iscell(dataForLME.(vs{ivs})) || iscategorical(dataForLME.(vs{ivs})) || isstring(dataForLME.(vs{ivs})) || all(ismember(dataForLME.(vs{ivs}),[0 1]))
            vsIsFactor(ivs) = true;
        end
    end
    if all(vsIsFactor==1)
        rscript=lme_interaction(vs{1},vs{2});
        interaction_file1 = fullfile(outputpath,['contrast_emmeans_' vs{1} '_by_' vs{2} '_' TAG '.csv']);        
    else
        rscript=lme_interaction_continous(vs{~vsIsFactor},vs{vsIsFactor});
        interaction_file1 = fullfile(outputpath,['interaction_' vs{1} '_by_' vs{2} '_' TAG '_simple_slop.csv']);
    end
    Rscript = [Rscript;rscript];
    % interaction_file2 = fullfile(outputpath,['contrast_emmeans_' vs{2} '_by_' vs{1} '_' TAG '.csv']);
else
    interaction_file1  = '';%interaction_file2 = '';
end

Rscript = strrep(Rscript,'dataForLME.csv',dfile);
Rscript = strrep(Rscript,'FORMULA',formula);
Rscript = strrep(Rscript,'DETAILVAR',detailVar);
Rscript = strrep(Rscript,'TAG',TAG);
 
try
    RUN_R_Script(Rscript, scriptFile )
end
warning off
lme = readtable(fullfile(outputpath,['mixed_effects_model_' TAG '.csv']));
tt = readtable(fullfile(outputpath,['tvalues_' TAG '.csv']));
stats.T = tt;
stats.F = lme;
if isempty(detailVar)
    lsmeans=[]; contrast=[];  pMat=[]; interaction=[];
    warning on

else

    lsmeans  = readtable(fullfile(outputpath,['lsmeans_' detailVar '_' TAG '.csv']));
    contrast = readtable(fullfile(outputpath,['contrast_' detailVar '_' TAG '.csv']));


    if size(contrast,1)>=1&&isnumeric(contrast.estimate)

        uni_cond = unique(dataForLME.(detailVar));
        if isnumeric(uni_cond)
            uni_cond = num2str(uni_cond);
            uni_cond = strcat(detailVar,uni_cond);
        end
        uni_cond = cellstr(uni_cond);


        contrast_pairs = split(contrast.contrast,' - ');
        if size(contrast,1)==1, contrast_pairs = contrast_pairs';end
        if size(contrast_pairs,2)~=2,error('your condition names possibly contains " - "');end

        pMat=nan(size(contrast,1));
        for ii = 1:size(contrast_pairs,1)
            VarPairs = contrast_pairs(ii,:);
            r = ismember(uni_cond,VarPairs{1});
            c = ismember(uni_cond,VarPairs{2});
            pMat(r,c)  = contrast.p_value(ii);
            pMat(c,r)  = contrast.p_value(ii);
        end
        pMat(all(isnan(pMat),2),:) =[];
        pMat(:,all(isnan(pMat),1)) =[];
        pMat = array2table(pMat,"RowNames",uni_cond,'VariableNames',uni_cond);

        if ~isempty(intr_str)
            vs1 = readtable(fullfile(outputpath,['contrast_' vs{1} '_by_' vs{2} '_' TAG '.csv']));
            vs2 = readtable(fullfile(outputpath,['contrast_' vs{2} '_by_' vs{1} '_' TAG '.csv']));            
            interaction.(vs{1}) = vs1;
            interaction.(vs{2}) = vs2;

            vs1.Properties.VariableNames{3} = 'condition';
            vs2.Properties.VariableNames{3} = 'condition';
            vs1.contrastVar = repmat(vs(1),size(vs1,1),1);
            vs2.contrastVar = repmat(vs(2),size(vs2,1),1);
            vs12 = [vs1;vs2];
            vs12 = movevars(vs12,'contrastVar','After','condition');
            interaction.contrast = vs12;
        else
            interaction = {};
        end
    else
        pMat = [];

    end

    if ~isempty(interaction_file1)
        temp1 = readtable(interaction_file1);
        %     temp2=  readtable(interaction_file2);
        %     interaction = [];
        eemeans = temp1; %
        interaction.eemeans = eemeans;
        try
        interaction.pMat = extract_P_value(interaction,vs); 
        interaction = rmfield(interaction,vs{1});
        interaction = rmfield(interaction,vs{2});
        end
        if contains(interaction_file1,'simple_slop')
            Results.simpleSlop   = readtable(interaction_file1);
            Results.simpleEffect = readtable(strrep(interaction_file1,'simple_slop','fit'));
        end
    else
        interaction = [];
    end

end



Results.lme = lme;
Results.lsmeans  = lsmeans;
Results.contrast = contrast;
Results.pMat     = pMat;
Results.interaction = interaction;
Results.T  = stats.T;
Results.F  = stats.F;

warning on

printFields(Results, 'Results')
end
function   printFields(s, name)
    if nargin < 2
        name = 's';  
    end

    if ~isstruct(s)
        fprintf('%s:\n ', name);
        disp(s);
        
    else
        fields = fieldnames(s);
        for i = 1:numel(s)
            for j = 1:length(fields)
                fieldName = fields{j};
                value = s(i).(fieldName);
                fullName = sprintf('%s(%d).%s', name, i, fieldName);
                printFields(value, fullName);   
            end
        end
         
    end
end


function P = extract_P_value(interaction,vs)
try
    ee = interaction.eemeans ; ee.pairs = strcat(ee.(vs{1}),'-',ee.(vs{2}));
catch
    P =[];
    return
end
try
contra{1} = interaction.(vs{1});
contra{1}.Properties.VariableNames(ismember(contra{1}.Properties.VariableNames,vs))={'cond'};
contra{2} = interaction.(vs{2});
contra{2}.Properties.VariableNames(ismember(contra{2}.Properties.VariableNames,vs))={'cond'};
contra = cat(1,contra{:});

for i =1:size(contra,1)
    contra.pairs{i} = strcat(split(contra.contrast{i},' - ')','-',contra.cond{i});
end
P= nan(size(ee.pairs,1));
for i =1:size(ee.pairs,1)
    for j =1:size(ee.pairs,1)
        tar={ee.pairs{i},ee.pairs{j}};
        ind = findPairs(contra.pairs,tar);
        if ~isempty(ind)
            P(i,j) = contra.p_value(ind);
        end
    end
end
pp(:,:,1) = P;pp(:,:,2) = P';
P = squeeze(nanmean(pp,3));
catch
    P = [];
end


end

function ind = findPairs(pairs,tar)
ct = sortrows([split(tar{1},'-'),split(tar{2},'-')]);
ind=[];
for i = 1:numel(pairs)
    cp = sortrows([split(pairs{i}{1},'-'),split(pairs{i}{2},'-')]) ;
    if isequal(cp,ct)||isequal(cp,ct(:,[2 1]))
        ind =i;
    end
end
end
function RUN_R_Script(R_script,scriptFile)
fid = fopen(scriptFile,'w');
for i = 1:numel(R_script)
    fprintf(fid,'%s\n',R_script{i});
end
fclose(fid);

status = system(['R -f "' scriptFile '"']);
if status>0
    error('R scripts run failed')
end
end


function Rscript = R_s_tmplate
Rscript = ...
    {
    '#Antomicalanalysis-mixedeffectsmodel:(thisisatemplate,don''tmodify)'
    '#'
    '#loadrequiredpackages,ifnotinstalled,theninstallthem.'
    'if(!require(''lme4''))'
    '  install.packages("lme4")'
    'if(!require(''lsmeans''))'
    '  install.packages("lsmeans")'
    'if(!require(''afex''))'
    '  install.packages("afex")'
    '#'
    '#library(multcompView)'
    '#library(multcomp)'
    '#'
    '#Loaddata:(selectmanually)'
    'fileName=''dataForLME.csv'''
    'stat="TAG";'
    '#SETTHECORRECTPATHTODATAFILESCREATEDINMATLAB(MANUALLY):'
    '#Definethecorrectpath:'
    'parentfolder=getwd();'
    'data_file<-paste(parentfolder,fileName,sep="")'
    'D<-read.table(data_file,header=TRUE,sep=",")'
    '#'
    '#Outfiles:'
    'path<-paste(parentfolder,sep="")'
    'file1<-paste("mixed_effects_model_",stat,".csv",sep="")'
    'file2<-paste("tvalues_",stat,".csv",sep="")'
    'file3<-paste("lsmeans_DETAILVAR_",stat,".csv",sep="")'
    'file4<-paste("contrast_DETAILVAR_",stat,".csv",sep="")'
    'dir.create(path)'
    'file1=file.path(path,file1)'
    'file2=file.path(path,file2)'
    'file3=file.path(path,file3)'
    'file4=file.path(path,file4)'
    'setwd(path)'
    '#'
    'model<-mixed(FORMULA,data=D,test_intercept=FALSE)'
    'write.csv(model$anova_table,file1)'
    'summary(model$full.model)'
    'anova(model)'
    '#'
    'coef(summary(model))'
    'write.csv(coef(summary(model)),file2)'
    '#'
    'output3<-lsmeans(model,pairwise~DETAILVAR,adjust="None",mode="kenward-roger")'
    'write.csv(summary(output3$lsmeans),file3)'
    'write.csv(summary(output3$contrast),file4)'
    '#=================================='
    };
end

function rscript=lme_interaction(VAR1,VAR2)
rscript_tmplate={'file5<-paste("contrast_VAR1_by_VAR2_",stat,".csv",sep="")'
    'file5=file.path(path,file5)'
    'file6<-paste("contrast_emmeans_VAR1_by_VAR2_",stat,".csv",sep="")'
    'file6=file.path(path,file6)'
    'emm_i1 <- emmeans(model, "VAR1", by = c("VAR2"))'
    'output=update(pairs(emm_i1), by = NULL, adjust = "holm")'
    'write.csv(summary(output),file5)'
    'write.csv(summary(emm_i1),file6)'
    };

rscript1 = strrep(rscript_tmplate,'VAR1',VAR1);
rscript1 = strrep(rscript1,'VAR2',VAR2);

rscript2 = strrep(rscript_tmplate,'VAR1',VAR2);
rscript2 = strrep(rscript2,'VAR2',VAR1);

rscript = [rscript1;rscript2];

end
 

function rscript=lme_interaction_continous(varCont,varFactor)
rscript_tmplate={
'library(interactions)'
'library(effects)'
'# mixed model'
'model <- lmer(FORMULA , D)'
'write.csv(coef(summary(model)), paste0(''mixed_model_lmer_'', ''TAG'', ''.csv''))'
'print(paste("Model coefficients for", "TAG", ":"))'
'print(coef(summary(model)))'
'# interaction'
'sim_slopes_results <- sim_slopes(model, pred = "varCont", modx = "varFactor", johnson_neyman = FALSE, cond.int = TRUE)'
'write.csv(sim_slopes_results$ints, paste0(''interaction_varCont_varFactor_int.csv''))'
'write.csv(sim_slopes_results$slopes, paste0(''interaction_varCont_varFactor_slope.csv''))'
't1 <- sim_slopes_results$ints'
't2 <- sim_slopes_results$slopes'
't1$type <- "intercept"'
't2$type <- "slop"'
'combined_table <- rbind(t1, t2)'
'write.csv(combined_table, paste0(''interaction_varCont_by_varFactor_TAG_simple_slop.csv''))'
'print(paste("Simulated slopes for", "TAG", ":"))'
'print(sim_slopes_results$slopes)'
'# effect'
'min_value <- min(D$varCont)'
'max_value <- max(D$varCont)'
'npoint <- 20'
'generated_values <- seq(min_value, max_value, length.out = npoint)'
'inter.sd <- effect(c(''varCont*varFactor''), mod = model, xlevels = list(varFactor = c(1, 0), varCont = generated_values ))'
'inter.sd <- as.data.frame(inter.sd)'
'write.csv(inter.sd, paste0(''interaction_varCont_by_varFactor_TAG_fit.csv''))'};

% 'file6<-paste("contrast_emmeans_VAR1_by_VAR2_",stat,".csv",sep="")'

rscript1 = strrep(rscript_tmplate,'varCont',varCont);
rscript1 = strrep(rscript1,'varFactor',varFactor);
 
rscript = [rscript1;];

end





