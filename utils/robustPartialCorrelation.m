function [r,p,stat]=robustPartialCorrelation(x,y,z)
% author, zhibing xiao.
wkdir = pwd;
if ~istable(z)
    z = array2table(z);
end
zvars = z.Properties.VariableNames;
zvars = strjoin(zvars,'+');
tb = [table(x,y), z];
tb(sum(isnan(tb.Variables),2)>0,:)=[];
writetable(tb,'tmp_R.csv')

s = r_script(wkdir,zvars);
RUN_R_Script(s,'tmp_run.R')

stat = readtable('tmp_Routput.csv');
stat = table2struct(stat);
residual_x = readtable('tmp_res_x.csv');
residual_y = readtable('tmp_res_y.csv');
stat.residual_x = residual_x.x;
stat.residual_y = residual_y.x;
r = stat.cor;
p = stat.p_value;
if p<0.001,flag='***';elseif p<0.01,flag='**';elseif p<0.05,flag='*';else,flag='';end
fprintf('n=% 5d  ||  %s = %7.3f, p = %.3f  %s\n',stat.df+1,'p',r,p,flag)

delete('tmp_R.csv')
delete("tmp_Routput.csv")
delete("tmp_res_x.csv")
delete("tmp_res_y.csv")

delete('tmp_run.R')


function RUN_R_Script(R_script,scriptFile)
fid = fopen(scriptFile,'w');
for i = 1:numel(R_script)
    fprintf(fid,'%s\n',R_script{i});
end
fclose(fid);

status = system(['R -f "' scriptFile '"  > /dev/null 2>&1']);
if status>0
    status = system(['R -f "' scriptFile '"']);
    error('R scripts run failed')
end



function s = r_script(wkdir,zvars)
s = {'library(robustbase)'
'library(ppcor)'
'WkDir = ''WORKDIRECTORY'''
'setwd(WkDir)'
'temp_data = read.table(''tmp_R.csv'',header=TRUE,sep=",")'
'robust_model_x <- lmrob(x ~ ZVARS, data = temp_data,   setting = "KS2014")'
'robust_model_y <- lmrob(y ~ ZVARS, data = temp_data,   setting = "KS2014")'
'res_x <- residuals(robust_model_x)'
'res_y <- residuals(robust_model_y)'
'partial_r <- cor.test(res_x, res_y)'
'REPORT = data.frame('
'  t_value = partial_r$statistic,'
'  df = partial_r$parameter,'
'  cor = partial_r$estimate,'
'  p_value = partial_r$p.value'
')'
'write.csv(res_x,''tmp_res_x.csv'')'
'write.csv(res_y,''tmp_res_y.csv'')'
'write.csv(REPORT,''tmp_Routput.csv'')'
};
s = strrep(s,'WORKDIRECTORY',wkdir);
s = strrep(s,'ZVARS', zvars);
