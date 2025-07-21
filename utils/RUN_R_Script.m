function RUN_R_Script(R_script,scriptFile)
% Author, Zhibing Xiao
% update # add option, when R_script is empty, run scriptFile directly
 
disp('R: running...')

if ~isempty(R_script)
    fid = fopen(scriptFile,'w');
    for i = 1:numel(R_script)
        fprintf(fid,'%s\n',R_script{i});
    end
    fclose(fid);
end

scriptDir = char(fileparts(scriptFile));
% status = system(['R -f "' char(scriptFile) '"  ']);
status = system(['cd "' scriptDir '";' ...   
            'R -f "' char(scriptFile) '"  ']); % without output
% status = system(['cd "' scriptDir '";' ...   
%             'R -f "' char(scriptFile) '"  > /dev/null 2>&1 ']); % without output
if status>0
    error('R scripts run failed')
end
disp('R: success')


end