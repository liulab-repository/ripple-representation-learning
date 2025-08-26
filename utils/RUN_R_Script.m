function RUN_R_Script(R_script,scriptFile)
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
status = system(['cd "' scriptDir '";' ...   
            'R -f "' char(scriptFile) '"  ']); % without output
if status>0
    error('R scripts run failed')
end
disp('R: success')


end