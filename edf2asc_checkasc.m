function edf2asc_checkasc(dataPath,exeFilePath)
% This script can convert .EDF files to .asc files.
% The processing might need a period of time. If you have a significant
% number of files need to processing, please well arrange your time 
% before running this script.
%
% By BYC 2018-10-8

if ~exist('dataPath','var')
    error('There is no file in this path or there is a wrong path.');
end

if ~exist('exeFilePath','var')
    error('There is no file in this path or there is a wrong path.');
elseif strfind(exeFilePath,'edf2asc.exe')
    exeFileName = [exeFilePath ' -ntime_check'];
else
    exeFileName = [fullfile(exeFilePath,'edf2asc.exe') ' -ntime_check'];
end

cd (dataPath);
datasfile= dir([dataPath '\*.EDF']);
datasNum = length(datasfile);
status=[];
check_operation = 0; % how to deal the existed asc files? 1. overwrite; 0. skip. You could set 0 after the second process.

if check_operation == 1
    delete *.asc
end

for i = 1: datasNum

       oldFileName = fullfile(dataPath,datasfile(i).name);
       if strfind(oldFileName,'.EDF')
           check_ascName = strrep(datasfile(i).name,'.EDF','.asc');
       elseif strfind(oldFileName,'.edf')
           check_ascName = strrep(datasfile(i).name,'.edf','.asc');
       else
           error([datasfile(i).name ' is not a EDF file']);
       end
       
       if check_operation == 0
           check_result = fopen(check_ascName);
           if check_result > 0
               status = cat(1,status, [-1 i]);
               continue
           end
       end
           
       cmd = [exeFileName 32 oldFileName];
       [~, cmdout]= system(cmd);
       status = cat(1,status, [contains(cmdout,'success') i]);
end

% %% reverse for the error report
% error_report = status(status(:,1)~=1,:);
% if ~isempty(error_report)
%     report_savename = [datestr(now,'yymmddHHMM') '_edf2asc.mat'];
%     save(report_savename,'status','error_report');
% end

end