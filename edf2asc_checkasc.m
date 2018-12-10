function edf2asc_checkasc(dataPath,exeFilePath)
% This script can convert .EDF files to .asc files.
% The processing might need a period of time. If you have a significant
% number of files need processing, please well arrange your time 
% before running this script.
%
% By BYC 2018-10-8

if ~exist('dataPath','var')
    error('Please input a data path.');
end

if ~exist('exeFilePath','var')
    error('Please input the path for edf2asc.exe.');
% elseif strfind(exeFilePath,'edf2asc.exe')
%     exeFileName = [exeFilePath ' -ntime_check'];
else
    exeFileName = [fullfile(exeFilePath,'edf2asc.exe') ' -ntime_check -y -vel -miss 0 -failsafe'];
end

cd (dataPath);
datasfile1= dir([dataPath '\*.EDF']);
datasfile2 = dir([dataPath '\*.edf']);

if length(datasfile1) > length(datasfile2)
    datasfile = datasfile1;
else
    datasfile = datasfile2;
end

datasNum = length(datasfile);
status=[];

for i = 1: datasNum

       oldFileName = fullfile(dataPath,datasfile(i).name);
       if strfind(oldFileName,'.EDF')
           check_ascName = strrep(datasfile(i).name,'.EDF','.asc');
       elseif strfind(oldFileName,'.edf')
           check_ascName = strrep(datasfile(i).name,'.edf','.asc');
       else
           error([datasfile(i).name ' is not a EDF file']);
       end
           
       cmd = [exeFileName 32 oldFileName];
       [~, cmdout]= system(cmd);
       status = cat(1,status, [strfind(cmdout,'success') i]);
end

% %% reverse for the error report
% error_report = status(status(:,1)~=1,:);
% if ~isempty(error_report)
%     report_savename = [datestr(now,'yymmddHHMM') '_edf2asc.mat'];
%     save(report_savename,'status','error_report');
% end

end