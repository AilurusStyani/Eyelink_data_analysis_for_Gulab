fileNamei = 'D:\ION\2018_rotation\ZXY\data\eyedatas_VT_B_test1810171809.asc';

fid=fopen(fileNamei);
fseek(fid,0,'eof');
numline=ftell(fid);
fclose(fid);
rawdata=importdata(fileNamei,' ',numline);
%extrace data

ind=strfind(rawdata,'...');
ind=~cellfun(@isempty,ind);
data=rawdata(ind);
clear rawdata;
data=strrep(data,'...','');
data=cellfun(@str2num,data, 'UniformOutput',false);
positionData=cell2mat(data);

eyeDatas = positionData(:,1:4);
clear data;
clear positionData;
clear ind;

figure(1000)
plot(eyeDatas(:,2),eyeDatas(:,3))