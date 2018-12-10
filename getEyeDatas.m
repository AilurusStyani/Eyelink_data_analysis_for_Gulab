function getEyeDatas(dataPath,error_set,pupil_baseline_duration)

data_check = 0; % 1. check data validation, poor block will be skipped; otherwise all data will be processed.

if ~exist('dataPath','var') || isempty(dir(fullfile(dataPath,'*.asc')))
    error('There is no asc files in this path or there is a wrong path.');
else
    filespath = dataPath;
end

status = [];
cd(filespath);

% find all the asc files in this folder
fileName = dir('*.asc');
fileNum = length(fileName);

for i =1:fileNum
    fileNamei = fileName(i).name;
    trialInfoi = strrep(fileNamei,'.asc','.mat');
    
    % extract the time on file names as the error marker
    num = regexp(fileNamei, '(\d+)','tokens');
    NUM = str2double(num{end});
    
    saveName = ['converted_' trialInfoi];
  
    % extrace data
    fid=fopen(fileNamei);
    fseek(fid,0,'eof');
    numline=ftell(fid);
    fclose(fid);
    rawdata=importdata(fileNamei,' ',numline);
    
    
    fo = regexp(rawdata(contains(rawdata,'fixation onset','IgnoreCase',true)),'(\d+)','tokens');
    fixation_onset = nan(length(fo),2);
    for j = 1 : length(fo)
        if ~isempty(fo{j})
            fixation_onset(j,:) = [str2double(fo{j}{1}) str2double(fo{j}{2})]; % fixation on
        end
    end
    fo_del =diff(fixation_onset(:,2));
    fixation_onset(fo_del==0,:) = [];
    
    pa = regexp(rawdata(contains(rawdata,'pupil adapting','IgnoreCase',true)),'(\d+)','tokens');
    pupil_adapting = nan(length(pa),2);
    for j = 1 : length(pa)
        if ~isempty(pa{j})
            pupil_adapting(j,:) = [str2double(pa{j}{1}) str2double(pa{j}{2})]; % pupil fitting
        end
    end
    pa_del = diff(pupil_adapting(:,2));
    pupil_adapting(pa_del==0,:) = [];
    
    so = regexp(rawdata(contains(rawdata,'stimulus onset','IgnoreCase',true)),'(\d+)','tokens');
    stimulus_onset = nan(length(so),2);
    for j = 1 : length(so)
        if ~isempty(so{j})
            stimulus_onset(j,:) = [str2double(so{j}{1}) str2double(so{j}{2})]; % stimulus on
        end
    end
    so_del = diff(stimulus_onset(:,2));
    stimulus_onset(so_del==0,:) = [];
    
    st = regexp(rawdata(contains(rawdata,'stimulus term','IgnoreCase',true)),'(\d+)','tokens');
    stimulus_term = nan(length(st),2);
    for j = 1 : length(st)
        if ~isempty(st{j})
            stimulus_term(j,:) = [str2double(st{j}{1}) str2double(st{j}{2})]; % stimulus term
        end
    end
    st_del = diff(stimulus_term(:,2));
    stimulus_term(st_del==0,:) = [];
    
    sc = regexp(rawdata(contains(rawdata,'start choice','IgnoreCase',true)),'(\d+)','tokens');
    start_choice = nan(length(sc),2);
    for j = 1 : length(sc)
        if ~isempty(sc{j})
            start_choice(j,:) = [str2double(sc{j}{1}) str2double(sc{j}{2})]; % start choice
        end
    end
    sc_del = diff(start_choice(:,2));
    start_choice(sc_del==0,:) = [];
    
    dm = regexp(rawdata(contains(rawdata,'decision made','IgnoreCase',true)),'(\d+)','tokens');
    decision_made = nan(length(dm),2);
    for j = 1 : length(dm)
        if ~isempty(dm{j})
            decision_made(j,:) = [str2double(dm{j}{1}) str2double(dm{j}{2})]; % decision making
        end
    end
    dm_del = diff(decision_made(:,2));
    decision_made(dm_del==0,:) = [];
    
    te = regexp(rawdata(contains(rawdata,'trial end','IgnoreCase',true)),'(\d+)','tokens');
    trial_end = nan(length(te),2);
    for j = 1 : length(te)
        if ~isempty(te{j})
            trial_end(j,:) = [str2double(te{j}{1}) str2double(te{j}{2})]; % trial end
        end
    end
    te_del = diff(trial_end(:,2));
    trial_end(te_del==0,:) = [];
    
    nc = regexp(rawdata(contains(rawdata,'no choice making','IgnoreCase',true)),'(\d+)','tokens');
    no_choice = nan(length(nc),2);
    for j = 1 : length(nc)
        if ~isempty(nc{j})
            no_choice(j,:) = [str2double(nc{j}{1}) str2double(nc{j}{2})]; % no choice making, especiall for no choice version as trial end
        end
    end
    nc_del = diff(no_choice(:,2));
    no_choice(nc_del==0,:) = [];
    
    ind=strfind(rawdata,'...');
    ind=~cellfun(@isempty,ind);
    data=rawdata(ind);
    clear rawdata;
    data=strrep(data,'...','');
    data=cellfun(@str2num,data, 'UniformOutput',false);
    positionData=cell2mat(data);
    
    eyeData = positionData(:,1:4);
    clear data;
    clear positionData;
    clear ind;
    
%     % check if it is 2000 Hz, convert to 1000 Hz
%     time=eyeData(:,1);
%     d_data = [];
%     if ~isempty(find(unique(diff(time))==0, 1))
%         for k=1:2
%             if time(k)==time(k+1)
%                 d_data = eyeData(k:2:end,:);
%             end
%         end
%         eyeData = d_data;
%         time = eyeData(:,1);
%         status = cat(1,status,[0 NUM]);
%     end
    
    % replace all the mission data as zero
    time = eyeData(:,1);
    dt=mode(diff(time));
    ind=(time-time(1))/dt+1;
    ndata=zeros(ind(end),size(eyeData,2))*(-1);
    ndata(:,1)=time(1):dt:time(end);
    ndata(ind,2:end)=eyeData(:,2:end);
    
    clear eyeData
    
    figure(50);clf
    plot(ndata(:,1),ndata(:,2),'r');
    hold on
    plot(ndata(:,1),ndata(:,3),'g');
    plot(ndata(:,1),ndata(:,4),'b');
    
    if data_check == 1
        [data_purified,~,errorflagj,purifyNUMj] = BlinkNoisePurify_NaN(ndata,4,dt); % remove the blink data, replaced by NaN

        clear ndata
        error_set(error_set == 7) = [];
        if ~isempty(intersect(errorflagj,error_set))  || purifyNUMj > 3
    %         if intersect(errorflagj,error_set) == 7
    %             fprintf(2,[fileName(i).name ' base line calculation skipped by eye blink. \n'])
    %         else
                fprintf(2,[fileName(i).name ' block skipped by error ' num2str(intersect(errorflagj,error_set)') ' ,  by checked ' num2str(purifyNUMj) ' times.\n'])
    %         end
            delete(fullfile(filespath,saveName));
            continue
        end
    else
        [data_purified,~,~,~] = BlinkNoisePurify_NaN(ndata,4,dt); % remove the blink data, replaced by NaN
    end
    
    pupil_eyedata = cell(length(decision_made),1);
    saccade_eyedata = cell(length(decision_made),1);
    baseline = nan(length(decision_made),1);
    base_std = nan(length(decision_made),1);
    base_se = nan(length(decision_made),1);
    
    for j = decision_made(:,2)'
        %% collect eye data for pupil size
        trialST = find(data_purified(:,1) >= stimulus_onset(stimulus_onset(:,2) == j),1);
        if isempty(trial_end)
            trialEND = find(data_purified(:,1) <= no_choice(stimulus_onset(:,2) == j),1,'last'); % if no choice version
        else
            trialEND = find(data_purified(:,1) <= trial_end(trial_end(:,2) == j),1,'last');
        end
        pupil_eyedata{j} = data_purified(trialST:trialEND,:);
        
%         take_point = find(data_purified(:,1) >= decision_made(j,1),1);
%         point_before = 1500; % ms
%         point_after = 500; % ms
%         pupil_eyedata{j} = data_purified(take_point- point_before/dt :take_point+ point_after/dt,:);
        
        %% collect eye data for saccade
        saccadeST = find(data_purified(:,1) >= stimulus_onset(stimulus_onset(:,2) == j),1);
        if isempty(decision_made)
            saccadeEND = find(data_purified(:,1) <= no_choice(no_choice(:,2) == j),1,'last'); % if no choice version
        elseif isnan(decision_made(j,1))
            continue
%             saccadeEND = find(data_purified(:,1) <= no_choice(no_choice(:,2) == j),1,'last'); % if no choice made
        else
            saccadeEND = find(data_purified(:,1) <= trial_end(trial_end(:,2) == j),1,'last');
        end
        saccade_eyedata{j} = data_purified(saccadeST:saccadeEND,:);
        
        %% calculation for pupil size base line
        
        pre_baseT = pupil_baseline_duration / dt; % how many ms for base line calculation
        
        if stimulus_onset(stimulus_onset(:,2) == j) - pre_baseT <= pupil_adapting(pupil_adapting(:,2) == j)
            error('Not able to yield credible pupil base line, please check the pupil adapting duration and base line sampling instant')
        end
        
        baseEND = find(data_purified(:,1) <= stimulus_onset(stimulus_onset(:,2) == j),1,'last');
        
        baseset = data_purified(baseEND-pre_baseT:baseEND,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if sum(isnan(baseset)) > 0
            fprintf(2,['Base line calculation skipped in trial' num2str(j) fileName(i).name ' for eye blinking. \n'])
            continue
        else
            baseline(j) = nanmean(baseset(:,4));
            base_std(j) = nanstd(baseset(:,4));
            base_se(j) = base_std(j) / sqrt(pre_baseT / dt);
        end
%         fprintf(1,['S.E / Mean for Trial ' num2str(j) ' in ' fileName(i).name ' is ' num2str(base_se(j)) '/' num2str(baseline(j)) '\n'])
    end
        
    
    save(saveName,'pupil_eyedata','saccade_eyedata','baseline','base_std','dt') %,'fixation_on','pupil_fit','stimulus_on','stimulus_term','start_choice','decision_making','trial_end','no_choice');
    
    status = cat(1,status,[1 NUM]);

    clear time
    clear ndata
    clear eyeData baseline base_std dt fixation_onset pupil_adapting stimulus_onset stimulus_term start_choice decision_made trial_end no_choice
end

% %% reverse for the error report
% cd(filespath);
% error_report = status(status(:,1)~=1,:);
% if ~isempty(error_report)
%     report_savename = [datestr(now,'yymmddHHMM') '_getEyeDataFV.mat'];
%     save(report_savename,'status','error_report')
% end

