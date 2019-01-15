function getEyeDatas(dataPath,error_set,no_choice_version)

if ~exist('dataPath','var')
    error('Please import the path where is your data.');
else
    cd(dataPath);
    % find all the asc files in this folder
    fileName = dir('*.asc');
    if isempty(fileName)
        error('There is no asc file in this path.');
    end
end

status = [];
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
    fixation_on = nan(length(fo),2);
    for j = 1 : length(fo)
        if ~isempty(fo{j})
            fixation_on(j,:) = [str2double(fo{j}{1}) str2double(fo{j}{2})]; % fixation on
            if j > 1
                if fixation_on(j,2) == fixation_on(j-1,2)
                    fixation_on(j-1,:) = NaN;
                end
            end
        end
    end
    fixation_on(isnan(fixation_on(:,1)),:) = [];
%     
    pf = regexp(rawdata(contains(rawdata,'pupil adapting','IgnoreCase',true)),'(\d+)','tokens');
% pf = regexp(rawdata(contains(rawdata,'pupil fitting','IgnoreCase',true)),'(\d+)','tokens');
    pupil_fit = nan(length(pf),2);
    for j = 1 : length(pf)
        if ~isempty(pf{j})
            pupil_fit(j,:) = [str2double(pf{j}{1}) str2double(pf{j}{2})]; % pupil fitting
            if j > 1
                if pupil_fit(j,2) == pupil_fit(j-1,2)
                    pupil_fit(j-1,:) = NaN;
                end
            end
        end
    end
    pupil_fit(isnan(pupil_fit(:,1)),:) = [];
    
    so = regexp(rawdata(contains(rawdata,'stimulus onset','IgnoreCase',true)),'(\d+)','tokens');
    stimulus_on = nan(length(so),2);
    for j = 1 : length(so)
        if ~isempty(so{j})
            stimulus_on(j,:) = [str2double(so{j}{1}) str2double(so{j}{2})]; % stimulus on
            if j > 1
                if stimulus_on(j,2) == stimulus_on(j-1,2)
                    stimulus_on(j-1,:) = NaN; 
                end
            end
        end
    end
    stimulus_on(isnan(stimulus_on(:,1)),:) = [];
        
    st = regexp(rawdata(contains(rawdata,'stimulus term','IgnoreCase',true)),'(\d+)','tokens');
    stimulus_term = nan(length(st),2);
    for j = 1 : length(st)
        if ~isempty(st{j})
            stimulus_term(j,:) = [str2double(st{j}{1}) str2double(st{j}{2})]; % stimulus term
            if j > 1
                if stimulus_term(j,2) == stimulus_term(j-1,2)
                    stimulus_term(j-1,:) = NaN;
                end
            end
        end
    end
    stimulus_term(isnan(stimulus_term(:,1)),:) = [];
    
    sc = regexp(rawdata(contains(rawdata,'start choice','IgnoreCase',true)),'(\d+)','tokens');
    start_choice = nan(length(sc),2);
    for j = 1 : length(sc)
        if ~isempty(sc{j})
            start_choice(j,:) = [str2double(sc{j}{1}) str2double(sc{j}{2})]; % start choice
            if j > 1
                if start_choice(j,2) == start_choice(j-1,2)
                    start_choice(j-1,:) = NaN;
                end
            end
        end
    end
    start_choice(isnan(start_choice(:,1)),:) = [];
%     
    dm = regexp(rawdata(contains(rawdata,'decision made','IgnoreCase',true)),'(\d+)','tokens');
%  dm = regexp(rawdata(contains(rawdata,'decision making','IgnoreCase',true)),'(\d+)','tokens');
    decision_making = nan(length(dm),2);
    for j = 1 : length(dm)
        if ~isempty(dm{j})
            decision_making(j,:) = [str2double(dm{j}{1}) str2double(dm{j}{2})]; % decision making
        end
    end
    decision_making(isnan(decision_making(:,1)),:) = [];
    dm_index = find(diff(decision_making(:,2)) == 1);
    decision_making = decision_making([1;dm_index+1],:);
    
    te = regexp(rawdata(contains(rawdata,'trial end','IgnoreCase',true)),'(\d+)','tokens');
    trial_end = nan(length(te),2);
    for j = 1 : length(te)
        if ~isempty(te{j})
            trial_end(j,:) = [str2double(te{j}{1}) str2double(te{j}{2})]; % trial end
            if j > 1
                if trial_end(j,2) == trial_end(j-1,2)
                    trial_end(j-1,:) = NaN;
                end
            end
        end
    end
    trial_end(isnan(trial_end(:,1)),:) = [];
    
    nc = regexp(rawdata(contains(rawdata,'no choice making','IgnoreCase',true)),'(\d+)','tokens');
    no_choice = nan(length(nc),2);
    for j = 1 : length(nc)
        if ~isempty(nc{j})
            no_choice(j,:) = [str2double(nc{j}{1}) str2double(nc{j}{2})]; % no choice making, especiall for no choice version as trial end
%             if j > 1
%                 if no_choice(j,2) == no_choice(j-1,2)
%                     no_choice(j-1,:) = NaN;
%                 end
%             end
        end
    end
     no_choice(isnan(no_choice(:,1)),:) = [];
    
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
    
    % check if 2000 Hz, convert as 1000 Hz
    time=eyeData(:,1);
    d_data = [];
    if ~isempty(find(unique(diff(time))==0, 1))
        for k=1:2
            if time(k)==time(k+1)
                d_data = eyeData(k:2:end,:);
            end
        end
        eyeData = d_data;
        time = eyeData(:,1);
        status = cat(1,status,[0 NUM]);
    end
    
    % replace all the mission data as 0
    dt=mode(diff(time));
    ind=(time-time(1))/dt+1;
    ndata=zeros(ind(end),size(eyeData,2))*(-1);
    ndata(:,1)=time(1):dt:time(end);
    ndata(ind,2:end)=eyeData(:,2:end);
    
    clear eyeData
    
    pupil_eyedata = cell(length(pupil_fit),1);
    saccade_eyedata = cell(length(pupil_fit),1);
    baseline = nan(length(pupil_fit),1);
    base_std = nan(length(pupil_fit),1);
    base_se = nan(length(pupil_fit),1);
    skipped_trial_num = 0;
    
    triali = decision_making(:,2)';
    for j = triali(~isnan(triali))
        %% collect eye data for pupil size analyses
        % set the begining point of your analysis period in the trials for pupil size.
        begining_point_pupil = stimulus_on;
        
        pupilST = find(ndata(:,1) >= begining_point_pupil(begining_point_pupil(:,2) == j,1),1);
        pupilST = pupilST + 0; % 改动前后多少ms换这个数字
        
        % set the ending point of your analysis period in the trials for pupil size.
        ending_point_pupil = trial_end;
        
        if isempty(ending_point_pupil)
%             trialEND = find(ndata(:,1) <= no_choice(no_choice(:,2) == j,1),1,'last'); % for no choice version
            error('there is something wrong with your ending point for your analysis.')
        else
            pupilEND = find(ndata(:,1) <= ending_point_pupil(ending_point_pupil(:,2) == (j),1),1,'last'); % end point
%              trialEND = find(ndata(:,1) <= decision_making(decision_making(:,2) == j,1),1,'last'); % end point
            pupilEND = pupilEND + 0; % 改动前后多少ms换这个数字
        end
        pupil_eyedata{j} = ndata(pupilST:pupilEND,:);
        
        %% collect eye data for saccade / micro-saccade analyses
        % set the begining point of your analysis period in the trials for saccade.
        begining_point_saccade = stimulus_on;
        
        saccadeST = find(ndata(:,1) >= begining_point_saccade(begining_point_saccade(:,2) == j,1),1); % start point
        saccadeST= saccadeST + 0; % 改动前后多少ms换这个数字
        
        % set the ending point of your analysis period in the trials for saccade.
        ending_point_saccade = trial_end;
        
        if isempty(ending_point_saccade)
%             saccadeEND = find(ndata(:,1) <= no_choice(no_choice(:,2) == j,1),1,'last'); % if no choice version
            error('there is something wrong with your ending point for your analysis.')
        else
             saccadeEND = find(ndata(:,1) <= ending_point_saccade(ending_point_saccade(:,2) == (j),1),1,'last'); % end point
             saccadeEND = saccadeEND + 0; % 改动前后多少ms换这个数字
        end
        saccade_eyedata{j} = ndata(saccadeST:saccadeEND,:);
        
        %% calculation for pupil size base line
        baseEND = find(ndata(:,1) >= stimulus_on(stimulus_on(:,2) == j,1),1);
        
        pre_baseT = 200; % how many ms for base line calculation
        
        baseset = ndata(baseEND-pre_baseT:baseEND,:);
        [baseset_purified,~,errorflagj,purifyNUMj] = BlinkNoisePurify_NaN(baseset,dt,error_set,4); % remove the blink data, replaced by NaN
        if ~isempty(intersect(errorflagj,error_set))  || purifyNUMj > 3
            if intersect(errorflagj,error_set) == 7
                fprintf(2,['Trial ' num2str(j) ' in ' fileName(i).name ' base line calculation skipped by eye blink. \n'])
            else
                fprintf(2,['Trial ' num2str(j) ' in ' fileName(i).name ' base line calculation skipped by error ' num2str(intersect(errorflagj,error_set)) '\n'])
            end
            continue
        else
            baseline(j) = nanmean(baseset_purified(:,4));
            base_std(j) = nanstd(baseset_purified(:,4));
            base_se(j) = base_std(j) / sqrt(pre_baseT / dt);
%             fprintf(1,['S.E / Mean for Trial ' num2str(j) ' in ' fileName(i).name ' is ' num2str(base_se(j)) '/' num2str(baseline(j)) '\n'])
        end
    end
        
    
    save(saveName,'pupil_eyedata','saccade_eyedata','baseline','base_std','dt','skipped_trial_num','decision_making','stimulus_on','stimulus_term','start_choice','decision_making','trial_end','no_choice');
    
    status = cat(1,status,[1 NUM]);

    clear time
    clear ndata
    clear eyeData baseline base_std dt fixation_on pupil_fit stimulus_on stimulus_term start_choice decision_making trial_end no_choice
end

% %% reverse for the error report
% cd(filespath);
