% function [treated_data,blink_data,resultFlag]=blink_treatment_pupil1(eyedata,pupilcol,NUM)
function [blink_time_data,blink_time,resultFlag]=blink_treatment_pupil1(eyedata,pupilcol,dt)
% This function can eliminate eye blinking data or squinting data based on pupil size.
% the first input is for eye data in column.
% Column 1 is for time, column 2 is for x pixel, column 3 is for y pixel, column 4 is for pupil size, in default.
%
% BY BYC 06/SEP/2018
%
% meaning of result flags:
error_0 = 0; % 0. no eye blink;
error_1 = 1; % 1. small pupil size, the experiments may not operated in a strict darkness room;
error_2 = 2; % 2. eye closeing in end;
error_3 = 3; % 3. eye closeing in recording begining;
error_4 = 4; % 4. sudden increase / decrease;
error_5 = 5; % 5. at least one peak within/just before/soon after the eye blink;
error_6 = 6; % 6. not firmly closeing in eye blink or just squinting;
error_7 = 7; % 7. A long time with eye closing (1/5);
error_8 = 8; % 8. have a minimum of pupil size & blink detected;
error_9 = 9; % 9. Still have unknow rifts in x and y pixel, that should be caused by not rigorous operation or other unknown conditions.
error_10=10; % 10.Still lots of points outside the screen. It may be caused by no correctly calibration & validation, or may be the subject squinting for too long (>1s)
error_11=11; % 11.Still have unknown noisy point (sudden decrease to zero for 1 or 2 samples, discontinuously), this maybe caused by Eyelink.
error_12=12; % 12.Input data too short

pupildata = eyedata(:,pupilcol);
x_col = 2;
y_col = 3;

global SCREEN; % global screen in weight, height, distance, smooth

if ~exist('SCREEN.width','var')
    SCREEN.width = 1280;
end

if ~exist('SCREEN.height','var')
    SCREEN.height = 1024;
end

resultFlag = [];
onset_set=[];
term_set=[];
blink_time = [];
mean_pupil = nanmean(pupildata);
std_pupil = nanstd_median(pupildata);
time_adding_onset = 150; % eliminate 150 ms more before blink onset
time_adding_term = 150; % eliminate 150 ms more after blink termination

if mean_pupil <= 1000  % reference needed for the pupil size under light versus a strict darkness room
    resultFlag = cat(1,resultFlag,error_1); % this value only refer to ~60 cm as head to screen distance.
end

% remove noise which sudden to 0 in 1 ms
for i = 3:length(pupildata)-2
    outlier_check = isoutlier(pupildata(i-2:i+2));
    if outlier_check(3) == 1
        pupildata(i) = (pupildata(i-2) + 2*pupildata(i-1) + 2*pupildata(i+1) + pupildata(i+2)) / 6;
        resultFlag = cat(1,resultFlag,error_11);
    end
end

%% calculate the smoothed clopes
onset_set = zeros(length(pupildata),1);

for i = 1 : length(pupildata)
    if i == 1 || i == length(pupildata)
        onset_set(i,:) = 0;
    elseif i == 2 || i == length(pupildata) - 1
        onset_set(i) = (pupildata(i-1) - pupildata(i+1)) / 2;
    elseif i == 3 || i == length(pupildata) - 2
        onset_set(i) = (pupildata(i-2) + pupildata(i-1) - pupildata(i+1) - pupildata(i+2)) / 6;
    else
        onset_set(i) = (pupildata(i-3) + pupildata(i-2) + pupildata(i-1) - pupildata(i+1) - pupildata(i+2) - pupildata(i+3)) / 12;
%     else
%         onset_set(i,:) = (eyedata(i+4,[col_x,col_y]) + eyedata(i+3,[col_x,col_y]) + eyedata(i+2,[col_x,col_y]) + eyedata(i+1,[col_x,col_y]) - eyedata(i-1,[col_x,col_y]) - eyedata(i-2,[col_x,col_y]) - eyedata(i-3,[col_x,col_y]) - eyedata(i-4,[col_x,col_y])) / 20;
    end
end
term_set = rot90(-onset_set,2);

% check for length of data
if length(onset_set) <= 50
    blink_time_data = [];
    blink_time = [];
    resultFlag = cat(1,resultFlag,error_12);
    return
end

% peak_thres = nanstd_median(onset_set);
[~,~,peak_thres,~] = isoutlier(onset_set,'median'); % this function employ 3 * SD as threshold

% figure(2000);clf
% plot(term_set,'c')
% hold on
% plot([0 length(term_set)],[peak_thres*3 peak_thres*3],'--k')
% plot(onset_set,'k')
% plot(eyedata(:,2),'r')
% plot(eyedata(:,3),'g')
% plot(eyedata(:,4),'b')


% try to find the blink point every 50ms,
[~,index_onset] = findpeaks(onset_set,'minPeakHeight',peak_thres*3,'minPeakDistance',50/dt);
[~,index_term_r] = findpeaks(term_set,'minPeakHeight',peak_thres*3,'minPeakDistance',50/dt);
index_term = rot90((length(term_set)-index_term_r),2);
% pks_term = rot90(pks_term_r,2);

%% chack for error and validation
if ~isempty(index_onset) && ~isempty(index_term)
% check for 3 and mark in resultFlag: eye closeing in recording start;
    if sum(pupildata(1:min(index_term)) == 0) ~= 0 || min(index_term) <= time_adding_term
        index_onset = [1; index_onset];
        resultFlag = cat(1,resultFlag,error_3);
    elseif sum(diff(pupildata(1:min(index_term))) < 0) == 0
        index_onset = [1; index_onset];
        resultFlag = cat(1,resultFlag,error_3);
    end

    % check for 2 and mark in resultFlag: eye closeing in end;
    if sum(pupildata(max(index_term):end) == 0) ~= 0
        index_term = [index_term; length(pupildata)];
        resultFlag = cat(1,resultFlag,error_2);
    elseif sum(diff(pupildata(max(index_term):end)) > 0) == 0
        index_term = [index_term; length(pupildata)];
        resultFlag = cat(1,resultFlag,error_2);
    end

    % pair eye blink points
    paird_term = zeros(size(index_onset));
    for i = 1:length(index_onset)
        paird_index_term = find(index_term > index_onset(i),1);
        if ~isempty(paird_index_term)
            paird_term(i) = index_term(paird_index_term);
        end
    end
    paird_term = paird_term(paird_term~=0);

    % check for 4 and mark in resultFlag: sudden increase
    if sum(diff(paird_term)==0) > 0
        paird_term_r = paird_term(end:-1:1);
        paird_term_delIndex = length(paird_term)-find(diff(paird_term_r)==0); % u can record increase point if needed
        paird_term(paird_term_delIndex) = [];
        resultFlag = cat(1,resultFlag,error_4);
    end
    
    % pair eye blink points
    paird_onset = zeros(size(paird_term));
    for i=1:length(paird_term)
        paird_index_onset = find(index_onset < paird_term(i),1,'last');
        if ~isempty(paird_index_onset)
            paird_onset(i) = index_onset(paird_index_onset);
        end
    end
    paird_onset = paird_onset(paird_onset~=0);

    % check for 4 and mark in resultFlag: sudden decrease
    if sum(diff(paird_onset)==0) > 0
        paird_onset(diff(paird_onset)==0) = [];
        resultFlag = cat(1,resultFlag,error_4);
    end
    
    paird_onset = paird_onset - time_adding_onset;
    paird_term = paird_term + time_adding_term;
  
    for j = 1:length(paird_onset)
        if paird_onset(j) < 1
            paird_onset(j) = 1;
        end
    end

    for j = 1:length(paird_term)
        if paird_term(j) > length(pupildata)
            paird_term(j) = length(pupildata);
        end
    end

    paird_onset_fin =[];
    paird_term_fin =[];
    for i = 1:length(paird_term)
        if sum(pupildata(paird_onset(i):paird_term(i))==0)>0 
            paird_onset_fin = cat(1,paird_onset_fin,paird_onset(i));
            paird_term_fin = cat(1,paird_term_fin,paird_term(i));
        elseif paird_term(i) - paird_onset(i) <= 1200 % more reference needed for short squinting
            paird_onset_fin = cat(1,paird_onset_fin,paird_onset(i));
            paird_term_fin = cat(1,paird_term_fin,paird_term(i));
            resultFlag = cat(1,resultFlag,error_6); % not firmly closeing in eye blink or squint;
        else 
            if std(paird_onset(i)+100:paird_term(i)-100) < std_pupil/2
                paird_onset_fin = cat(1,paird_onset_fin,paird_onset(i));
                paird_term_fin = cat(1,paird_term_fin,paird_term(i));
                resultFlag = cat(1,resultFlag,error_8); % not sure, squinting for more than 1s or a pair of sudden increase/decrease, need additional check
            end
            continue
        end
    end
    
    if isempty(paird_onset_fin) || isempty(paird_term_fin)
        resultFlag = cat(1,resultFlag,error_0); % 0. no eye blink;
        blink_index_pre = [];
        blink_pair = [];
    else
        
        % record index for blink and treated data
        blink_pair = zeros(length(paird_term_fin),2);
        blink_index_pre = [];
        for i = 1:length(paird_term_fin)
            if i == 1
                if paird_onset_fin(1) == 1
                    blink_index_pre = cat(1,blink_index_pre,(1:paird_term_fin(1))');
                    blink_pair(i,:) = [paird_onset_fin(i) paird_term_fin(i)];
                    continue
                end
            end
            blink_index_pre = cat(1,blink_index_pre,(paird_onset_fin(i):paird_term_fin(i))');
            blink_pair(i,:) = [paird_onset_fin(i) paird_term_fin(i)];  
        end
        
        blink_time = [eyedata(blink_pair(:,1),1) eyedata(blink_pair(:,2),1)];
        
        % check for: at least one peak within/just before/soon after the eye blink;
        for i = 1:size(blink_pair,1)
            if sum(pupildata(blink_pair(i,1):blink_pair(i,2)) > max(pupildata(blink_pair(i,:))) + std_pupil) > 0 % change std_pupil here to chang the threshold
            resultFlag = cat(1,resultFlag,error_5);
            end
        end
    end

elseif isempty(index_term) && ~isempty(index_onset)
    if sum(pupildata(max(index_onset):end)==0) > 0
        resultFlag = cat(1,resultFlag,error_2); % 2. eye closeing in end;
        index_onset = max(index_onset - time_adding_onset,1);
        blink_index_pre =  (max(index_onset) : length(pupildata))';
        blink_pair = [max(index_onset) length(pupildata)];
    elseif sum(diff(pupildata(max(index_onset) : end)) > 0) == 0
        blink_index_pre = (max(index_onset) - time_adding_onset : length(pupildata))';
        index_onset = max(index_onset - time_adding_onset,1);
        blink_index_pre =  (max(index_onset) : length(pupildata))';
        blink_pair = [max(index_onset) length(pupildata)];
        resultFlag = cat(1,resultFlag,error_2);
    else
        resultFlag = cat(1,resultFlag,error_0); % 0. no eye blink;
        blink_index_pre = [];
        blink_pair = [];
    end
elseif isempty(index_onset) && ~isempty(index_term)
    if sum(pupildata(1:min(index_term))==0) > 0
        resultFlag = cat(1,resultFlag,error_3); % 3. eye closeing in recording start;
        index_term = min(index_term + time_adding_term,length(pupildata));
        blink_index_pre = (1:min(index_term))';
        blink_pair = [1 min(index_term)];
    elseif sum(diff(pupildata(1 : min(index_term))) < 0) == 0
        index_term = min(index_term + time_adding_term,length(pupildata));
        resultFlag = cat(1,resultFlag,error_3);
        blink_index_pre = (1:min(index_term))';
        blink_pair = [1 min(index_term)];
    else
        resultFlag = cat(1,resultFlag,error_0); % 0. no eye blink;
        blink_index_pre = [];
        blink_pair = [];
    end
else
    resultFlag = cat(1,resultFlag,error_0); % 0. no eye blink;
    blink_index_pre = [];
    blink_pair = [];
end 



pupil_raw = pupildata;
pupil_raw(blink_index_pre) = [];
additional_del = find(pupil_raw == 0);
if ~isempty(additional_del)
    blink_index = cat(1,blink_index_pre,additional_del);
    blink_index = unique(blink_index);
    resultFlag = cat(1,resultFlag,error_11);
else
    blink_index = blink_index_pre;
end

x_pixel = eyedata(:,x_col);
y_pixel = eyedata(:,y_col);
x_pixel(blink_index) = [];
y_pixel(blink_index) = [];

if sum(diff(x_pixel) > 800) + sum(diff(y_pixel) > 800) > 0
    resultFlag = cat(1,resultFlag,error_9);
end

if sum(or(or(or(x_pixel < -50, x_pixel > SCREEN.width+50), y_pixel < -50), y_pixel > SCREEN.height+50)) > 0
    resultFlag = cat(1,resultFlag,error_10);
end

treated_index = 1:length(pupildata);
treated_index(blink_index) = [];
treated_data = eyedata(treated_index,1);
resultFlag = unique(resultFlag);
blink_time_data = eyedata(blink_index,1);
blink_pair_index = [];

if ~isempty(blink_pair)
    for i = 1 : size(blink_pair,1)
        blink_pair_index = cat(1,blink_pair_index,[eyedata(blink_pair(i,1),1) eyedata(blink_pair(i,2),1)]);
    end
end
        

if length(treated_index) <= length(pupildata) / 5 * 4
    resultFlag = cat(1,resultFlag,error_7);
end


% if sum(~ismember(resultFlag(1),0))
%     %% plot for test
%     figure(2500); clf; set(gcf,'color','white');
%     plot(pupildata);
%     hold on;
%     plot(onset_set,'r-');
%     plot(term_set,'g-');
%     plot(eyedata(:,2),'r');
%     plot(eyedata(:,3),'g')
%     plot(eyedata(:,4),'b');
%     if ~isempty(blink_pair)
%         plot(blink_pair(:,1),pupildata(blink_pair(:,1)),'ro');
%         plot(blink_pair(:,2),pupildata(blink_pair(:,2)),'go');
%         ylim([-200,inf]);
%     %     title(NUM);
%         plot(blink_index,ones(size(blink_index))*(-100),'r.');
%     end
%     plot(treated_index,ones(size(treated_index))*(-120),'k.');
%     % plot(blink_index,ones(size(blink_index)),'bs');
%     % plot(treated_index,ones(size(treated_index)),'ks');
%     hold off
%     % test end
% end
end

function output = nanstd_median(input)
output = sqrt( nansum( power( abs(input - median(input,'omitnan')),2 )) / (sum(~isnan(input))-1) );
end
    