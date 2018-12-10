function main
% First edition BY BYC at 2018/10/16
%
% meaning of error flags:
% 0. no eye blink;
% 1. small pupil size, the experiments may not operated in a strict 
%     darkness room;
% 2. eye closeing in end;
% 3. eye closeing in recording begining;
% 4. sudden increase / decrease;
% 5. at least one peak within/just before/soon after the eye blink;
% 6. not firmly closeing in eye blink or just squinting;
% 7. A long time with eye closing (1/5);
% 8. have a minimum of pupil size & blink detected;
% 9. still have unknow rifts in x and y pixel, that should be caused
%     by not regorous operation or other unknown conditions.
% 10.Still lots of points outside the screen. It may be caused by 
%     no correctly calibration & validation, or may be the subject squinting for too long (>1s)
% 11.Still have unknown noisy point (sudden decrease to zero for ~1
%     to 20ms, discontinuously), this maybe caused by Eyelink.
%
% BY BYC SEP/2018
%
% Last updated at 2018/11/16 by BYC

close all
% clear all
warning off

data_path = ''; % the path for candidate data
exeFilePath = ''; % the path for exe file 'edf2asc.exe'
fileSaveName = 'sample_';

error_set = [7 8 9 10]; % for more detail, please read BlinkNoisePurify_NaN.m
maxtrialNUM = 500; % the maximum trial numbers in all blocks
smooth_Hz = 60; % the aim value for data smooth / low pass
fixation_degree = 3; % fixation degree
saccade_thres = 30; % degree/s, 30 as example
saccade_thres_cap = 100; % degree/s, 100 as example
% microsaccade_thres = 5; % times of S.D as threshold.
pupil_baseline_duration = 200; % ms; the duration of time for calculating base line for pupil size before the period of stimulus onset.

skip_files_convert = 1; % 0: need do converting; 1: skip files conversion. if EDF files have already been converted, set 1 to saving your time
skip_eyedata_collect = 1; % 0: need do converting; 1: skip eye data collection. if eye data have already been converted, set 1 to saving your time

%% calculation for fixation windows
global SCREEN
SCREEN.screenWidth = 1280; % pixel
SCREEN.screenHeight = 1024; % pixel
SCREEN.view_distance = 60; % cm
SCREEN.screenWidth_real = 37.5; % cm
SCREEN.screenHeight_real = 30; % cm
pixel_x = tand(fixation_degree) * SCREEN.view_distance / SCREEN.screenWidth_real * SCREEN.screenWidth;
pixel_y = tand(fixation_degree) * SCREEN.view_distance / SCREEN.screenHeight_real * SCREEN.screenHeight;

if pixel_x ~= pixel_y
    error('可能未用正屏显示器采样，或参数输入错误')
end
% view_range = [screenWidth/2-pixel_x screenWidth/2+pixel_x screenHeight/2-pixel_y screenHeight/2+pixel_y];

%% set config for subplots
set(figure(1),'pos',[27 63 1849 892],'Name','Pupil size');clf;
subplot1 = tight_subplot(4,4,[0.03 0.03]);
% set(subplot1,'ytick',[]);
suptitle('Pupil size');

set(figure(2),'pos',[27 63 1849 892],'Name','Saccade track');clf;
subplot2 = tight_subplot(4,4,[0.03 0.03]);
% set(subplot2,'ytick',[]);
suptitle('Saccade track');

set(figure(3),'pos',[27 63 1849 892],'Name','Saccade rate');clf;
subplot3 = tight_subplot(4,4,[0.03 0.03]);
% set(subplot2,'ytick',[]);
suptitle('Saccade rate');

%% convert EDF files to readable eye data mat files, it may takes a long time
if skip_files_convert == 0
    edf2asc_checkasc(data_path,exeFilePath);
end
if skip_eyedata_collect == 0
    getEyeDatas(data_path,error_set,pupil_baseline_duration);
end

% read the converted mat files
cd(data_path);
files = dir('converted_*');
fileNUM = length(files);
blinkNUM = NaN(maxtrialNUM,fileNUM);


for i = 1 : fileNUM % for each block
    load(files(i).name)
    trialInfo = strrep(files(i).name,'converted_','');
    trialInfo = load(trialInfo);

    pupil_plot_set = nan(length(pupil_eyedata),max(max(cell2mat(cellfun(@size,pupil_eyedata,'UniformOutput',false)))));%生成NaN矩阵(trial数,最多项的列数)
    plot_x_P = (1: max(cell2mat(cellfun(@size,pupil_eyedata,'UniformOutput',0)))) * dt;%时间轴
    
%     trial_choice = trialInfo.choice == 2; % 1:close; 2: far
%     trial_time = trialInfo.time == 3/2; % 用分数表示
%     combined_c_t = and(trial_choice,trial_time)
%     trial_distance = trialInfo.distance == ?; % 用分数表示
%     trial_velocity = ( trialInfo.distance ./ trialInfo.time == 0.15 ); % 可用小数表示
    
    test = cell2mat(pupil_eyedata);
    if ~isempty(test)
        figure(50);clf;
        plot(test(:,4));
        hold on
        plot(test(:,2),'r');
        plot(test(:,3),'g');
    end
    
    %% data pre-processing for pupil data
    for j = 1 : length(pupil_eyedata) % for each trial
%     for j = find(trial_choice == 2) % for chosen trial 上面算出来的是逻辑值，所以永远都是取 == 1，只要改前面的参数名。在saccade哪里也要做对应更改
        trial_Peyedata = pupil_eyedata{j};
        smooth_level = 1000 / mean(diff(trial_Peyedata(:,1))) / smooth_Hz;
        
%         figure(5);clf;
%         plot(trial_Peyedata(:,4),'b');
%         hold on
%         plot(trial_Peyedata(:,2),'r');
%         plot(trial_Peyedata(:,3),'g');
                
        if ~isempty(trial_Peyedata)
            % check for data validation
            if sum(isnan(trial_Peyedata(:,2))) > size(trial_Peyedata,1) / 5  % more than 1/5 of data were eliminated
                disp(['Trial ' num2str(j) ' in ' files(i).name ' skipped.'])
                continue
            end
        else
            disp(['Something wrong here, no data for trial ' num2str(j) ' in ' file(i).name ])
            continue
        end
         
%         figure(51);clf;
%         plot(trial_Peyedata(:,4));
%         hold on
%         plot(trial_Peyedata(:,2),'r');
%         plot(trial_Peyedata(:,3),'g');
        
%         trial_Peyedata_a = trial_Peyedata;
%         trial_Peyedata(:,2) = smoothdata(trial_Peyedata(:,2),'gaussian',5,'includenan');
%         trial_Peyedata(:,3) = smoothdata(trial_Peyedata(:,3),'gaussian',5,'includenan');
        
%         figure(52);clf;
%         plot(trial_Peyedata_a(:,4));
%         hold on
%         plot(trial_Peyedata_a(:,2),'r');
%         plot(trial_Peyedata_a(:,3),'g');
%         
%         test = trial_Peyedata_a(:,2:3) == trial_Peyedata(:,2:3);
%         find(test == 1)
%         % count for blink times
%         [blinkPNUM(j,i),~] = size(blinktimesj); 
        
        % prepare for plot
        if ~isnan(baseline(j))
%             pupil_plot_set(j,1:length(trial_Peyedata(:,4))) = smooth(trial_Peyedata(:,4),smooth_level)'-baseline(j);
            pupil_plot_set(j,1:length(trial_Peyedata(:,4))) = ( trial_Peyedata(:,4)'-baseline(j) ) / baseline(j);
        end         
    end % trial end
    
    %% plot graphs with shaded error bar for pupil size
    
    if isempty(baseline)
        disp([files(i).name ' is skipped.'])
        continue
    end
    
    axes(subplot1(i));
    shadedErrorBar(plot_x_P,pupil_plot_set,{@nanmedian,@nanstd},'lineprops', '-b');
%     shadedErrorBar([],pupil_plot_set,{@nanmedian,@(plot_x_P) nanstd(plot_x_P)*1.96},'lineprops',{'b-'}); %/ sqrt(size(pupil_plot_set,1))*1.96},'lineprops',{'b-'});
    hold on
%     yticks('auto');
%     xticks('auto');
%     set(subplot1(i),'ytick',[-500 -250 0 250 500 1000],'yticklabel',{'-500','-250','0','250','500','1000'},'ylim',[-750 1000]);
    set(subplot1(i),'ytick',[-0.25 0 0.25 ],'yticklabel',{'-0.25','0','0.25'},'ylim',[-0.25 0.25]);
    set(subplot1(i),'xtick',[0 1000 2000 3000 4000],'xticklabel',{'0','1s','2s','3s','4s'})
    title(subplot1(i),files(i).name);
    
    % save pupil data, saccade data and optic flow duration
    save([fileSaveName files(i).name],'pupil_plot_set','baseline','trialInfo')
    
    %% data pre-processing for micro-saccade
    degree_plot_set = nan(length(saccade_eyedata),max(max(cell2mat(cellfun(@size,saccade_eyedata,'UniformOutput',false)))));
    plot_x_S = (1: max(cell2mat(cellfun(@size,saccade_eyedata,'UniformOutput',0)))) * dt;
    
    microsaddade_pair = cell(length(saccade_eyedata),1);
    microsaccade_meanV = cell(length(saccade_eyedata),1);
    microsaccade_meanAMP = cell(length(saccade_eyedata),1);

    
    for j = 1 : length(saccade_eyedata) % 可以改为 for j = find(trial_choice == 1) 之类
%     for j = find(trial_choice == 1)
        trial_Seyedata = saccade_eyedata{j};
        smooth_level = 1000 / mean(diff(trial_Seyedata(:,1))) / smooth_Hz;
        if ~isempty(trial_Seyedata)            
            % check for data validation
            if sum(isnan(trial_Seyedata(:,2))) > size(trial_Seyedata,1) / 5
                disp(['Trial ' num2str(j) ' in ' files(i).name ' skipped.'])
                continue
            end
        else
            disp(['Something wrong here, no data for trial ' num2str(j) ' in ' file(i).name ])
            continue
        end
            
        trial_Seyedata(:,2) = smoothdata(trial_Seyedata(:,2),'gaussian',5,'includenan');
        trial_Seyedata(:,3) = smoothdata(trial_Seyedata(:,3),'gaussian',5,'includenan');
        
        % calculate for microsaccade
        micro_sac_index = findMicrosaccade(trial_Seyedata);
        if ~isempty(micro_sac_index)
            if i == 1
                micro_pop_set = zeros(1,length(trial_Seyedata)+1000);
                micro_pop_set(micro_sac_index(:,1)) = 1;
            else
                micro_pop_set(micro_sac_index(:,1)) = micro_pop_set(micro_sac_index(:,1)) + 1;
            end
        end
        if ~isempty(micro_sac_index)
            for k = 1:size(micro_sac_index,1)
                axes(subplot2(i));
                hold on
                plot(trial_Seyedata(micro_sac_index(k,1):micro_sac_index(k,2),2),trial_Seyedata(micro_sac_index(k,1):micro_sac_index(k,2),3),'color',[0.8/k 0.8/j 0.8]);
            end
        end
    end
    axes(subplot2(i));
    hold on
    plot(SCREEN.screenWidth/2,SCREEN.screenHeight/2,'g*')
    view_window = 200;
    set(subplot2(i),'ylim',[SCREEN.screenHeight/2 - view_window , SCREEN.screenHeight/2 + view_window],'xlim',[SCREEN.screenWidth/2 - view_window,  SCREEN.screenWidth/2 + view_window]); 

    % 画1°的圆
    circle_r1 = calculate_r(1); 
    plot_circle(SCREEN.screenWidth/2,SCREEN.screenHeight/2,circle_r1,[0 0 0],subplot2(i))

    % 画3°的圆
    circle_r2 = calculate_r(3);
    plot_circle(SCREEN.screenWidth/2,SCREEN.screenHeight/2,circle_r2,[0.5 0.5 0.5],subplot2(i))
    title(subplot2(i),files(i).name);
        
        
    %% data pre-processing for saccade
    degree_plot_set = nan(length(saccade_eyedata),max(max(cell2mat(cellfun(@size,saccade_eyedata,'UniformOutput',false)))));
    plot_x_S = (1: max(cell2mat(cellfun(@size,saccade_eyedata,'UniformOutput',0)))) * dt;
    
    microsaddade_pair = cell(length(saccade_eyedata),1);
    microsaccade_meanV = cell(length(saccade_eyedata),1);
    microsaccade_meanAMP = cell(length(saccade_eyedata),1);
    
    for j = 1 : length(saccade_eyedata)
        trial_Seyedata = saccade_eyedata{j};
        smooth_level = 1000 / mean(diff(trial_Seyedata(:,1))) / smooth_Hz;
        if ~isempty(trial_Seyedata)            
            % check for data validation
            if sum(isnan(trial_Seyedata(:,2))) > size(trial_Seyedata,1) / 5
                disp(['Trial ' num2str(j) ' in ' files(i).name ' skipped.'])
                continue
            end
        else
            disp(['Something wrong here, no data for trial ' num2str(j) ' in ' file(i).name ])
            continue
        end
        
        % count for blink times
%         [blinkSNUM(j,i),~] = size(blinktimesj);
        
        % check for view window validation
%         del_index = or(or(or(Seyedataj_purified(:,2) < view_range(1) , Seyedataj_purified(:,2) > view_range(2) ),Seyedataj_purified(:,3) < view_range(3)),Seyedataj_purified(:,3) > view_range(4));
%         Seyedataj_purified(del_index,2:4) = NaN;
        
        % convert x-y pixel data to degree/s
        degree_saccade = pixel2degreev(trial_Seyedata,1,2,3);
        
        % prepare for plot
        degree_plot_set(j,1:length(trial_Seyedata(:,4))) = smooth(degree_saccade(:,2),smooth_level)';
        
        % calculate for saccade
        [microsaddade_pair{j},microsaccade_meanV{j},microsaccade_meanAMP{j}] = findSaccade(degree_saccade,saccade_thres,saccade_thres_cap);
    end
    
    axes(subplot3(i));
%     shadedErrorBar(plot_x_S,degree_plot_set,{@nanmean,@nanstd},'lineprops', '-b');
%     shadedErrorBar([],degree_plot_set,{@nanmedian,@(plot_x_S) nanstd(plot_x_S) / sqrt(size(degree_plot_set,1))*1.96},'lineprops',{'b-'});
%     hold on
%     set(subplot2(i),'ytick',[]);

    
%     figure(10);clf
    micro_pop = microsaccade_rate_cal(micro_pop_set,50,10);
    plot(micro_pop(:,1)*10,micro_pop(:,2),'k');
%     title('Microsaccade rate')
    title(subplot3(i),files(i).name);  
    
    clear eyedata startTime endTime trialInfo pupil_plot_set degree_plot_set plot_x
end
end

function [] = plot_circle( circle_x,circle_y,circle_r,circle_color,figure_handle)
theta = 0:0.1:2*pi;
circle1 = circle_x + circle_r * cos(theta);
circle2 = circle_y + circle_r * sin(theta);
axes(figure_handle);
hold on
plot(circle1,circle2,'--','Color',circle_color);
axis equal
end

function circle_r = calculate_r(circle_degree)
global SCREEN
circle_r = tand(circle_degree) * SCREEN.view_distance * SCREEN.screenWidth / SCREEN.screenWidth_real;
end

function micro_rate = microsaccade_rate_cal(microNUM,bin_width,step_length)
micro_rate = [];
for i = 1: floor(length(microNUM) / step_length ) - ceil(bin_width / step_length)
    if  i == 1
        micro_rate = cat(1,micro_rate,[i nansum(microNUM(1:bin_width)) / bin_width]);
    else
        micro_rate = cat(1,micro_rate,[i nansum(microNUM( step_length*i : step_length*i + bin_width )) / bin_width]);
    end
end
end

function y =nanse(input)
y=nanstd(input)*1.96/sqrt(sum(~isnan(input)));
end
