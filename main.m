
function main
% By BYC at 2018/10/08
close all
clear all
warning off

file_path = 'D:\ION\2018_rotation\CJW\data0105\heading';
subject{1} = 'CJW';
subject{2} = 'YZX';
subject{2} = 'PPY';
subject{3} = 'DQS';

no_choice_version = 0; % 1 for no choice version; other number for choise version
error_set = [7 8 9 10 12]; % for more detail, please read BlinkNoisePurify_NaN.m
maxtrialNUM = 400; % the maximum trial numbers in all blocks
fixation_degree = 3; % degree measure as an angle ( not radian ), the data outside are invalide as NaN
saccade_thres = 3; % degree/s
saccade_thres_cap = 100; % degree/s
micro_rate_bin_width = 50; % ms
micro_rate_step_length = 5; % ms
mic_sac_dir_group_num = 8; % how many different degree groups you wanna to assort with?

edf2asc = 1; % 0:convert data from EDF files to asc files; 1: skip this process;
data_convert = 1; % 0:convert data to mat files; 1: skip this process;

exeFilePath = 'D:\ION\2018_rotation\CJW\analyses\edf2asc'; % Where you put the 'edf2asc.exe'
fileSaveName = 'PupilAdapt2TrialEnd_'; % change prefix based on the period you are interested

%% calculation for fixation windows
global SCREEN
SCREEN.screenWidth = 1280; % pixel
SCREEN.screenHeight = 1024; % pixel
SCREEN.view_distance = 60; % cm
SCREEN.screenWidth_real = 37.5; % cm
SCREEN.screenHeight_real = 30; % cm
checkWindow_degree = 1.5; % °
dx = tand(checkWindow_degree) * SCREEN.view_distance / SCREEN.screenWidth_real * SCREEN.screenWidth;
dy = tand(checkWindow_degree) * SCREEN.view_distance / SCREEN.screenHeight_real * SCREEN.screenHeight;
view_range = [SCREEN.screenWidth/2-dx SCREEN.screenWidth/2+dx SCREEN.screenHeight/2-dy SCREEN.screenHeight/2+dy];

    
%% set config for subplots
set(0,'defaultfigurecolor','w');
set(figure(1),'pos',[20 20 1849 900],'Name','Pupil size');clf;
subplot1 = tight_subplot(2,3,[0.1 0.035],[0.1 0.1]);

set(0,'defaultfigurecolor','w');
set(figure(2),'pos',[20 20 1849 900],'Name','Microsaccade dots');clf;
subplot2 = tight_subplot(2,3,[0.1 0.035],[0.1 0.1]);

set(0,'defaultfigurecolor','w');
set(figure(3),'pos',[20 20 1849 900],'Name','Microsaccade rate');clf;
subplot3 = tight_subplot(2,3,[0.1 0.035],[0.1 0.1]);

set(0,'defaultfigurecolor','w');
set(figure(4),'pos',[20 20 1849 900],'Name','Microsaccade direction');clf;
subplot4 = tight_subplot(2,3,[0.1 0.035],[0.1 0.1]);

saccade_maxV = [];
saccade_AMP = [];

    
for sub = 1:length(subject)
    data_path = fullfile(file_path,subject{sub});


    %% convert EDF files to ASC files, it may takes a long time
    if edf2asc == 0
        edf2asc_checkasc(data_path,exeFilePath);
    end

    %% extract eyedata trial by trial from the whole block
    if data_convert == 0
        getEyeDatas(data_path,error_set,no_choice_version);
    end

    % read the converted MAT files
    cd(data_path)
    files = dir('converted_*');
    fileNUM = length(files);
    blinkNUM = NaN(maxtrialNUM,fileNUM);

    % combine files from different blocks
    combinedSac = [];
    combined_dm = [];
    combinedSC = [];
    combined_choice = [];
    combinedSO = [];
    combinedST = [];
    combined_head = [];
    combined_pupil = [];
    combined_baseline = [];
    trial_no_micro = 0;
        
    for i = 1 : fileNUM % for each block
        load(files(i).name)
        trialInfo = strrep(files(i).name,'converted_','');
        load(trialInfo)
        combinedSac = [combinedSac;saccade_eyedata];
        combined_dm = [combined_dm;decision_making];
        combinedSC = [combinedSC; start_choice];
        combinedSO = [combinedSO;stimulus_on];
        combinedST = [combinedST;stimulus_term];
        combined_choice = [combined_choice;choice'];
        combined_head = [combined_head;head'];
        combined_pupil = [combined_pupil;pupil_eyedata];
        combined_baseline = [combined_baseline;baseline];
    end

    %%% seperate the trial into different categories
    trialhead = unique(combined_head);
    easy_head = trialhead( [1:2 , end-2:end] ); % defind easy trials
    hard_head = trialhead( length(trialhead) /2 - 1 : length(trialhead) /2 + 2 ); % defind hard trials
    medium_head = trialhead( ~ismember(trialhead, [easy_head;hard_head] ));
    
    % medium_head1 = trialhead(find(trialhead>=91 & trialhead<180));
    % medium_head2 = trialhead(find(trialhead<=89 & trialhead>0));
%     medium_head = [medium_head1;medium_head2];


    result_temp = [];
    axes(subplot1(sub));
    hold on
    axes(subplot2(sub));
    hold on
    axes(subplot3(sub));
    hold on
    axes(subplot4(sub));
    hold on
    axis equal

    for d = 1:3
        
        if d == 1
            [Lia,Locb] = ismember(combined_head,hard_head); % for hard tasks
            colormap = [220 83 86]/255;
        elseif d==2
            [Lia,Locb] = ismember(combined_head,easy_head); % for easy tasks
            colormap = [142 195 167]/255;
        else
            [Lia,Locb] = ismember(combined_head,medium_head); % for medium difficulty tasks
            colormap = [240 203 105]/255;
        end
        
        % reset each vaule based on different difficulties
        saccade_eyedata = combinedSac(Lia,:);
        decision_making = combined_dm(Lia,:);
        start_choice = combinedSC(Lia,:);
        trial_choice = combined_choice(Lia,:);
        stimulus_on = combinedSO(Lia,:);
        stimulus_term = combinedST(Lia,:);
        trial_head = combined_head(Lia,:);
        pupil_eyedata = combined_pupil(Lia,:);
        baseline = combined_baseline(Lia,:);
        degree_group = -pi : 2*pi/mic_sac_dir_group_num : pi;
        trial_num_dir = 0;
        
        pupil_plot_set = nan(length(pupil_eyedata),max(max(cell2mat(cellfun(@size,pupil_eyedata,'UniformOutput',false)))));%生成NaN矩阵(trial数,最多项的列数)
        plot_x_P = (1: max(cell2mat(cellfun(@size,pupil_eyedata,'UniformOutput',0)))) * dt;%时间轴
    
%         % plot for debug
%         test = cell2mat(pupil_eyedata);
%         if ~isempty(test)
%             test_plot = figure(50);
%             clf(test_plot);
%             
%             plot(test(:,4));
%             hold on
%             plot(test(:,2),'r');
%             plot(test(:,3),'g');
%         end
    
        %% data pre-processing for pupil data
        for j = 1 : sum(Lia) % for each trial
            trial_Peyedata = pupil_eyedata{j};
        
            soIndex = find( trial_Peyedata(:,1) == stimulus_on(j,1),1);
            stIndex = find( trial_Peyedata(:,1)==stimulus_term(j,1),1);
            scIndex = find( trial_Peyedata(:,1)==start_choice(j,1),1);
            dmIndex = find( trial_Peyedata(:,1)==decision_making(j,1),1);
            
    %%% select the period of time for analyses       
%             trial_Peyedata = trial_Peyedata((dcIndex-x):(dcIndex+500),:);
%             trial_Peyedata = trial_Peyedata((soIndex):(scIndex+500),:);
            trial_Peyedata = trial_Peyedata((soIndex + 0):(scIndex + 500),:);
            
%             % plot for test
%             figure(500);clf;
%             plot(trial_Peyedata(:,4),'b');
%             hold on
%             plot(trial_Peyedata(:,2),'r');
%             plot(trial_Peyedata(:,3),'g');
%             % test end
            
            if isempty(trial_Peyedata)      
                disp(['Something wrong here, no data for trial ' num2str(j) ' in ' files(i).name ])
                continue
            else
                
                % elimate blink data by NaN value
                [Peyedataj_purified,~,errorflagj,purifyNUMj] = BlinkNoisePurify_NaN(trial_Peyedata,dt,error_set,4);
               if isempty(Peyedataj_purified)   
                fprintf(2,['Trial ' num2str(j) ' in ' files(i).name ' pupil calculation skipped by too much NaN. \n'])
                continue
               end

                % check for Nan data
               if sum(isnan(Peyedataj_purified(:,2))) >= size(Peyedataj_purified,1) / 5
                    fprintf(2,['Trial ' num2str(j) ' in ' files(i).name ' pupil calculation skipped by too much NaN. \n'])
                    continue
               end

                % check for data validation
                if ~isempty(intersect(errorflagj,error_set))  || purifyNUMj > 3
                    if intersect(errorflagj,error_set) == 7
                        fprintf(2,['Trial ' num2str(j) ' in ' files(i).name ' pupil calculation skipped by eye blink. \n'])
                        continue
                    else
                        disp(['Trial ' num2str(j) ' in ' files(i).name ' excluded by error ' num2str(intersect(errorflagj,error_set))])
                        continue
                    end
                end
            end
        
            % prepare for plot
            if ~isnan(baseline(j))
    %             pupil_plot_set(j,1:length(trial_Peyedata(:,4))) = smooth(trial_Peyedata(:,4),smooth_level)'-baseline(j);
                pupil_plot_set(j,1:length(Peyedataj_purified(:,4))) = ( Peyedataj_purified(:,4)'-baseline(j) ) / baseline(j);
            end         
        end % trial end
    
        %% plot graphs with shaded error bar for pupil size

        if isempty(baseline)
            disp([subject{sub} ' is skipped.'])
            continue
        end
        axes(subplot1(sub));
        pupil_mean = nanmean(pupil_plot_set,1);
        pupil_se = nanse(pupil_plot_set,1);
        shadedErrorBar(plot_x_P,pupil_mean,pupil_se,'lineprops', {'-','color',colormap});
%         shadedErrorBar(plot_x_P,pupil_plot_set,{@nanmedian,@(plot_x_P)nanstd(plot_x_P)/ sqrt(sum(~isnan(pupil_plot_set(:,1))))},'lineprops', {'-','color',colormap});
%     yticks('auto');
%     xticks('auto');
%     set(subplot1(i),'ytick',[-500 -250 0 250 500 1000],'yticklabel',{'-500','-250','0','250','500','1000'},'ylim',[-750 1000]);
        set(subplot1(sub),'ytick',[-0.05 0 0.2],'yticklabel',{'-0.05','0','0.2'},'ylim',[-0.05 0.2]);
        set(subplot1(sub),'xtick',[0 1000 2000 3000 4000],'xticklabel',{'0','1s','2s','3s','4s'})
        title(subplot1(sub),subject{sub});
    
%         % save pupil data, saccade data and optic flow duration
%         save([fileSaveName subject],'pupil_plot_set','baseline','trialInfo')
    
        
        %% data pre-processing for saccade analyses
        degree_plot_set = nan(length(saccade_eyedata),max(max(cell2mat(cellfun(@size,saccade_eyedata,'UniformOutput',false)))));
        plot_x_S = (1: max(cell2mat(cellfun(@size,saccade_eyedata,'UniformOutput',0)))) * dt;

        Saccade_pair = cell(length(saccade_eyedata),1);
        acceleration_sac_xy = cell(length(saccade_eyedata),1);
        acceleration_sac = cell(length(saccade_eyedata),1);
        degrees_sac_xy = cell(length(saccade_eyedata),1);
        meanCI_S = [];


        % microsaccade_rate = zeros(1,(1300+x+1200));
        microsaccade_trialNUM=0;
        microsaccade_rate = zeros(1,max(max(cell2mat(cellfun(@size,saccade_eyedata,'UniformOutput',false)))));
%         microsaccade_rate = zeros(1,1600);
        
        if (mic_sac_dir_group_num) < 4
            error('It should have at least 4 direction to analyse micro-saccade directions.')
        end
        
        dir_assort = zeros(mic_sac_dir_group_num + 1,1);

        
        for j =1 : sum(Lia)
            far_choice = 0;
            near_choice = 0;
            trial_Seyedata = saccade_eyedata{j};
            if isempty(saccade_eyedata{j})
                 disp(['Something wrong here, no data for trial ' num2str(j) ' in ' files(i).name ])
                continue
            end

            soIndex = find(trial_Seyedata(:,1)==stimulus_on(j,1),1);
            stIndex = find(trial_Seyedata(:,1)==stimulus_term(j,1),1);
            scIndex = find(trial_Seyedata(:,1)==start_choice(j,1),1);
            dmIndex = find(trial_Seyedata(:,1)==decision_making(j,1),1);
            
    %%% select the period of time for analyses       
%             trial_Seyedata = trial_Seyedata((dcIndex-x):(dcIndex+500),:);
            trial_Seyedata = trial_Seyedata((soIndex):(scIndex+500),:);
%             trial_Seyedata = trial_Seyedata((soIndex):(dmIndex),:);
%             trial_Seyedata = trial_Seyedata((dmIndex-1500):(dmIndex+10),:);
            
%             % plot for test
%             figure(6);clf;
%             plot(trial_Seyedata(:,4),'b');
%             hold on
%             plot(trial_Seyedata(:,2),'r');
%             plot(trial_Seyedata(:,3),'g');
%             % test end
            
            if isempty(trial_Seyedata)      
                disp(['Something wrong here, no data for trial ' num2str(j) ' in ' files(i).name ])
                continue
            else
                
                % elimate blink data by NaN value
                [Seyedataj_purified,~,errorflagj,purifyNUMj] = BlinkNoisePurify_NaN(trial_Seyedata,dt,error_set,4);
               if isempty(Seyedataj_purified)   
                fprintf(2,['Trial ' num2str(j) ' in ' files(i).name ' saccade calculation skipped by too much NaN. \n'])
                continue
               end

                % check for Nan data
               if sum(isnan(Seyedataj_purified(:,2))) >= size(Seyedataj_purified,1) / 5
                    fprintf(2,['Trial ' num2str(j) ' in ' files(i).name ' saccade calculation skipped by too much NaN. \n'])
                    continue
               end

                % check for data validation
                if ~isempty(intersect(errorflagj,error_set))  || purifyNUMj > 3
                    if intersect(errorflagj,error_set) == 7
                        fprintf(2,['Trial ' num2str(j) ' in ' files(i).name ' saccade calculation skipped by eye blink. \n'])
                        continue
                    else
                        disp(['Trial ' num2str(j) ' in ' files(i).name ' excluded by error ' num2str(intersect(errorflagj,error_set))])
                        continue
                    end
                end
            end

            % convert x-y pixel data to degree/s
            degree_saccade = pixel2degreexy(Seyedataj_purified,1,2,3); % convert x-y pixel data to linear velocity data

            % detection micro-saccade
            microSaccade_index = findMicrosaccade_ellipse(degree_saccade);
            
            % 每个microsaccade起止时序
            if ~isempty(microSaccade_index)
                trial_num_dir = trial_num_dir +1;
                plot(subplot2(sub),microSaccade_index(:,1), ones(length(microSaccade_index(:,1)),1)*5+j,'.k');
                set(subplot2(sub),'YTickLabelMode','auto','XTickLabelMode','auto');
                title(subplot2(sub),subject{sub});
                microsaccade_rate(1,microSaccade_index(:,1)) = microsaccade_rate(1,microSaccade_index(:,1))+1;
                microsaccade_trialNUM = microsaccade_trialNUM+1;
                
                saccade_AMP = cat(1,saccade_AMP,atand(sqrt(power(Seyedataj_purified(microSaccade_index(:,2),2) - Seyedataj_purified(microSaccade_index(:,1),2),2) ...
                    + power(Seyedataj_purified(microSaccade_index(:,2),3) - Seyedataj_purified(microSaccade_index(:,1),3),2)) ./ SCREEN.screenWidth .* SCREEN.screenWidth_real ./ SCREEN.view_distance));
                
                for k = 1:size(microSaccade_index,1)
                    microsac_eyedata = degree_saccade(microSaccade_index(k,1):microSaccade_index(k,2),2:3);
                    micro_center = ( degree_saccade(microSaccade_index(k,1),2:3) + degree_saccade(microSaccade_index(k,2),2:3) ) / 2;
                    delv2c = sqrt(power(microsac_eyedata(:,1) - micro_center(1) , 2) + power(microsac_eyedata(:,2) - micro_center(2) , 2));
                    max_ampNUM = find(delv2c == max(delv2c),1);
                    
                    saccade_maxV = cat(1,saccade_maxV,max(sqrt(degree_saccade(microSaccade_index(k,1):microSaccade_index(k,2),2).^2 + degree_saccade(microSaccade_index(k,1):microSaccade_index(k,2),3).^2)));

    %           figure(120) % plot for debug
    %           plot(microsac_eyedata(:,1),microsac_eyedata(:,2));
    %           hold on
    %           plot(microsac_eyedata(max_ampNUM,1),microsac_eyedata(max_ampNUM,2),'r*');
    %           plot(micro_center(1),micro_center(2),'g*');
    
                % if k == size(microSaccade_index,1)%the last microsaccade
                    if microsac_eyedata(max_ampNUM,1) > micro_center(1)
                        far_choice = far_choice+1;   
                    else
                        near_choice = near_choice+1;
                    end
                end
    %         result_temp = [result_temp;far near trial_head(j) trial_choice(j) d];
    
            % preparing for drawing micro-saccade direction
                dir_set = [trial_Seyedata(microSaccade_index(:,2),2) - trial_Seyedata(microSaccade_index(:,1),2), trial_Seyedata(microSaccade_index(:,2),3) - trial_Seyedata(microSaccade_index(:,1),3)];
                dir_degree = atan2(dir_set(:,2),dir_set(:,1));
                dir_abs_set = abs(degree_group - dir_degree);
                [~,dir_assort_index] = min(dir_abs_set,[],2);
                
                if ~isempty(dir_assort_index)
                    for mic_dir = 1:length(dir_assort_index)
                        dir_assort(dir_assort_index(mic_dir)) = dir_assort(dir_assort_index(mic_dir)) + 1;
                    end
                end
            else
                trial_no_micro = trial_no_micro + 1;
            end
        end

        % plot the micro-saccade rate
        micro_rate_plot = microsaccade_rate_cal(microsaccade_rate,micro_rate_bin_width,micro_rate_step_length/dt,microsaccade_trialNUM);
        plot(subplot3(sub),(micro_rate_plot(:,1)*micro_rate_step_length*dt),micro_rate_plot(:,2),'color',colormap,'LineWidth',0.8);
%         plot(subplot3(sub),[1500 1500],[0 ceil(max(micro_rate_plot(:,2))/2)*2],'b','LineWidth',1);
    %     plot(subplot3(sub),[(1300+TIME(t)*1000) (1300+TIME(t)*1000)],[0 ceil(max(micro_rate_plot(:,2))/2)*2],'b','LineWidth',1);
        set(subplot3(sub),'YTickLabel',[0:0.5:ceil(max(micro_rate_plot(:,2))/2)*2],'YTick',[0:0.5:ceil(max(micro_rate_plot(:,2))/2)*2],'YLim',[0 ceil(max(micro_rate_plot(:,2))/2)*2]);
%         set(subplot3(sub),'XTick',[0:500:2000],'XTickLabel',[-1500:500:500]);
        set(subplot3(sub),'XTick',[0:1000:4000],'XTickLabel',{'1s','2s','3s','4s'});
    %      set(subplot3(f),'XTickLabelMode','auto','YTickLabelMode','auto');
        xlabel('time (ms)');
        ylabel('microsaccade rate (Hz)');
        set(subplot3(sub),'LineWidth',1);
        
        % plot for micro-saccade direction
        dir_plot_set = [dir_assort' .* cos(degree_group) ; dir_assort' .* sin(degree_group)];
        dir_plot_set_1st = dir_plot_set(:,1) + dir_plot_set(:,end);
        dir_plot_set1 = [dir_plot_set_1st,dir_plot_set(:,2:end-1),dir_plot_set_1st] ./ trial_num_dir;
        plot(subplot4(sub),dir_plot_set1(1,:),dir_plot_set1(2,:),'color',colormap,'LineWidth',0.8);
%         plot(subplot4(sub),0,0,'k*');
        for mic_dir = 1:size(dir_plot_set1,2)
            plot(subplot4(sub),[0,dir_plot_set1(1,mic_dir)],[0,dir_plot_set1(2,mic_dir)],'--k');
        end

        disp(['For subject' subject{sub} ',' num2str(trial_no_micro) ' trials were no micro-saccade detected.']);
        clear eyedata startTime endTime choice distance head time pupil_plot_set degree_plot_set plot_x



    end
    saveas(3,'4.tif');
    
    % csvwrite('microsaccade_direction6.csv',result_temp);
    figure(5);
    title('Max Velocity to Amplitude for micro-saccades');
    xlabel('Max Velocity');
    ylabel('Amplitude');
    hold on;
    plot(saccade_maxV,saccade_AMP,'b.');
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

function data_out = nanse(import,direction)
if direction == 1
    data_out = nanstd(import,direction) ./ sqrt(sum(~isnan(import(:,1))));
elseif direction == 2
    data_out = nanstd(import,direction) ./ sqrt(sum(~isnan(import(1,:))));
else
    error('invalid value of direction')
end
end
