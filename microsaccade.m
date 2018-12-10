data_path = 'D:\学习工具\ION\GuLab\distance_perception\data';
status = [];
cd(data_path);
fileName = dir('converted_*.mat'); % find the files you want
fileNum = length(fileName);
smooth_Hz = 60; % the aim value for data smooth / low pass
fixation_degree = 3; % fixation degree
saccade_thres = 30; % degree/s, 30 as example
saccade_thres_cap = 100; % degree/s, 100 as example
microsaccade_thres = 6; % times of S.D as threshold.
% scale = length(pupil_plot_set); % data length
for i=1:fileNum
    load(fileName(i).name);
%     xlswrite('grand_mean.xlsx',pupil_plot_set,j);
    status = cat(1,status,saccade_eyedata); % get your datas
%     trialInfo = strrep(files(i).name,'converted_','');
%     trialInfo = load(trialInfo);
%     trial_choice = trialInfo.choice == 1; % 1:close; 2: far  
    fclose('all');


end

set(figure(2),'pos',[27 63 1849 892],'Name','Saccade track');clf;
subplot2 = tight_subplot(4,4,[0.03 0.03]);
% set(subplot2,'ytick',[]);
suptitle('Saccade track');

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
else
    pixel_d = pixel_x;
end
% view_range = [screenWidth/2-pixel_x screenWidth/2+pixel_x screenHeight/2-pixel_y screenHeight/2+pixel_y];


%% data pre-processing for micro-saccade
        degree_plot_set = nan(length(status),max(max(cell2mat(cellfun(@size,status,'UniformOutput',false)))));
    plot_x_S = (1: max(cell2mat(cellfun(@size,status,'UniformOutput',0)))) * dt;
    
    microsaddade_pair = cell(length(status),1);
    microsaccade_meanV = cell(length(status),1);
    microsaccade_meanAMP = cell(length(status),1);
    axes(subplot2(i));
   
    for j = 1 : length(status) % 可以改为 for j = find(trial_choice == 1) 之类
%     for j = find(trial_choice == 1)
        trial_Seyedata = status{j};
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
        
        % calculate for microsaccade
        micro_sac_index = findMicrosaccade(trial_Seyedata);
 
        if ~isempty(micro_sac_index)
            hold on
            for k = 1:size(micro_sac_index,1)
                plot(trial_Seyedata(micro_sac_index(k,1):micro_sac_index(k,2),2),trial_Seyedata(micro_sac_index(k,1):micro_sac_index(k,2),3),'color',[0.8/k 0.8/j 0.8]);
            end
        end
    end
        plot(SCREEN.screenWidth/2,SCREEN.screenHeight/2,'g*')
        view_window = 200;
        set(figure,'pos',[27 63 1849 892],'Name','Saccade track');clf;
        set(subplot2(i),'ylim',[SCREEN.screenHeight/2 - view_window , SCREEN.screenHeight/2 + view_window],'xlim',[SCREEN.screenWidth/2 - view_window,  SCREEN.screenWidth/2 + view_window]); 
        
        % 画1°的圆
        circle_r1 = calculate_r(0.1); 
        plot_circle(SCREEN.screenWidth/2,SCREEN.screenHeight/2,circle_r1,[0 0 0])
        
        % 画3°的圆
        circle_r2 = calculate_r(1.2);
        plot_circle(SCREEN.screenWidth/2,SCREEN.screenHeight/2,circle_r2,[0.5 0.5 0.5])
%         title(plot,files(i).name);
        
        
%     %% data pre-processing for saccade
%     degree_plot_set = nan(length(status),max(max(cell2mat(cellfun(@size,status,'UniformOutput',false)))));
%     plot_x_S = (1: max(cell2mat(cellfun(@size,status,'UniformOutput',0)))) * dt;
%     
%     microsaddade_pair = cell(length(status),1);
%     microsaccade_meanV = cell(length(status),1);
%     microsaccade_meanAMP = cell(length(status),1);
%     
%     for j = 1 : length(status)
%         trial_Seyedata = status{j};
%         smooth_level = 1000 / mean(diff(trial_Seyedata(:,1))) / smooth_Hz;
%         if ~isempty(trial_Seyedata)            
%             % check for data validation
%             if sum(isnan(trial_Seyedata(:,2))) > size(trial_Seyedata,1) / 5
%                 disp(['Trial ' num2str(j) ' in ' files(i).name ' skipped.'])
%                 continue
%             end
%         else
%             disp(['Something wrong here, no data for trial ' num2str(j) ' in ' file(i).name ])
%             continue
%         end
%         
%         % count for blink times
%         [blinkSNUM(j,i),~] = size(blinktimesj);
%         
%         % check for view window validation
%         del_index = or(or(or(Seyedataj_purified(:,2) < view_range(1) , Seyedataj_purified(:,2) > view_range(2) ),Seyedataj_purified(:,3) < view_range(3)),Seyedataj_purified(:,3) > view_range(4));
%         Seyedataj_purified(del_index,2:4) = NaN;
%         
%         % convert x-y pixel data to degree/s
%         degree_saccade = pixel2degreev(trial_Seyedata,1,2,3);
%         
%         % prepare for plot
%         degree_plot_set(j,1:length(trial_Seyedata(:,4))) = smooth(degree_saccade(:,2),smooth_level)';
%         
%         % calculate for saccade
%         [microsaddade_pair{j},microsaccade_meanV{j},microsaccade_meanAMP{j}] = findSaccade(degree_saccade,saccade_thres,saccade_thres_cap);
%     end
%     
%     axes(subplot2(i));
%     shadedErrorBar(plot_x_S,degree_plot_set,{@nanmean,@nanstd},'lineprops', '-b');
%     hold on
% %     set(subplot2(i),'ytick',[]);
%     title(subplot2(i),files(i).name);
%     
    
    
    
    clear eyedata startTime endTime trialInfo pupil_plot_set degree_plot_set plot_x

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