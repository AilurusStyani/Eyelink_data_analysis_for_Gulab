% function saccade_index = findMicrosaccade(eyedata)
function valid_microsac = findMicrosaccade_ellipse(eyedata)
% col: 1 time; 2 velocity; 3 x pixel; 4 y pixel;5 pupil size
max_degree = 1.5; % set the maximum degree for the microsaccade
microsaccade_interval = 20; % ms
dt = mode(diff(eyedata(:,1))) ./ 1000;
% saccade_index = [];
valid_microsac = [];
%% calculate for 5 * SD (median)
thres_x = nanstd_median(eyedata(:,2)) * 5;
thres_y = nanstd_median(eyedata(:,3)) * 5;

index =  power(eyedata(:,2),2) ./ thres_x^2 + power(eyedata(:,3),2) ./ thres_y^2 > 1;
sac_on = find(diff(index) == 1);
sac_term = find(diff(index) == -1);

sac_term_index = zeros(size(sac_on));
for i = 1:length(sac_on)
    sac_term_index(i) = sac_term(find(sac_term > sac_on(i),1));
end
sac_term_index(sac_term_index == 0) = [];
sac_term = sac_term_index;

sac_on_index = zeros(size(sac_term));
for i = 1:length(sac_term)
    sac_on_index(i) = sac_on(find(sac_on < sac_term(i),1,'last'));
end
sac_on_index(sac_on_index == 0) = [];
sac_on = sac_on_index;

sac_del_index = [];
for i = 2:length(sac_on)
    if sac_on(i) < sac_term(i-1) + microsaccade_interval
        sac_del_index = cat(1,sac_del_index,i);
    end
end

sac_on(sac_del_index) = [];
sac_term(sac_del_index - 1) =[];

% figure(52);clf;
% plot(eyedata(:,2));
% hold on
% plot([0 size(eyedata,1)], [(std_median(eyedata(:,2)) * 1.5),(std_median(eyedata(:,2)) * 1.5)],'r');

%% calculate for amplitude and check for 6ms as minimum duration
for i = 1:length(sac_on)
    if sac_term(i) - sac_on(i) >= 6
        degree = pixel2deg(eyedata(sac_on(i):sac_term(i),4:5),dt);
        if degree <= max_degree
            valid_microsac = cat(1,valid_microsac,[sac_on(i), sac_term(i)]);
        end  
    end
end
%% check for 8 ms as minimum duration
% if length(sac_on) >= length(sac_term)
%     for i = 1:length(sac_term)
%         if (sac_term(i) - sac_on(find(sac_on < sac_term(i),1,'last'))) * dt > 8
%             if nansum(eyedata(sac_on(find(sac_on < sac_term(i),1,'last')):sac_term(i),2))/1000 < max_degree
%                 saccade_index = cat(1,saccade_index,[sac_on(find(sac_on < sac_term(i),1,'last')) sac_term(i)]);
%             end
%         end
%     end
% elseif length(sac_on) < l
%     ength(sac_term)
%     for i = 1 : length(sac_on)
%         if (sac_term(find(sac_term > sac_on(i),1))) * dt > 8
%             if nansum(eyedata(sac_on(i):sac_term(find(sac_term > sac_on(i),1)),2))/1000 < max_degree
%                 saccade_index = cat(1,saccade_index,[sac_on(i) sac_term(find(sac_term > sac_on(i),1))]);
%             end
%         end
%     end
% end
end

function output = nanstd_median(input)
output = sqrt( nansum( power( abs(input - median(input,'omitnan')),2 )) / (sum(~isnan(input))-1) );
end

function degree = pixel2deg(eyedata,dt)
global SCREEN
if ~exist('SCREEN.height','var')
    screen_h = 30; % cm
else
    screen_h = SCREEN.height;
end

if ~exist('SCREEN.width','var')
    screen_w = 37.5; % cm
else
    screen_w = SCREEN.width;
end

if ~exist('SCREEN.height_pixel','var')
    screen_hp = 1024; % pixel
else
    screen_hp = SCREEN.height_pixel;
end

if ~exist('SCREEN.width_pixel','var')
    screen_wp = 1280; % pixel
else
    screen_wp = SCREEN.width_pixel;
end

% set the distance from eye to screen, 60 cm as default
if ~exist('SCREEN.viewdistance','var') && ~exist('SCREEN.distance','var')
    view_dis = 60; % cm
elseif exist('SCREEN.viewdistance','var')
    view_dis = SCREEN.viewdistance;
else
    view_dis = SCREEN.distance;
end
smooth_x = eyedata(:,1);
smooth_y = eyedata(:,2);

max_pixel_dis = max(sqrt(power(smooth_x - smooth_x(1),2) + power(smooth_y - smooth_y(1),2)));
degree = atand(max_pixel_dis ./ screen_wp .* screen_w ./ view_dis);

% point2np2 = power((smooth_x(2:end) - smooth_x(1:end-1)),2) + power((smooth_y(2:end) - smooth_y(1:end-1)),2);
% point2c = sqrt(power((smooth_x-screen_wp/2)./screen_wp.*screen_w,2) + power((smooth_y-screen_hp/2)./screen_hp.*screen_h,2)); % convert pixel to distance to center
% 
% degree = acosd( ( power(point2c(2:end),2) + view_dis^2 + power(point2c(1:end-1),2) + view_dis^2 - point2np2) ./ (2 .* sqrt( power(point2c(2:end),2) + view_dis^2) .* sqrt( power(point2c(1:end-1),2) + view_dis^2)) );
% degree = sqrt(power((smooth_x(2:end) - smooth_x(1:end-1)),2) + power((smooth_y(2:end) - smooth_y(1:end-1)),2)) ./ (tand(1) .* view_dis ./ screen_w .* screen_wp);
end
        
% function pixel_distance = calculate_degree(degree)
% global SCREEN
% pixel_distance = tand(degree) * SCREEN.view_distance / SCREEN.screenWidth_real * SCREEN.screenWidth;
% end