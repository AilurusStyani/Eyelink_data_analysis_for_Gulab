function degree_data=pixel2degreev(eyedata,col_t,col_x,col_y)
% this function can convert pixel data of x/y to saccade degree in every ms
%
% BYC Sep 2018

screen_h = 30; % cm
screen_hp = 1024; % pixel
screen_w = 37.5; % cm
screen_wp = 1280; % pixel
view_dis = 60; % cm

if ~exist('clo_t','var') || isempty(col_t); col_t = 1; end
if ~exist('clo_x','var') || isempty(clo_x); col_x = 2; end
if ~exist('clo_y','var') || isempty(clo_y); col_y = 3; end

time = eyedata(:,col_t);
x_pixel = eyedata(:,col_x);
y_pixel = eyedata(:,col_y);

v_d = [];
for i = 1 : length(time)
    
    if i == 1
        v_d = cat(1,v_d,0);
        continue
    end
    
    timed = time(i)-time(i-1);
    
    % covert all the pixel to real distance
    ni2c = sqrt(((x_pixel(i)-screen_wp/2)/screen_wp*screen_w)^2 + ((y_pixel(i)-screen_hp/2)/screen_hp*screen_h)^2); % present point to screen center
    pi2c = sqrt(((x_pixel(i-1)-screen_wp/2)/screen_wp*screen_w)^2 + ((y_pixel(i-1)-screen_hp/2)/screen_hp*screen_h)^2);
    pi2ni = sqrt(((x_pixel(i)-x_pixel(i-1))/screen_wp*screen_w)^2 + ((y_pixel(i)-y_pixel(i-1))/screen_hp*screen_h)^2);
    eye2ni = sqrt(ni2c^2 + view_dis^2);
    eye2pi = sqrt(pi2c^2 + view_dis^2);
    degreei = acosd((eye2ni^2 + eye2pi^2 - pi2ni^2)/(2 * eye2ni * eye2pi));

    v_d = cat(1,v_d,degreei/timed*1000);
end

degree_data = eyedata;
degree_data(:,col_x) = v_d;
degree_data(:,col_y) = [];
end
    