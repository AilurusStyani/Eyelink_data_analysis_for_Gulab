function degree_data=pixel2degreev(eyedata,col_t,col_x,col_y)
% this function can convert pixel data of x/y to saccade degree in every ms
% eyedata shoule be in at least three column, first for time marker, second for x-axis position in pixel, third for y-axis
% position in pixel. col_t, col_x and col_y present which column is for time, x-axis and y-axis.
%
% BYC Sep 2018

% checking for screen parameter, 1280*1024 pixel and 37.5 * 30 cm as default. 
% The unit for height and width could differ from cm, but should be the same to eye-screen distance.
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

if ~exist('clo_t','var') || isempty(col_t); col_t = 1; end
if ~exist('clo_x','var') || isempty(clo_x); col_x = 2; end
if ~exist('clo_y','var') || isempty(clo_y); col_y = 3; end

dt = mode(eyedata(:,col_t));
index_length = length(eyedata(:,1));

if index_length <= 4
    error('Input data are too short or you may need to transpose the matrix.')
end

index_p = zeros(index_length,2);

for i = 1 : index_length
    if i == 1
        index_p(i,:) = [eyedata(1,col_x) eyedata(1,col_y)];
    elseif i == 2
        index_p(i,:) = [eyedata(1,col_x) eyedata(1,col_y)] + [eyedata(3,col_x) eyedata(3,col_y)] - 2 * [eyedata(2,col_x) eyedata(2,col_y)];
    elseif i < index_length - 2
        index_p(i,:) = [eyedata(i+2,col_x) eyedata(i+2,col_y)] + [eyedata(i-2,col_x) eyedata(i-2,col_y)] - 2 * [eyedata(i,col_x) eyedata(i,col_y)];
    elseif i == index_length - 1
        index_p(i,:) = [eyedata(i+1,col_x) eyedata(i+1,col_y)] + [eyedata(i-1,col_x) eyedata(i-1,col_y)] - 2 * [eyedata(i,col_x) eyedata(i,col_y)];
    elseif i == index_length
        index_p(i,:) = [eyedata(index_length,col_x) eyedata(index_length,col_y)];
    end
end

% convert x-y in pixel to x-y in real distance 
index_cm = [index_p(:,1) ./ screen_wp .* screen_w , index_p(:,2) ./ screen_hp .* screen_h];

% convert x-y distance to vector distance
index_vd = sqrt(power(index_cm(:,1),2) + power(index_cm(:,2),2));

% convert to velocity as degree/time
index_v = atand(index_vd) ./ dt;

[~,col_num] = size(eyedata);
unchanged_col = 1:col_num;
unchanged_col([col_t col_x col_y]) = [];

degree_data = [eyedata(:,1) , index_v, eyedata(:,unchanged_col)];
end
    