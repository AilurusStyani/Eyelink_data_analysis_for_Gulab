function [x_v,y_v] = pixel2velocity_ind(pixelData,x_col,y_col,t_col)
% This function converts pixel data to velocity(degree/s) data for individual x-axis and y-axis.
% x_col marks for the column for x-axis, y_col marks for y-axis and t_col for the time(ms)
% gives output data x_v for velocity in x-axis and y_v for velocity in y-axis.
%
%
% BYC at 22 OCT 2018

global SCREEN;

if ~exist('SCREEN.width','var')
    SCREEN.width = 1280; % pixel
end
if ~exist('SCREEN.height','var')
    SCREEN.height = 1024; % pixel
end
if ~exist('SCREEN.width_real','var')
    SCREEN.width_real = 37.5; % cm
end
if ~exist('SCREEN.height_real','var')
    SCREEN.height_real = 30; % cm
end
if ~exist('SCREEN.distance','var')
    SCREEN.distance = 60; % cm
end

central_point = [SCREEN.width/2 SCREEN.height/2];

dt = pixelData(2:end,t_col)-pixelData(1:end-1,t_col);


