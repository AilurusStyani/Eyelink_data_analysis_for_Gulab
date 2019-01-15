function [purifiedData,blinktimes,errorflag,purify_times] = BlinkNoisePurify_NaN(importdata,dt,error_set,pupilcol)
% This function can treat blink data based on pupil size.
% the first import is for eye data in column, the last import indicates
% which column is for pupil size.
% Or column 1 for time, column 2 for x pixel, column 3 for y pixel,
% column 4 for pupil size in default.
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
% 8. have a minimum of pupil size when the blink detected, which probably be the semi-blink;
% 9. still have unknow rifts in x and y pixel, that should be caused
%     by not regorous operation or other unknown conditions.
% 10.Still lots of points outside the screen. It may be caused by 
%     not correctly calibration & validation, or may be the subject squinting for too long (>1s)
% 11.Still have unknown noisy point (sudden decrease to zero for ~1
%     to 20ms, discontinuously), this maybe caused by Eyelink.
% 12.Input data too short
%
% BY BYC SEP/2018

if ~exist('pupilcol','var')
    pupilcol = 4;
end

purify_times = 0;
blinktimes = [];
errorflag = [];

if median(importdata(:,pupilcol)) == 0
    errorflag = 7;
    purifiedData = [];
    blinktimes = [];
    purify_times = [];
    return
end
    
for i=1:100 % check 100 times as maximum, more than which number should be consider as unreliable
    if i == 1
        [blink_time_data,blink_timesi,Flagi] = blink_treatment_pupil1(importdata,pupilcol,dt);
        if ~isempty(Flagi)
            errorflag = cat(1,errorflag,Flagi);
            errorflag = unique(errorflag);
        end
        if ~isempty(blink_timesi)
            blinktimes = cat(1,blinktimes,blink_timesi);
        end
        if ~isempty(blink_time_data)
            [~,delete_index,~] = intersect(importdata(:,1),blink_time_data);
            checked_data = importdata;
            checked_data(delete_index,:) = [];
            purify_times = purify_times + 1;
        else
            break
        end
    elseif ~isempty(checked_data) || isempty(intersect(errorflag,error_set))
        [blink_time_datai,blink_timesi,Flagi] = blink_treatment_pupil1(checked_data,pupilcol,dt);
        
%         % plot for debug
%         figure(2000);clf;
%         plot(inputdata(:,1),inputdata(:,3),'r')
%         hold on
%         plot(inputdata(:,1),inputdata(:,2),'g')
%         plot(checked_data(:,1),checked_data(:,2),'*g')
%         plot(inputdata(:,1),inputdata(:,4),'b')
%         % end

        if ~isempty(Flagi)
            errorflag = cat(1,errorflag,Flagi);
            errorflag = unique(errorflag);
        end
        if ~isempty(blink_timesi)
            blinktimes = cat(1,blinktimes,blink_timesi);
        end
        if ~isempty(blink_time_datai)
            blink_time_data = cat(1,blink_time_data,blink_time_datai);
            [~,delet_indexi,~] = intersect(checked_data(:,1),blink_time_datai);
            checked_data(delet_indexi,:) = [];
            purify_times = purify_times + 1;
        else
            break
        end
    end
end
blink_time_data = unique(blink_time_data);
[~,delete_index,~] = intersect(importdata,blink_time_data);
purifiedData = importdata;
purifiedData(delete_index,2:4) = NaN;

if size(purifiedData,1) < size(importdata,1) /5*4
    errorflag = cat(1,errorflag,7);
    unique(errorflag);
end
errorflag=errorflag';
end