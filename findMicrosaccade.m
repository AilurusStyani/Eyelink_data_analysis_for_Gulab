function saccade_index = findMicrosaccade(eyedata)
max_degree = 1.8; % set the maximum degree for the microsaccade

dt = mode(diff(eyedata(:,1)));
saccade_index = [];

eyedata = pixel2degreev(eyedata,1,2,3);

%% chcek for 5 * SD (median)
index =  eyedata(:,2) > (std_median(eyedata(:,2)) * 5); 
sac_on = find(diff(index) == 1);
sac_term = find(diff(index) == -1);

%% check for 8 ms as minimun duration
if length(sac_on) >= length(sac_term)
    for i = 1:length(sac_term)
        if (sac_term(i) - sac_on(find(sac_on < sac_term(i),1,'last'))) * dt > 10
            if nansum(eyedata(sac_on(find(sac_on < sac_term(i),1,'last')):sac_term(i),2))/1000 < max_degree
                saccade_index = cat(1,saccade_index,[sac_on(find(sac_on < sac_term(i),1,'last')) sac_term(i)]);
            end
        end
    end
elseif length(sac_on) < length(sac_term)
    for i = 1 : length(sac_on)
        if (sac_term(find(sac_term > sac_on(i),1))) * dt > 10
            if nansum(eyedata(sac_on(i):sac_term(find(sac_term > sac_on(i),1)),2))/1000 < max_degree
                saccade_index = cat(1,saccade_index,[sac_on(i) sac_term(find(sac_term > sac_on(i),1))]);
            end
        end
    end
end
end

function output = std_median(input)
output = sqrt( sum( power( abs(input - median(input,'omitnan')),2 )) / (length(input)-1) );
end

% function pixel_distance = calculate_degree(degree)
% global SCREEN
% pixel_distance = tand(degree) * SCREEN.view_distance / SCREEN.screenWidth_real * SCREEN.screenWidth;
% end
