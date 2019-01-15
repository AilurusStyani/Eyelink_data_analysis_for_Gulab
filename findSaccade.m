function [saccade_pair,saccade_meanV,mean_amplitude] = findSaccade(degree_data,saccade_thres,saccade_cap)
% this function can find saccades on the threshold

sac_bin = find(degree_data(:,2) >= saccade_thres);
sac_bin_r = rot90(sac_bin,2);
if ~isempty(sac_bin)
    sac_pair = sac_bin(diff(sac_bin) ~= 1);
    sac_pair = cat(1,sac_pair,sac_bin_r(diff(sac_bin_r) ~= -1 ));
    sac_pair = cat(1,sac_pair,sac_bin(1));
    sac_pair = sort(sac_pair);
else
    saccade_pair = [];
    saccade_meanV = [];
    mean_amplitude = [];
    return
end

saccade_pair = [];
mean_index = [];
amplitude = [];
del_set = [];
if exist('saccade_cap','var')
    [~,del_index] = findpeaks(degree_data(:,2),'MinPeakHeight',saccade_cap,'MinPeakDistance',50);
    for i = 1:length(del_index)
        del_set = cat(1,del_set,find(sac_pair > del_index(i),1));
        del_set = cat(1,del_set,find(sac_pair < del_index(i),1,'last'));
    end
    unique(del_set);
    sac_pair(del_set) = [];
end

for i = 1 : length(sac_pair) - 1
    if all(degree_data(sac_pair(i):sac_pair(i+1),2) > saccade_thres) && sac_pair(i+1) - sac_pair(i) > 8 % saccade maintained for at least 4ms
        saccade_pair = cat(1,saccade_pair,[ degree_data(sac_pair(i),1) degree_data(sac_pair(i+1),1) ] );
        amplitude = cat(1,amplitude,sum( degree_data(sac_pair(i):sac_pair(i+1),2) ));
        mean_index = cat(1,mean_index,[sac_pair(i):sac_pair(i+1)]');
    end
end

saccade_meanV = nanmean(degree_data(mean_index,2));
mean_amplitude = nanmean(amplitude)/1000;