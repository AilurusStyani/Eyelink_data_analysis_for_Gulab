function micro_rate = microsaccade_rate_cal(microNUM,bin_width,step_length,trial_num)
% coded by BYC at 2019/1/11
% please import data, bin width, step length and trial number

micro_rate = [];
for i = 1: floor(length(microNUM) ./ step_length ) - ceil(bin_width ./ step_length)
    if  i == 1
        micro_rate = cat(1,micro_rate,[i nansum(microNUM(1:bin_width)) ./ trial_num]);
    else
        micro_rate = cat(1,micro_rate,[i nansum(microNUM( step_length*i : step_length*i + bin_width )) ./ trial_num]);
    end
end
micro_rate(:,2) = micro_rate(:,2) ./bin_width .* 1000;
        
        
        