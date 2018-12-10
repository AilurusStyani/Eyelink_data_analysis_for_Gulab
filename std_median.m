function output = std_median(input)
output = sqrt( nansum( power( abs(input - median(input,'omitnan')),2 )) / (length(input(~isnan(input)))-1) );