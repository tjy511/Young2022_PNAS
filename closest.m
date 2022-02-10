function [n,idx] = closest(vec,target)

% [n,idx] = closest(vec,target)
% Quick function to find the closest value in an array, given a target
% number. Target can be a single number or a linear array. 
% NOTE: Will return first number if there are multiple minima. 
% 
% TJ Young
% 31 August 2016

for ii = 1:length(target)
    if isnan(nanmean(target))
        n(ii) = NaN;
        idx(ii) = NaN; 
    else
    tmp = abs(vec-target(ii));
    [~,idx(ii)] = min(tmp); % Index location of closest value in array
    n(ii) = vec(idx(ii)); % Value of closest value in array
    end
end