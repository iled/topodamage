function [x_low, y_low, x_high, y_high, ranges_low, ranges_high] = range_analysis(study_var, fixed_vars, log_names)
% [x_low, y_low, x_high, y_high, ranges_low, ranges_high] = range_analysis(study_var, fixed_vars, log_names)
% Compute the coordinates corresponding to the points where the values of
% the study variable are within a low and high range, and all of the fixed_vars
% are within the same range (low, median, or high).
% study_var is a GRIDobj
% fixed_vars is a cell array of GRIDobj's
% log_names (optional) is a cell array with the names of the variables that have
% a log-normal distribution. Then the percentile ranges are computed for those
% cases.
% x_low, y_low are the (x, y) coordinates of the points for low values of the
% study_var
% x_high, y_high are the (x, y) coordinates of the points for low values of the
% study_var
% ranges_low is a row vector that contains integer IDs that tell from which range
% the selected points come from, for low values of the study_var
% ranges_high is a row vector that contains integer IDs that tell from which range
% the selected points come from, for high values of the study_var
% ranges IDs: 1 = low, 2 = medium, 3 = high.
% Julio Caineta Nov 6 2017

%low percentile max
prcmax1 = 20;
%median percentile range
prcmin2 = 35;
prcmax2 = 65;
%high percentile min
prcmin3 = 80;

% default value for log_flag is all false
nfixed = numel(fixed_vars);
nflags = nfixed + 1;
log_flag = false(1, nflags);
if nargin > 2
    nlogs = numel(log_names);
    if nlogs
        var_names = cell(1, nflags);
        for i = 1:nflags
            if i == 1
                var_names{i} = study_var.name;
            else
                var_names{i} = fixed_vars{i - 1}.name;
            end
        end
        log_flag = ismember(var_names, log_names);
    end
end

% set up percentiles
ps = [prcmax1, prcmin2, prcmax2, prcmin3];
perc = zeros(nflags, 4);
for p = 1:nflags
    if p == 1
        x = study_var.Z(:);
    else
        x = fixed_vars{p - 1}.Z(:);
    end
    
    if log_flag(p)
        
        perc(p, :) = lognperc(x, ps);
    else
        perc(p, :) = prctile(x, ps);
    end
end

% low and high percentiles conditions for the study variable
study_low_cond = study_var.Z(:) < perc(1, 1);
study_high_cond = study_var.Z(:) > perc(1, 4);

%%%%% LOW %%%%%
% build condition based on low range percentiles for the fixed variables
fixed_cond = true(size(study_var.Z(:)));
for i = 1:nfixed
    fixed_cond = fixed_cond & (fixed_vars{i}.Z(:) < perc(i + 1, 1));
end

% find locations of low values of study_var where fixed_vars are in the low range
i_small = find(study_low_cond & fixed_cond);
[x_low1, y_low1] = ind2coord(study_var, i_small);
% find locations of high values of study_var where fixed_vars are in the low range
i_high = find(study_high_cond & fixed_cond);
[x_high1, y_high1] = ind2coord(study_var, i_high);

%%%%% MEDIAN %%%%%
% build condition based on median range percentiles for the fixed variables
fixed_cond = true(size(study_var.Z(:)));
for i = 1:nfixed
    fixed_cond = fixed_cond & (fixed_vars{i}.Z(:) > perc(i + 1, 2) & ...
        fixed_vars{i}.Z(:) < perc(i + 1, 3));
end

% find locations of low values of study_var where fixed_vars are close to median
i_small = find(study_low_cond & fixed_cond);
[x_low2, y_low2] = ind2coord(study_var, i_small);
% find locations of high values of study_var where fixed_vars are close to median
i_high = find(study_high_cond & fixed_cond);
[x_high2, y_high2] = ind2coord(study_var, i_high);

%%%%% HIGH %%%%%
% build condition based on high range percentiles for the fixed variables
fixed_cond = true(size(study_var.Z(:)));
for i = 1:nfixed
    fixed_cond = fixed_cond & (fixed_vars{i}.Z(:) > perc(i + 1, 4));
end

% find locations of low values of study_var where fixed_vars are in the high range
i_small = find(study_low_cond & fixed_cond);
[x_low3, y_low3] = ind2coord(study_var, i_small);
% find locations of high values of study_var where fixed_vars are in the high range
i_high = find(study_high_cond & fixed_cond);
[x_high3, y_high3] = ind2coord(study_var, i_high);

%%%%% TOTAL %%%%%
x_low = [x_low1; x_low2; x_low3];
y_low = [y_low1; y_low2; y_low3];
x_high = [x_high1; x_high2; x_high3];
y_high = [y_high1; y_high2; y_high3];

% output a code to identify from which range the points come from
% 1 = low, 2 = medium, 3 = high
ranges_low = [repmat(1, 1, numel(x_low1)), ...
    repmat(2, 1, numel(x_low2)), ...
    repmat(3, 1, numel(x_low3))];
ranges_high = [repmat(1, 1, numel(x_high1)), ...
    repmat(2, 1, numel(x_high2)), ...
    repmat(3, 1, numel(x_high3))];
end