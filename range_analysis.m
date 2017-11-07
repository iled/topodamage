function [x_low, y_low, x_high, y_high] = range_analysis(study_var, fixed_vars)
% [x_low, y_low, x_high, y_high] = range_analysis(study_var, fixed_vars)
% Compute the coordinates corresponding to the points where the values of
% the study variable are within a low and high range, and all of the fixed_vars
% are within the same range (low, median, or high).
% study_var is a GRIDobj
% fixed_vars is a cell array of GRIDobj's
% x_low, y_low are the (x, y) coordinates of the points for low values of the
% study_var
% x_high, y_high are the (x, y) coordinates of the points for low values of the
% study_var
% Julio Caineta Nov 6 2017

%low percentile max
prcmax1 = 20;
%median percentile range
prcmin2 = 35;
prcmax2 = 65;
%high percentile min
prcmin3 = 80;

% low and high percentiles for the study variable
study_p = prctile(study_var.Z(:), [prcmax1, prcmin3]);
study_low_cond = study_var.Z(:) < study_p(1);
study_high_cond = study_var.Z(:) > study_p(2);

%%%%% LOW %%%%%
% build condition based on low range percentiles for the fixed variables
fixed_cond = true(size(study_var.Z(:)));
for i = 1:numel(fixed_vars)
    fixed_p = prctile(fixed_vars{i}.Z(:), prcmax1);
    fixed_cond = fixed_cond & (fixed_vars{i}.Z(:) < fixed_p);
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
for i = 1:numel(fixed_vars)
    fixed_p = prctile(fixed_vars{i}.Z(:), [prcmin2, prcmax2]);
    fixed_cond = fixed_cond & (fixed_vars{i}.Z(:) > fixed_p(1) & ...
        fixed_vars{i}.Z(:) < fixed_p(2));
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
for i = 1:numel(fixed_vars)
    fixed_p = prctile(fixed_vars{i}.Z(:), prcmin3);
    fixed_cond = fixed_cond & (fixed_vars{i}.Z(:) > fixed_p);
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
end