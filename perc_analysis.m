function [x_low, y_low, x_high, y_high] = perc_analysis(study_var, fixed_vars)

%low percentile max
prcmax1 = 20;
%median percentile range
prcmin2 = 35;
prcmax2 = 65;
%high percentile min
prcmin3 = 80;

% low and high percentiles for the study variable
study_p = prctile(study_var.Z, [prcmax1, prcmin3]);
% build condition based on median range percentiles for the fixed variables
%fixed_p = {}
fixed_cond = true(size(fixed_vars));
for i = 1:numel(fixed_vars)
    %fixed_p{i} = prctile(fixed_vars{i}, [prcmin2, prcmax2]);
    fixed_p = prctile(fixed_vars{i}(:), [prcmin2, prcmax2]);
    fixed_cond = fixed_cond & (fixed_vars{i}(:) > fixed_p(1)  & fixed_vars{i}(:) < fixed_p(2));
end

% find locations of low values of study_var where fixed_vars are close to median
i_small = find(study_var.Z(:) < study_p(1) & fixed_cond);
[x_low, y_low] = ind2coord(study_var, i_small);
% find locations of high values of study_var where fixed_vars are close to median
i_high = find(study_var.Z(:) > study_p(2) & fixed_cond);
[x_high, y_high] = ind2coord(study_var, i_high);
end