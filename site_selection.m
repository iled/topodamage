%% DATA: load data
% topographic metrics
DEM = GRIDobj('resources/DEM30_gauss_filled.tif');
drainage_area = GRIDobj('resources/drainage_area_mdf.tif');
drainage_density = GRIDobj('resources/drainage_density_fixed.tif');
slope = GRIDobj('resources/slope_gauss.tif');
wetness_index = GRIDobj('resources/wetness_index.tif');
% non-topographic metrics
housing_age = GRIDobj('resources/agehouses_tmp3.tif');
canopy  = GRIDobj('resources/canopyclipped.tif');
impervious = GRIDobj('resources/impervious_final.tif');
soil = GRIDobj('resources/soils_13.tif');

% fix impervious map
impervious.refmat = DEM.refmat;
impervious.cellsize = DEM.cellsize;
% fix age of housing
housing_age.Z(housing_age.Z < 1755) = NaN;
housing_age.refmat = DEM.refmat;
housing_age.cellsize = DEM.cellsize;

%% RUN: compute averages over 3 standard city blocks (330 m x 330 m)
% excluded, have onw filters below: drainage density, housing age, soil
DEM_avg = DEM;
drainage_area_avg = drainage_area;
slope_avg = slope;
wetness_index_avg = wetness_index;
canopy_avg = canopy;
impervious_avg = impervious;

n_px = 11;
h = ones(n_px, n_px) / n_px ^ 2;

grids = {DEM, drainage_area, slope, wetness_index, ...
    canopy, impervious};
grids_avg = {'DEM_avg', 'drainage_area_avg', 'slope_avg', 'wetness_index_avg', ...
    'canopy_avg', 'impervious_avg'};

for i = 1:numel(grids)
    tempval = eval(grids_avg{i});
    tempval.Z = imfilter(grids{i}.Z, h);
    assignin('base', grids_avg{i}, tempval);
end

%% RUN: specific filter for drainage density: max filter with edge n_px
n_px = 5;
drainage_density_avg = drainage_density;
drainage_density_avg.Z = fillmissing(drainage_density_avg.Z, 'constant', -1);
drainage_density_avg.Z = ordfilt2(drainage_density_avg.Z, n_px ^ 2, ones(n_px, n_px));
drainage_density_avg.Z(drainage_density_avg.Z == -1) = NaN;

%% FIG: plot original drainage density vs filtered
subplot(1, 2, 1)
imagesc(drainage_density), colorbar
subplot(1, 2, 2)
imagesc(drainage_density_avg), colorbar

%% RUN: specific filter for age of housing: min filter with edge n_px
n_px = 5;
housing_age_avg = housing_age;
housing_age_avg.Z = fillmissing(housing_age_avg.Z, 'constant', 9999);
housing_age_avg.Z = ordfilt2(housing_age_avg.Z, 1, ones(n_px, n_px));
housing_age_avg.Z(housing_age_avg.Z == 9999) = NaN;
housing_age_avg.Z(housing_age_avg.Z == 0) = NaN;

%% DEBUG: plot original housing age vs filtered
subplot(1, 2, 1)
imagesc(housing_age), colorbar
subplot(1, 2, 2)
imagesc(housing_age_avg), colorbar

%% RUN: specific filter for soil: mode filter with edge n_px
% SLOW: this is slow, may take a few minutes, depending on n_px
n_px = 5;  % didn't work well with 11
soil_avg = soil;
soil_avg.Z = mode_filter(soil_avg.Z, n_px);

%% DEBUG: plot original soil vs filtered
figure
colormap(jet(numel(unique(soil_avg.Z))))
subplot(1, 2, 1)
imagesc(soil), colorbar
subplot(1, 2, 2)
imagesc(soil_avg), colorbar

%% RUN: update grids cell arrays
grids{1, end + 1} = drainage_density;
grids_avg{1, end + 1} = 'drainage_density_avg';
grids{1, end + 1} = housing_age;
grids_avg{1, end + 1} = 'housing_age_avg';
grids{1, end + 1} = soil;
grids_avg{1, end + 1} = 'soil_avg';

%% DEBUG: save plot
print('age_houses_minfilter5x5', '-djpeg', '-r300')

%% RUN: crop averaged grids to fix borders after filtering
for i = 1:numel(grids)
    tempval = eval(grids_avg{i});
    tempval.Z = imcrop(tempval.Z, [9, 8, 572 - 17, 574 - 16]);
    tempval.size = size(tempval.Z);
    assignin('base', grids_avg{i}, tempval);
end

%% DEBUG: cell size and refmat
clc
for i = 1:numel(grids)
    tempval = eval(grids_avg{i});
    fprintf('grid: %s \t cellsize: %.2f \n refmat: \n', grids_avg{i}, ...
        tempval.cellsize)
    tempval.refmat(end, :)
    fprintf('size: %d, %d\n\n', size(tempval.Z, 1), size(tempval.Z, 2))
end

%% FIG: plot the averaged grids
for i = 1:numel(grids)
    subplot(3, 3, i)
    imagesc(eval(grids_avg{i}));
    title(grids_avg{i}, 'Interpreter', 'none');
end

%% DEBUG: plot age of housing distribution
histogram(housing_age_avg.Z)

%% RUN: analyze each one of the topographic metrics, fix non-topographic
%% v1: using the percentiles criteria, only fixing nontopo, excluding soil
test_v = 'v1';
topo = {drainage_area_avg, drainage_density_avg, slope_avg, wetness_index_avg};
nontopo = {housing_age_avg, canopy_avg, impervious_avg};
locations = cell(1, numel(topo));
log_metrics = {wetness_index_avg.name};
for i = 1:numel(topo)
    [xl, yl, xh, yh] = perc_analysis(topo{i}, nontopo, log_metrics);
    locations{i} = {xl, yl, xh, yh};
end

%% v2: using the similar ranges criteria, only fixing nontopo, excluding soil
test_v = 'v2';
topo = {drainage_area_avg, drainage_density_avg, slope_avg, wetness_index_avg};
nontopo = {housing_age_avg, canopy_avg, impervious_avg};
locations = cell(1, numel(topo));
log_metrics = {wetness_index_avg.name};
for i = 1:numel(topo)
    [xl, yl, xh, yh] = range_analysis(topo{i}, nontopo, log_metrics);
    locations{i} = {xl, yl, xh, yh};
end

%% v3: using the similar ranges criteria, fixing topo
test_v = 'v3';
topo = {drainage_area_avg, drainage_density_avg, slope_avg, wetness_index_avg};
nontopo = {housing_age_avg, canopy_avg, impervious_avg, soil_avg};
locations = cell(1, numel(topo));
log_metrics = {wetness_index_avg.name};
%log_metrics = {drainage_area_avg.name, drainage_density_avg.name, wetness_index_avg.name};
for i = 1:numel(topo)
    [xl, yl, xh, yh] = range_analysis(topo{i}, topo(1:end ~= i), log_metrics);
    locations{i} = {xl, yl, xh, yh};
end

%% v4: using the similar ranges criteria, fixing topo and all nontopo, except soil
test_v = 'v4';
topo = {drainage_area_avg, drainage_density_avg, slope_avg, wetness_index_avg};
nontopo = {housing_age_avg, canopy_avg, impervious_avg};
locations = cell(1, numel(topo));
log_metrics = {wetness_index_avg.name};
%log_metrics = {drainage_area_avg.name, drainage_density_avg.name, wetness_index_avg.name};
for i = 1:numel(topo)
    [xl, yl, xh, yh] = range_analysis(topo{i}, [topo(1:end ~= i) nontopo], log_metrics);
    locations{i} = {xl, yl, xh, yh};
end

%% v5: using the percentiles criteria, fixing topo and all nontopo, except soil
test_v = 'v5';
topo = {drainage_area_avg, drainage_density_avg, slope_avg, wetness_index_avg};
nontopo = {housing_age_avg, canopy_avg, impervious_avg};
locations = cell(1, numel(topo));
log_metrics = {wetness_index_avg.name};
%log_metrics = {drainage_area_avg.name, drainage_density_avg.name, wetness_index_avg.name};
for i = 1:numel(topo)
    [xl, yl, xh, yh] = perc_analysis(topo{i}, [topo(1:end ~= i) nontopo], log_metrics);
    locations{i} = {xl, yl, xh, yh};
end

%% DEBUG: log drainage area and density
% this is a test to evaluate the effect of considering drainage density
% and drainage area as log-like distributed
% add this to the v-test: log_metrics = {log_drainage_area_avg.name};
log_drainage_area_avg = drainage_area_avg;
log_drainage_area_avg.Z = log(log_drainage_area_avg.Z);
log_drainage_density_avg = drainage_density_avg;
log_drainage_density_avg.Z = log(log_drainage_density_avg.Z);

%% v6: using the percentiles criteria, fixing topo and nontopo, exclude soil
% exclude soil and drainage density
test_v = 'v6';
topo = {drainage_area_avg, drainage_density_avg, slope_avg, wetness_index_avg};
nontopo = {housing_age_avg, canopy_avg, impervious_avg};
locations = cell(1, numel(topo));
log_metrics = {wetness_index_avg.name};
%log_metrics = {drainage_area_avg.name, drainage_density_avg.name, wetness_index_avg.name};
for i = 1:numel(topo)
    [xl, yl, xh, yh] = perc_analysis(topo{i}, ...
        [topo(~ismember(1:end, [i, 2])) nontopo], log_metrics);
    locations{i} = {xl, yl, xh, yh};
end

%% v7: using the similar ranges criteria, fixing topo and nontopo,
% exclude soil and drainage density
test_v = 'v7';
topo = {drainage_area_avg, drainage_density_avg, slope_avg, wetness_index_avg};
nontopo = {housing_age_avg, canopy_avg, impervious_avg}; %, soil_avg
locations = cell(1, numel(topo));
log_metrics = {wetness_index_avg.name};
%log_metrics = {drainage_area_avg.name, drainage_density_avg.name, wetness_index_avg.name};
for i = 1:numel(topo)
    [xl, yl, xh, yh] = range_analysis(topo{i}, ...
        [topo(~ismember(1:end, [i, 2])) nontopo], log_metrics);
    locations{i} = {xl, yl, xh, yh};
end

%% v8: using the similar ranges criteria, fixing topo and nontopo,
% exclude soil, canopy and drainage density
test_v = 'v8';
topo = {drainage_area_avg, drainage_density_avg, slope_avg, wetness_index_avg};
nontopo = {housing_age_avg, impervious_avg}; %, soil_avg canopy_avg,
locations = cell(1, numel(topo));
log_metrics = {wetness_index_avg.name};
%log_metrics = {drainage_area_avg.name, drainage_density_avg.name, wetness_index_avg.name};
for i = 1:numel(topo)
    [xl, yl, xh, yh] = range_analysis(topo{i}, ...
        [topo(~ismember(1:end, [i, 2])) nontopo], log_metrics);
    locations{i} = {xl, yl, xh, yh};
end

%% FIG: plot each of the topographic metrics, with fixed non-topographic
%close all
figure('NumberTitle', 'off', 'Name', test_v)
topo_title = {'Drainage area', 'Drainage density', 'Slope', 'Wetness index'};

for i = 1:numel(topo)
    subplot(2, 2, i)
    %imageschs(DEM)
    imagesc(topo{i});
    hold on
    coord = locations{i};
    plot(coord{1}, coord{2}, 'og', coord{3}, coord{4},'sr')
    legend('low', 'high')
    %imagesc(topo{i})
    title(topo_title{i});
    colorbar
    hold off
end

shg
%% DEBUG: save plot
print(['site_selection_' test_v], '-djpeg', '-r300')