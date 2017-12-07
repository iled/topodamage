%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topographic Analysis - HW5
% Sam Mark, Arielle Woods, Julio Caineta
% Topographic maps generation
% 
% Sections start with meaningful keywords:
% RUN: section is required to be run
% FIG: section only used to make some plot
% DEBUG: section used for debugging purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RUN: load data
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
cdu = GRIDobj('resources/cdu3final.tif');

% fix impervious map
impervious.refmat = DEM.refmat;
impervious.cellsize = DEM.cellsize;
% fix age of housing
housing_age.Z(housing_age.Z < 1755) = NaN;
housing_age.refmat = DEM.refmat;
housing_age.cellsize = DEM.cellsize;

% titles of the four topographic metrics (used to plot)
topo_title = {'Drainage area', 'Drainage density', 'Slope', 'Wetness index'};

%% RUN: compute averages over 3 standard city blocks (330 m x 330 m)
% excluded here, as they have onw filters below: drainage density, housing age, soil
DEM_avg = DEM;
drainage_area_avg = drainage_area;
slope_avg = slope;
wetness_index_avg = wetness_index;
canopy_avg = canopy;
impervious_avg = impervious;

% ***** attempt to 
CDU.Z=changem(CDU.Z,NaN,0);
%GRIDobj2geotiff(CDU,'resources/CDU_NaN.tif')
%CDU_avg = CDU;

n_px = 11;
h = ones(n_px, n_px) / n_px ^ 2;

grids = {DEM, drainage_area, slope, wetness_index, ...
    canopy, impervious, CDU};
grids_avg = {'DEM_avg', 'drainage_area_avg', 'slope_avg', 'wetness_index_avg', ...
    'canopy_avg', 'impervious_avg', 'CDU_avg'};

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
title('original drainage density')
subplot(1, 2, 2)
imagesc(drainage_density_avg), colorbar
title('filtered drainage density with max filter')

%% RUN: specific filter for age of housing: min filter with edge n_px
n_px = 5; % didn't work well with n_px = 11
housing_age_avg = housing_age;
housing_age_avg.Z = fillmissing(housing_age_avg.Z, 'constant', 9999);
housing_age_avg.Z = ordfilt2(housing_age_avg.Z, 1, ones(n_px, n_px));
housing_age_avg.Z(housing_age_avg.Z == 9999) = NaN;
housing_age_avg.Z(housing_age_avg.Z == 0) = NaN;

%% RUN: specific filter for age of housing: mode filter with edge n_px
n_px = 5;
housing_age_mode = housing_age;
housing_age_mode.Z(housing_age_mode.Z == 0) = NaN;
housing_age_mode.Z = colfilt(housing_age_mode.Z, [n_px, n_px], 'sliding', @mode);
housing_age_mode.Z(housing_age_mode.Z == 0) = NaN;

%% DEBUG: plot original housing age vs filtered min vs filtered mode
subplot(1, 3, 1)
imagesc(housing_age), colorbar
title('original housing age (year built)')
subplot(1, 3, 2)
imagesc(housing_age_avg), colorbar
title('filtered housing age with min filter')
subplot(1, 3, 3)
imagesc(housing_age_mode), colorbar
title('filtered housing age with mode filter')

%% DEBUG: plot histograms for original housing age vs filtered min vs filtered mode
subplot(1, 3, 1)
histogram(housing_age.Z)
title('original housing age (year built)')
subplot(1, 3, 2)
histogram(housing_age_avg.Z)
title('filtered housing age with min filter')
subplot(1, 3, 3)
histogram(housing_age_mode.Z)
title('filtered housing age with mode filter')

%% RUN: specific filter for CDU: min filter with edge n_px
n_px = 5; % didn't work well with n_px = 11
cdu_avg = cdu;
cdu_avg.Z(cdu_avg.Z == 0) = NaN;
cdu_avg.Z = fillmissing(cdu_avg.Z, 'constant', 9999);
cdu_avg.Z = ordfilt2(cdu_avg.Z, 1, ones(n_px, n_px));
cdu_avg.Z(cdu_avg.Z == 9999) = NaN;
cdu_avg.Z(cdu_avg.Z == 0) = NaN;

%% RUN: specific filter for CDU: mode filter with edge n_px
n_px = 5;
cdu.Z(cdu.Z == 0) = NaN;
cdu_mode = cdu;
cdu_mode.Z(cdu_mode.Z == 0) = NaN;
cdu_mode.Z = colfilt(cdu_mode.Z, [n_px, n_px], 'sliding', @mode);
cdu_mode.Z(cdu_mode.Z == 0) = NaN;

%% DEBUG: plot original CDU vs filtered min vs filtered mode
subplot(1, 3, 1)
imagesc(cdu), colorbar
title('original cdu')
subplot(1, 3, 2)
imagesc(cdu_avg), colorbar
title('filtered cdu with min filter')
subplot(1, 3, 3)
imagesc(cdu_mode), colorbar
title('filtered cdu with mode filter')

%% DEBUG: plot histograms for CDU vs filtered min vs filtered mode
subplot(1, 3, 1)
histogram(cdu.Z, 8)
title('original cdu')
subplot(1, 3, 2)
histogram(cdu_avg.Z, 8)
title('filtered cdu with min filter')
subplot(1, 3, 3)
histogram(cdu_mode.Z, 8)
title('filtered cdu with mode filter')

%% RUN: specific filter for soil: mode filter with edge n_px
n_px = 11;
soil_avg = soil;
soil_avg.Z = colfilt(soil.Z, [n_px, n_px], 'sliding', @mode);

%% DEBUG: plot original soil vs filtered
figure
colormap(jet(numel(unique(soil_avg.Z))))
subplot(1, 2, 1)
imagesc(soil), colorbar
title('original soils type')
subplot(1, 2, 2)
imagesc(soil_avg), colorbar
title('filtered soils type with mode filter')

%% RUN: update grids cell arrays
grids{1, end + 1} = drainage_density;
grids_avg{1, end + 1} = 'drainage_density_avg';
grids{1, end + 1} = housing_age;
grids_avg{1, end + 1} = 'housing_age_avg';
grids{1, end + 1} = soil;
grids_avg{1, end + 1} = 'soil_avg';

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

%% FIG: plot the filtered grids
for i = 1:numel(grids)
    subplot(3, 3, i)
    imagesc(eval(grids_avg{i}));
    title(grids_avg{i}, 'Interpreter', 'none');
end

%% DEBUG: save the filtered grids
GRIDobj2geotiff(DEM_avg, 'resources/DEM_filtered.tif')
GRIDobj2geotiff(drainage_area_avg, 'resources/drainage_area_mdf_filtered.tif')
GRIDobj2geotiff(slope_avg, 'resources/slope_filtered.tif')
GRIDobj2geotiff(wetness_index_avg, 'resources/wetness_index_filtered.tif')
GRIDobj2geotiff(drainage_density_avg, 'resources/drainage_density_filtered.tif')
GRIDobj2geotiff(housing_age_avg, 'resources/housing_age_filtered.tif')
GRIDobj2geotiff(impervious_avg, 'resources/impervious_filtered.tif')
GRIDobj2geotiff(soil_avg, 'resources/soil_filtered.tif')
GRIDobj2geotiff(canopy_avg, 'resources/canopy_filtered.tif')

%% DEBUG: plot age of housing distribution
histogram(housing_age_avg.Z)

%% DEBUG: log drainage area and density
% this is a test to evaluate the effect of considering drainage density
% and drainage area as log-like distributed
% add this to the v-test: log_metrics = {log_drainage_area_avg.name};
log_drainage_area_avg = drainage_area_avg;
log_drainage_area_avg.Z = log(log_drainage_area_avg.Z);
log_drainage_density_avg = drainage_density_avg;
log_drainage_density_avg.Z = log(log_drainage_density_avg.Z);

%%%%%%%%%%%%%%%%%%%%%%%%
% SITE SELECTION TESTS %
%%%%%%%%%%%%%%%%%%%%%%%%
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
    [xl, yl, xh, yh, ranges_l, ranges_h] = range_analysis(topo{i}, ...
        [topo(~ismember(1:end, [i, 2])) nontopo], log_metrics);
    locations{i} = {xl, yl, xh, yh, ranges_l, ranges_h};
end

%% FIG: plot the results of the last site selection test (vN)
%close all
test_name = test_v;
save_figure = 'on';
plot_site_selection(topo, locations, topo_title, test_name, save_figure)
shg

%% FIG: plot the result of the last site selection test (vN)
% these plots are against the soil raster
%close all
test_name = [test_v '_soil'];
save_figure = 'on';
plot_site_selection(soil_avg, locations, topo_title, test_name, save_figure)
shg

%% FIG: plot the result of the last site selection test (vN)
% these plots are against the canopy raster
%close all
test_name = [test_v '_canopy'];
save_figure = 'on';
plot_site_selection(canopy_avg, locations, topo_title, test_name, save_figure)
shg

%% DEBUG: save plot
print(['site_selection_' test_v], '-djpeg', '-r300')
test_v = test_bckp;

%% RUN: save selected points, used in the main file
save('site_selection', 'locations')

%% RUN: save selection of 6 sites, 2 per topographic metric
%% row 1 is high, row 2 is low
% the coordinates were taken from an interactive plot
drainage_area_coords = [585181.1115, 4477942.9546; 582121.5172, 4474555.5467];
slope_coords = [589278.7824, 4484034.8252; 589907.092, 4483843.6007];
drainage_density_coords = [583816.9, 4480538.1; 583350.8, 4480483.5];
coordinates = {drainage_area_coords, slope_coords, drainage_density_coords};
%save('coordinates', 'drainage_area_coords', 'slope_coords', 'drainage_density_coords');
save('coordinates', 'coordinates');
