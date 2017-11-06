%% load data
% topographic metrics
DEM = GRIDobj('resources/DEM30_gauss_filled.tif');
drainage_area = GRIDobj('resources/drainage_area_mdf.tif');
drainage_density = GRIDobj('resources/drainage_density.tif');
slope = GRIDobj('resources/slope_gauss.tif');
wetness_index = GRIDobj('resources/wetness_index.tif');
% non-topographic metrics
housing_age = GRIDobj('resources/agehouses_tmp3.tif');
canopy  = GRIDobj('resources/canopyclipped.tif');
impervious = GRIDobj('resources/impervious_final.tif');
soil = GRIDobj('resources/soils_final.tif');

% fix impervious map
impervious.refmat = DEM.refmat;
impervious.cellsize = DEM.cellsize;
% fix age of housing
housing_age.Z(housing_age.Z < 1755) = NaN;

%% compute averages over 3 standard city blocks (330 m x 330 m)
DEM_avg = DEM;
drainage_area_avg = drainage_area;
drainage_density_avg = drainage_density;
slope_avg = slope;
wetness_index_avg = wetness_index;
housing_age_avg = housing_age;
canopy_avg = canopy;
impervious_avg = impervious;
soil_avg = soil;

n_px = 11;
h = ones(n_px, n_px) / n_px ^ 2;

grids = {DEM, drainage_area, drainage_density, slope, wetness_index, ...
    housing_age, canopy, impervious, soil};
grids_avg = {'DEM_avg', 'drainage_area_avg', 'drainage_density_avg', 'slope_avg', 'wetness_index_avg', ...
    'housing_age_avg', 'canopy_avg', 'impervious_avg', 'soil_avg'};

for i = 1:numel(grids)
    tempval = eval(grids_avg{i});
    tempval.Z = imfilter(grids{i}.Z, h);
    assignin('base', grids_avg{i}, tempval);
end

%% plot the averaged grids
for i = 1:numel(grids)
    subplot(3, 3, i)
    imagesc(eval(grids_avg{i}));
    title(grids_avg{i});
end

%% plot age of housing distribution
histogram(housing_age_avg.Z)

%% analyze each one of the topographic metrics, fix non-topographic
topo = {drainage_area_avg, drainage_density_avg, slope_avg, wetness_index_avg};
nontopo = {housing_age, canopy_avg, impervious_avg, soil_avg}; % TODO: fix and use housing average
locations = cell(1, numel(topo));
for i = 1:numel(topo)
    [xl, yl, xh, yh] = perc_analysis(topo{i}, nontopo);
    locations{i} = {xl, yl, xh, yh};
end
%% plot each of the topographic metrics, with fixed non-topographic
close all

for i = 1:numel(topo)
    subplot(2, 2, i)
    %imageschs(DEM)
    imagesc(topo{i});
    hold on
    coord = locations{i};
    plot(coord{1}, coord{2}, 'ok', coord{3}, coord{4},'sr')
    %imagesc(topo{i})
    %title();
    hold off
end
