%% load data
% topographic metrics
DEM = GRIDobj('resources/DEM30_gauss_filled.tif');
drainage_area = GRIDobj('resources/drainage_area_mdf.tif');
drainage_density = GRIDobj('resources/drainage_density.tif');
slope = GRIDobj('resources/slope_gauss.tif');
wetness_index = GRIDobj('resources/wetness_index.tif');
% non-topographic metrics
housing_age = GRIDobj('resources/age_houses_PointToRaster_Cli71.tif');
canopy  = GRIDobj('resources/canopyclipped.tif');
impervious = GRIDobj('resources/impervious_final.tif');
soil = GRIDobj('resources/soils_final.tif');

% fix impervious map
impervious.refmat = DEM.refmat;
impervious.cellsize = DEM.cellsize;

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
h = ones(n_px, n_px);

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
histogram(housing_age.Z(housing_age.Z > 0))

%% analyze slope, fix others

