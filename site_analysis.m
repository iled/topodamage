%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topographic Analysis - HW7
% Sam Mark, Arielle Woods, Julio Caineta
% Analysis of infrastructure damage
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
drainage_density_f = GRIDobj('resources/drainage_density_filtered.tif');
slope = GRIDobj('resources/slope_gauss.tif');
wetness_index = GRIDobj('resources/wetness_index.tif');
% non-topographic metrics
housing_age = GRIDobj('resources/agehouses_tmp3.tif');
canopy  = GRIDobj('resources/canopyclipped.tif');
impervious = GRIDobj('resources/impervious_final.tif');
soil = GRIDobj('resources/soils_13.tif');

%%
% fix impervious map
impervious.refmat = DEM.refmat;
impervious.cellsize = DEM.cellsize;
% fix age of housing
housing_age.Z(housing_age.Z < 1755) = NaN;
housing_age.refmat = DEM.refmat;
housing_age.cellsize = DEM.cellsize;

% path to filtered files
% (DEM_avg, 'resources/DEM_filtered.tif')
% (drainage_area_avg, 'resources/drainage_area_mdf_filtered.tif')
% (slope_avg, 'resources/slope_filtered.tif')
% (wetness_index_avg, 'resources/wetness_index_filtered.tif')
% (drainage_density_avg, 'resources/drainage_density_filtered.tif')
% (housing_age_avg, 'resources/housing_age_filtered.tif')
% (impervious_avg, 'resources/impervious_filtered.tif')
% (soil_avg, 'resources/soil_filtered.tif')
% (canopy_avg, 'resources/canopy_filtered.tif'

%% load damage data
damage_fn = 'resources/TopoDataAnalysis_all.csv';
opts = detectImportOptions(damage_fn);
% fix numeric columns loaded as char
% [opts.VariableNames; opts.VariableTypes]
opts.VariableTypes(strcmp(opts.VariableNames, 'DD')) = {'double'};
opts.VariableTypes(strcmp(opts.VariableNames, 'DD_nearestNeighbor')) = {'double'};
% read table
damage = readtable(damage_fn, opts);
% rename columns
damage.Properties.VariableNames{'DamagedStructure_1_sidewalk_2_road_3_fence_4_house_wall_5'} = 'DamagedStructure';
damage.Properties.VariableNames{'magnitudeOfDamage______Categorical_1_low_5Hig_'} = 'DamageMagnitude';
damage.Properties.VariableNames{'MaterialType_1_concrete_2_brick_3_asphalt_4_stone_5_wood_'} = 'MaterialType';
damage.Properties.VariableNames{'RegionalSlopeMagnitude_degreesFromHorizontal_'} = 'LocalSlopeDeg';

%% Add variable to distinguish sites of high and low topo metric
highs = cellfun(@(x) ~isempty(x), strfind(damage.SiteType, '_H'));
lows = cellfun(@(x) ~isempty(x), strfind(damage.SiteType, '_L'));
highlow = {};
highlow(highs, 1) = repmat({'High'}, [sum(highs) 1]);
highlow(lows, 1) = repmat({'Low'}, [sum(lows) 1]);
damage.HighLow = highlow;

%% get topographic metric values from maps
damage.DrainageDensity = geotiffinterp('resources/drainage_density_fixed.tif', damage.Lat, damage.Long);
damage.DrainageArea = geotiffinterp('resources/drainage_area_mdf.tif', damage.Lat, damage.Long);
[damage.Slope, x, y] = geotiffinterp('resources/slope_gauss.tif', damage.Lat, damage.Long);
damage.WetnessIndex = geotiffinterp('resources/wetness_index.tif', damage.Lat, damage.Long);

%% Dummy interp just to fetch the projected coordinates
[~, damage.x, damage.y] = geotiffinterp('resources/slope_gauss.tif', damage.Lat, damage.Long);

%% PLOT: sample sites over DEM
imagesc(DEM)
hold on
%plot(x, y, 'r*')
gscatter(damage.x, damage.y, damage.SiteType, 'mgrb', '+', 8, 'off')
legend('DD High', 'DD Low', 'Slope High', 'Slope Low', ...
        'Location', 'northwestoutside')
hold off

%% Frequency table: Site type
tabulate(damage.SiteType)
%% Frequency table: Magnitude of damage
tabulate(damage.DamageMagnitude)
%% Frequency table: Damaged infrastructure
tabulate(damage.DamagedStructure)
%% Frequency table: Material
tabulate(damage.MaterialType)

%% Contigency table: Site type vs Magnitude of damage
[tbl, chi2, p, labels] = crosstab(damage.SiteType, damage.DamageMagnitude)
%% Contigency table: Site type vs Damaged infrastructure
[tbl, chi2, p, labels] = crosstab(damage.SiteType, damage.DamagedStructure)

%% Group stats: Damaged structure and magnitude per site type
% grpstats requires the Statistics and Machine Learning Toolbox
grpstats(damage, 'SiteType', {'min', 'max', 'mean', @mode}, 'DataVars', ... 
    {'DamagedStructure', 'DamageMagnitude'})

%% Group plot: Damaged structure and magnitude per site type
% not very clear
gscatter(damage.DamageMagnitude, damage.DamagedStructure, damage.SiteType)

%% Group stats: Magnitude per site type -- only sidewalks
sidewalks = damage.DamagedStructure == 1;
grpstats(damage(sidewalks, :), 'SiteType', {'min', 'max', 'mean', @mode}, 'DataVars', ... 
    {'DamageMagnitude'})

%% Group stats: Magnitude per site type -- only walls
walls = damage.DamagedStructure == 4;
grpstats(damage(walls, :), 'SiteType', {'min', 'max', 'mean', @mode}, 'DataVars', ... 
    {'DamageMagnitude'})

%% Group stats: Magnitude per site type -- only sidewalks and walls
sidewalls = sidewalks | walls;
grpstats(damage(sidewalls, :), 'SiteType', {'min', 'max', 'mean', @mode}, 'DataVars', ... 
    {'DamageMagnitude'})

%% Group stats: Magnitude per damaged structure
grpstats(damage, 'DamagedStructure', {'min', 'max', 'mean', @mode}, 'DataVars', ... 
    {'DamageMagnitude'})

%% Group stats: Topographic metrics + 2 non-topo per site type
grpstats(damage, 'SiteType', {'min', 'max', 'mean', @mode}, 'DataVars', ...
    {'DA', 'S', 'W_I_', 'DD', 'DD_nearestNeighbor', 'impervious', 'ageHousing'})

%% Group stats: Computed slope vs local slope per site type
grpstats(damage, 'SiteType', {'min', 'max', 'mean', @mode}, 'DataVars', ...
    {'S', 'LocalSlopeDeg', 'DamageMagnitude'})

%% Group plot: Computed slope vs local slope per site type
gscatter(damage.S, damage.LocalSlopeDeg, damage.SiteType)

%% Group stats: Local slope and damage magnitude per lows and highs
grpstats(damage, 'HighLow', {'min', 'max', 'mean', @mode}, 'DataVars', ...
    {'LocalSlopeDeg', 'DamageMagnitude'})

%% Regional slope vs damage direction
rs_dir = damage.RegionalSlopeDirection_Az_;
rs_deg = damage.RegionalSlopeMagnitude_degreesFromHorizontal_;
dmg_dir = damage.DamageMetricsOrientation_Az_;
tilt = damage.DamageMetricsTilt_degreeFromVertical_;

% select the points where the regional slope and the damage direction are more
% than 20 degrees apart
pick = find(abs(rs_dir - dmg_dir) > 20);
% same but select only points where some tilt was measured
tilted = find(abs(rs_dir - dmg_dir) > 20 & tilt > 0);
t = table(rs_dir(tilted), dmg_dir(tilted), rs_deg(tilted), tilt(tilted), ...
    'RowNames', damage.SiteID(tilted), ...
    'VariableNames', {'RS_dir', 'Tilt_dir', 'RS_deg', 'Tilt_deg'})

%% compare different age measures
gscatter(damage.ageHousing, damage.AbsoluteAge_yearBc_, damage.SiteType)
%%
gscatter(damage.categoricAge1_young_2_mid_3_old, damage.AbsoluteAge_yearBc_, damage.SiteType)
%%
gscatter(damage.categoricAge1_young_2_mid_3_old, damage.ageHousing, damage.SiteType)
%% histogram of each age measure
subplot(1, 3, 1)
histogram(damage.categoricAge1_young_2_mid_3_old)
subplot(1, 3, 2)
histogram(damage.AbsoluteAge_yearBc_)
subplot(1, 3, 3)
histogram(damage.ageHousing)

%% plot matrix
f = figure('Visible', 'on', 'NumberTitle', 'off', ...
    'units','normalized','outerposition',[0 0 1 1]);
topo_metrics = {'DA', 'DD', 'DD_nearestNeighbor', 'S', 'W_I_'};
damage_metrics = {'DamagedStructure', 'DamageMagnitude'};
gplotmatrix(table2array(damage(:, topo_metrics)), ...
    table2array(damage(:, damage_metrics)), damage.SiteType, ...
    [], 'o', 10, [], '', topo_metrics, damage_metrics)
saveas(f, 'plotmatrix_topo_vs_damage', 'png');