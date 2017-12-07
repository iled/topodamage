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
canopy = GRIDobj('resources/canopyclipped.tif');
impervious = GRIDobj('resources/impervious_final.tif');
soil = GRIDobj('resources/soils_13.tif');

%% load non cropped filtered versions
% topographic metrics
drainage_area_nocrop = GRIDobj('resources/NoCrop/drainage_area_mdf_filtered_noCrop.tif');
drainage_density_nocrop = GRIDobj('resources/NoCrop/drainage_density_filtered_noCrop.tif');
slope_nocrop = GRIDobj('resources/NoCrop/slope_filtered_noCrop.tif');
wetness_index_nocrop = GRIDobj('resources/NoCrop/wetness_index_filtered_noCrop.tif');
% non-topographic metrics
housing_age_nocrop = GRIDobj('resources/NoCrop/housing_age_filtered_noCrop.tif');
canopy_nocrop = GRIDobj('resources/NoCrop/canopy_filtered_noCrop.tif');
impervious_nocrop = GRIDobj('resources/NoCrop/impervious_filtered_noCrop.tif');
soil_nocrop = GRIDobj('resources/NoCrop/soil_filtered_noCrop.tif');

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

%% interp
dd = interp(drainage_density_nocrop, damage.x, damage.y);

%% Dummy interp just to fetch the projected coordinates
[~, damage.x, damage.y] = geotiffinterp('resources/DEM30_gauss_filled.tif', damage.Lat, damage.Long);

%% PLOT: sample sites over DEM
imagesc(DEM)
hold on
gscatter(damage.x, damage.y, damage.SiteType, 'mgrb', '+', 8, 'on')
legend('Location', 'northwestoutside')
hold off

%% PLOT: sample sites over DEM + site selection
sites = load('coordinates.mat');
figure
imagesc(drainage_density_nocrop)
cz = colorbar;
cz.Label.String = 'Elevation [masl]';
colors = {'g', 'r', 'm'};
hold on
gscatter(damage.x, damage.y, damage.SiteType, 'mgrb', '+', 8, 'off')
for i = 1:3
    plot(sites.coordinates{i}(:, 1), sites.coordinates{i}(:, 2), 'p', ...
        'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{i})
end
legend('DD_L', 'DD_H', 'DA_L', 'DA_H', 'S_H', 'S_L', ...
    'Site DA', 'Site Slope', 'Site DD', 'Location', 'northwestoutside')
title('Site locations over the DEM')
shg
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
gscatter(damage.S, damage.LocalSlopeDeg, damage.SiteType, ...
    [], [], [], 'on', 'Regional slope', ...
    'Local slope')

%% Group stats: Local slope and damage magnitude per lows and highs
grpstats(damage, 'HighLow', {'min', 'max', 'mean', @mode}, 'DataVars', ...
    {'LocalSlopeDeg', 'DamageMagnitude'})

%% Regional slope vs damage direction
rs_dir = damage.RegionalSlopeDirection_Az_;
rs_deg = damage.LocalSlopeDeg;
dmg_dir = damage.DamageMetricsOrientation_Az_;
tilt = damage.DamageMetricsTilt_degreeFromVertical_;
% select the points where the regional slope and the damage direction are more
% than 20 degrees apart
pick = find(abs(rs_dir - dmg_dir) > 20);
% also save the remaining sites
pick2 = abs(rs_dir - dmg_dir) <= 20;
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
subplot(1, 2, 1)
gscatter(damage.categoricAge1_young_2_mid_3_old, damage.ageHousing, ..., 
    damage.SiteType, [], [], [], 'on', 'Categorical age (1: Young ? 3: Old)', ...
    'Absolute age (mode per block)')
title('All types of structure')
subplot(1, 2, 2)
gscatter(damage.categoricAge1_young_2_mid_3_old(walls), damage.ageHousing(walls), ..., 
    damage.SiteType(walls), [], [], [], 'on', 'Categorical age (1: Young ? 3: Old)', ...
    'Absolute age (mode per block)')
title('Only considering walls')
%% histogram of each age measure
subplot(1, 3, 1)
histogram(damage.categoricAge1_young_2_mid_3_old)
subplot(1, 3, 2)
histogram(damage.AbsoluteAge_yearBc_)
subplot(1, 3, 3)
histogram(damage.ageHousing)

%% only walls
subplot(1, 2, 1)
histogram(damage.categoricAge1_young_2_mid_3_old(walls), 3)
title('Categorical age (1: Young ? 3: Old)')
subplot(1, 2, 2)
histogram(damage.ageHousing(walls))
title('Absolute age (mode per block)')

%% plot matrix
f = figure('Visible', 'on', 'NumberTitle', 'off', ...
    'units','normalized','outerposition',[0 0 1 1]);
topo_metrics = {'DA', 'DD_nearestNeighbor', 'S', 'W_I_'};
damage_metrics = {'DamagedStructure', 'DamageMagnitude'};
gplotmatrix(table2array(damage(:, topo_metrics)), ...
    table2array(damage(:, damage_metrics)), damage.SiteType, ...
    [], 'o', 10, [], '', topo_metrics, damage_metrics)
saveas(f, 'plotmatrix_topo_vs_damage', 'png');

%% same but only for sidewalks and walls
f = figure('Visible', 'on', 'NumberTitle', 'off', ...
    'units','normalized','outerposition',[0 0 1 1]);
topo_metrics = {'DA', 'DD', 'DD_nearestNeighbor', 'S', 'W_I_'};
damage_metrics = {'DamagedStructure', 'DamageMagnitude'};
gplotmatrix(table2array(damage(sidewalls, topo_metrics)), ...
    table2array(damage(sidewalls, damage_metrics)), damage.SiteType(sidewalls), ...
    [], 'o', 10, [], '', topo_metrics, damage_metrics)
saveas(f, 'plotmatrix_topo_vs_damage_sidewalks+walls', 'png');

%% plot matrix tilt and length vs topo metrics -- all sites
f = figure('Visible', 'on', 'NumberTitle', 'off', ...
    'units','normalized','outerposition',[0 0 1 1]);
topo_metrics = {'DA', 'DD_nearestNeighbor', 'S', 'W_I_'};
damage_metrics = {'DamageMetricsTilt_degreeFromVertical_', 'DamageMetricsLength_m_'};
gplotmatrix(log(table2array(damage(:, topo_metrics))), ...
    log(table2array(damage(:, damage_metrics))), damage.SiteType, ...
    [], 'o', 10, [], '', topo_metrics, damage_metrics)
saveas(f, 'plotmatrix_LOG_topo_vs_tilt_length_all_sites', 'png');

%% plot matrix tilt and length vs topo metrics -- only slope site
f = figure('Visible', 'on', 'NumberTitle', 'off', ...
    'units','normalized','outerposition',[0 0 1 1]);
topo_metrics = {'S', 'LocalSlopeDeg'};
damage_metrics = {'DamageMetricsTilt_degreeFromVertical_', 'DamageMetricsLength_m_'};
slopes = cellfun(@(x) ~isempty(x), strfind(damage.SiteType, 'S_'));
% only those where directions agree up to 20 deg
slopes2 = slopes & pick2;
gplotmatrix(table2array(damage(slopes, topo_metrics)), ...
    table2array(damage(slopes, damage_metrics)), damage.HighLow(slopes), ...
    [], 'o', 10, [], '', topo_metrics, damage_metrics)
saveas(f, 'plotmatrix_topo_vs_tilt_length_slope_sites', 'png');

%% Agg site type: damage measurement vs magnitude
grpstats(damage, 'Person_s_Collecting', {'mean'}, 'DataVars', ...
    {'DamageMetricsArea_m_2_', 'DamageMetricsDisplacement_m_', ...
    'DamageMetricsLength_m_', 'DamageMetricsTilt_degreeFromVertical_', ...
    'DamageMetricsWidth_m_', 'DamageMagnitude'})

%%
metrics = {'DamageMetricsArea_m_2_', 'DamageMetricsDisplacement_m_', ...
    'DamageMetricsLength_m_', 'DamageMetricsTilt_degreeFromVertical_', ...
    'DamageMetricsWidth_m_'};
gplotmatrix(table2array(damage(:, {'DamageMagnitude'})), ...
    table2array(damage(:, metrics)), damage.Person_s_Collecting(:), ...
    [], 'o', 10, [], '', 'DamageMagnitude', metrics)

%%
metrics = {'DamageMetricsLength_m_', 'DamageMetricsTilt_degreeFromVertical_'};
gplotmatrix(table2array(damage(:, {'DamageMagnitude'})), ...
    log(table2array(damage(:, metrics))), damage.Person_s_Collecting(:), ...
    [], 'o', 10, [], '', 'DamageMagnitude', metrics)

%% histograms per group
gplotmatrix(table2array(damage(:, {'DamageMagnitude'})), ...
    [], damage.Person_s_Collecting(:), ...
    [], 'o', 10, 'on', 'grpbars', 'DamageMagnitude')

%%
figure
grps = unique(damage.Person_s_Collecting);
for i = 1:4
    subplot(2, 2, i)
    histogram(damage.DamageMagnitude(strcmp(damage.Person_s_Collecting, grps{i})), 5)
    title(grps{i})
end

%% histograms of tilt and length
figure
grps = unique(damage.Person_s_Collecting);
for i = 1:4
    subplot(3, 4, i)
    histogram(damage.DamageMetricsTilt_degreeFromVertical_(strcmp(damage.Person_s_Collecting, grps{i})), 5)
    title([grps{i} ': Tilt'])
    subplot(3, 4, i + 4)
    histogram(damage.DamageMetricsLength_m_(strcmp(damage.Person_s_Collecting, grps{i})), 5)
    title([grps{i} ': Length'])
    subplot(3, 4, i + 8)
    histogram(damage.DamageMagnitude(strcmp(damage.Person_s_Collecting, grps{i})), 5)
    title([grps{i} ': Damage magnitude'])
end