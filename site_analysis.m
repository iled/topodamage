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

%% filter rows of topo damage table by site type

dalRows = damage(strcmp(damage.SiteType, {'DA_L'}), :);
dahRows = damage(strcmp(damage.SiteType, {'DA_H'}), :);
daRows = damage(strcmp(damage.SiteType, {'DA_H'}) | strcmp(damage.SiteType, {'DA_L'}), :);
ddlRows = damage(strcmp(damage.SiteType, {'DD_L'}), :);
ddhRows = damage(strcmp(damage.SiteType, {'DD_H'}), :);
ddRows = damage(strcmp(damage.SiteType, {'DD_H'}) | strcmp(damage.SiteType, {'DD_L'}), :);
slRows = damage(strcmp(damage.SiteType, {'S_L'}), :);
shRows = damage(strcmp(damage.SiteType, {'S_H'}), :);
sRows = damage(strcmp(damage.SiteType, {'S_H'}) | strcmp(damage.SiteType, {'S_L'}), :);
%remove wetness index "inf" rows
sRows_noinf = damage(damage.WetnessIndex < Inf & strcmp(damage.SiteType, {'S_L'}) | strcmp(damage.SiteType, {'S_H'}), :);

% filter out NaNs for age of housing data
dah_nan = damage(~isnan(damage.AgeHousing) & strcmp(damage.SiteType, {'DA_H'}), :);
ddl_nan = damage(~isnan(damage.AgeHousing) & strcmp(damage.SiteType, {'DD_L'}), :);
sl_nan = damage(~isnan(damage.AgeHousing) & strcmp(damage.SiteType, {'S_L'}), :);
sh_nan = damage(~isnan(damage.AgeHousing) & strcmp(damage.SiteType, {'S_H'}), :);
s_nan = vertcat(sl_nan,sh_nan);
ddh_nan = damage(~isnan(damage.DrainageDensity) & strcmp(damage.SiteType, {'DD_H'}), :);
%%

% compute mean values for each magnitude of damage category for DA-H and
% S-H data
dahmag1 = dahRows(dahRows.MagnitudeOfDamage == 1, :); dahmag1mean = mean(dahmag1.DrainageArea);
dahmag2 = daRows(dahRows.MagnitudeOfDamage == 2, :); dahmag2mean = mean(dahmag2.DrainageArea);
dahmag3 = daRows(dahRows.MagnitudeOfDamage == 3, :); dahmag3mean = mean(dahmag3.DrainageArea);
dahmag4 = daRows(dahRows.MagnitudeOfDamage == 4, :); dahmag4mean = mean(dahmag4.DrainageArea);
dahmag5 = daRows(dahRows.MagnitudeOfDamage == 5, :); dahmag5mean = mean(dahmag5.DrainageArea);
shmag1 = shRows(shRows.MagnitudeOfDamage == 1, :); shmag1mean = mean(shmag1.Slope);
shmag2 = shRows(shRows.MagnitudeOfDamage == 2, :); shmag2mean = mean(shmag2.Slope);
shmag3 = shRows(shRows.MagnitudeOfDamage == 3, :); shmag3mean = mean(shmag3.Slope);
shmag4 = shRows(shRows.MagnitudeOfDamage == 4, :); shmag4mean = mean(shmag4.Slope);
shmag5 = shRows(shRows.MagnitudeOfDamage == 5, :); shmag5mean = mean(shmag5.Slope);

% scatterplot magnitude of damage vs DA and S, and plot mean magnitudes for DA-H and S-H measurements 
subplot(1,2,1), gscatter(daRows.DrainageArea, daRows.MagnitudeOfDamage, daRows.SiteType, [0 0.8 0.9;0 0.1 0.7], 'oo'), ylim([0.5 5.5]), yticks([1 2 3 4 5]), ...
    xlabel('Drainage Area'), ylabel('Magnitude of Damage'), axis square
hold on, scatter(dahmag1mean,1,100,'p','r'), scatter(dahmag2mean,2,100,'p','r'), scatter(dahmag3mean,3,100,'p','r'), scatter(dahmag4mean,4,100,'p','r'), scatter(dahmag5mean,5,100,'p','r'), hold off
subplot(1,2,2), gscatter(sRows.Slope, sRows.MagnitudeOfDamage, sRows.SiteType, [0 0.8 0.9;0 0.1 0.7], 'oo'), xlabel('Slope'), yticks([1 2 3 4 5]), ...
    ylim([0.5 5.5]), axis square
hold on, scatter(shmag1mean,1,100,'p','r'), scatter(shmag2mean,2,100,'p', 'r'), scatter(shmag3mean,3,100,'p','r'), scatter(shmag4mean,4,100,'p','r'), scatter(shmag5mean,5,100,'p','r'), hold off

%% get rows containing tilt and damage length measurements for each site, for H, L, and both H + L

daltilt = damage(~isnan(damage.DamageMetricsTilt_degreeFromVertical_) & strcmp(damage.SiteType, {'DA_L'}), :);
dahtilt = damage(~isnan(damage.DamageMetricsTilt_degreeFromVertical_) & strcmp(damage.SiteType, {'DA_H'}), :);
datilt = vertcat(daltilt,dahtilt);
ddltilt = damage(~isnan(damage.DamageMetricsTilt_degreeFromVertical_) & strcmp(damage.SiteType, {'DD_L'}), :);
ddhtilt = damage(~isnan(damage.DamageMetricsTilt_degreeFromVertical_) & strcmp(damage.SiteType, {'DD_H'}), :);
ddtilt = vertcat(ddltilt,ddhtilt);
sltilt = damage(~isnan(damage.DamageMetricsTilt_degreeFromVertical_) & strcmp(damage.SiteType, {'S_L'}), :);
shtilt = damage(~isnan(damage.DamageMetricsTilt_degreeFromVertical_) & strcmp(damage.SiteType, {'S_H'}), :);
stilt = vertcat(sltilt,shtilt);

da_length = vertcat(damage(~isnan(damage.DamageMetricsLength_m_) & strcmp(damage.SiteType, {'DA_L'}), :), ...
    damage(~isnan(damage.DamageMetricsLength_m_) & strcmp(damage.SiteType, {'DA_H'}), :));
dd_length = vertcat(damage(~isnan(damage.DamageMetricsLength_m_) & strcmp(damage.SiteType, {'DD_L'}), :), ...
    damage(~isnan(damage.DamageMetricsLength_m_) & strcmp(damage.SiteType, {'DD_H'}), :));
s_length = vertcat(damage(~isnan(damage.DamageMetricsLength_m_) & strcmp(damage.SiteType, {'S_L'}), :), ...
    damage(~isnan(damage.DamageMetricsLength_m_) & strcmp(damage.SiteType, {'S_H'}), :));

%% plot damage length and tilt for DA, DD, S, comparing H vs L sites for each

subplot(2,3,1), scatter(da_length.Low_1_High_2, da_length.DamageMetricsLength_m_), ylim([0 5]), xlim([0.75 2.25]), ...
    ylabel('Damage Length (m)'), xticklabels({'L', '', 'H'}), xlabel('DA'), axis square
hold on, scatter(1,dal_length_mean, 100, 'p', 'r'), scatter(2,dah_length_mean, 100, 'p', 'r'), hold off
subplot(2,3,2), scatter(dd_length.Low_1_High_2, dd_length.DamageMetricsLength_m_), xlim([0.75 2.25]), ylim([0 22]), xlabel('DD'), xticklabels({'L', '', 'H'}), axis square
hold on, scatter(1,ddl_length_mean, 100, 'p', 'r'), scatter(2,ddh_length_mean, 100, 'p', 'r'), hold off
subplot(2,3,3), scatter(s_length.Low_1_High_2, s_length.DamageMetricsLength_m_), xlim([0.75 2.25]), xlabel('S'), xticklabels({'L', '', 'H'}), axis square
hold on, scatter(1,sl_length_mean, 100, 'p', 'r'), scatter(2,sh_length_mean, 100, 'p', 'r'), hold off
subplot(2,3,4), scatter(datilt.Low_1_High_2, datilt.DamageMetricsTilt_degreeFromVertical_), xlim([0.75 2.25]), ylim([0 55]), xlabel('DA'), ylabel('Tilt (degrees from vertical)'), xticklabels({'L', '', 'H'}), axis square
hold on, scatter(1,daltilt_mean, 100, 'p', 'r'), scatter(2,dahtilt_mean, 100, 'p', 'r'), hold off
subplot(2,3,5), scatter(ddtilt.Low_1_High_2, ddtilt.DamageMetricsTilt_degreeFromVertical_), xlim([0.75 2.25]), xlabel('DD'), ylim([0 45]), xticklabels({'L', '', 'H'}), axis square
hold on, scatter(1,ddltilt_mean, 100, 'p', 'r'), scatter(2,ddhtilt_mean, 100, 'p', 'r'), hold off
subplot(2,3,6), scatter(stilt.Low_1_High_2, stilt.DamageMetricsTilt_degreeFromVertical_), xlim([0.75 2.25]), xlabel('S'), xticklabels({'L', '', 'H'}), ylim([0 13]), axis square
hold on, scatter(1,sltilt_mean, 100, 'p', 'r'), scatter(2,shtilt_mean, 100, 'p', 'r'), hold off

%% isolate rows containing damage to sidewalks (type 1) and house-walls (type 4)

da_walk = vertcat(damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'DA_L'}), :), ...
    damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'DA_H'}), :));
dd_walk = vertcat(damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'DD_L'}), :), ...
    damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'DD_H'}), :));
s_walk = vertcat(damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'S_L'}), :), ...
    damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'S_H'}), :));
dal_walk = damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'DA_L'}), :);
dah_walk = damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'DA_H'}), :);
sl_walk = damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'S_L'}), :);
sh_walk = damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'S_H'}), :);
ddl_walk = damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'DD_L'}), :);
ddh_walk = damage(damage.DamagedStructure == 1 & strcmp(damage.SiteType, {'DD_H'}), :);

da_wall = vertcat(damage(damage.DamagedStructure == 4 & strcmp(damage.SiteType, {'DA_L'}), :), ...
    damage(damage.DamagedStructure == 4 & strcmp(damage.SiteType, {'DA_H'}), :));
dd_wall = vertcat(damage(damage.DamagedStructure == 4 & strcmp(damage.SiteType, {'DD_L'}), :), ...
    damage(damage.DamagedStructure == 4 & strcmp(damage.SiteType, {'DD_H'}), :));
s_wall = vertcat(damage(damage.DamagedStructure == 4 & strcmp(damage.SiteType, {'S_L'}), :), ...
    damage(damage.DamagedStructure == 4 & strcmp(damage.SiteType, {'S_H'}), :));

%% for sidewalk data and wall-tilt data respectively, plot damage magnitude vs damage length for DD, DA, S, with combined H + L sites

L = {'DA','DD','S'};
subplot(1,2,1),
scatter(da_walk.MagnitudeOfDamage, da_walk.DamageMetricsLength_m_, 60, 'b'),
hold on
scatter(dd_walk.MagnitudeOfDamage, dd_walk.DamageMetricsLength_m_, 60, 'm'),
hold on
scatter(s_walk.MagnitudeOfDamage, s_walk.DamageMetricsLength_m_, 60, [0 0.8 0.7]),
hold off
ylim([-1 14]), xlim([0.5 5.5]), legend(L), title('Damage Type: Sidewalks'), xlabel('Magnitude Of Damage'), ...
    ylabel('Damage Length (m)'), xticks([1 2 3 4 5]), axis square,
subplot(1,2,2),
L = {'DA','DD','S'};
scatter(da_wall.MagnitudeOfDamage, da_wall.DamageMetricsTilt_degreeFromVertical_, 60, 'b'),
hold on
scatter(dd_wall.MagnitudeOfDamage, dd_wall.DamageMetricsTilt_degreeFromVertical_, 60, 'm'),
hold on
scatter(s_wall.MagnitudeOfDamage, s_wall.DamageMetricsTilt_degreeFromVertical_, 60, [0 0.8 0.7]), ...
    xlim([0.5 5.5]), legend(L), title('Damage Type: Walls'), xlabel('Magnitude Of Damage'), ...
    ylabel('Tilt (degrees from vertical'), xticks([1 2 3 4 5]), axis square
hold off

%% for above plot, reverse axes

L = {'DA','DD','S'};
subplot(1,2,1),
scatter(da_walk.DamageMetricsLength_m_, da_walk.MagnitudeOfDamage, 60, 'b'),
hold on
scatter(dd_walk.DamageMetricsLength_m_, dd_walk.MagnitudeOfDamage, 60, 'm'),
hold on
scatter(s_walk.DamageMetricsLength_m_, s_walk.MagnitudeOfDamage, 60, [0 0.8 0.7]),
hold off
xlim([-1 14]), ylim([0.5 5.5]), legend(L), title('Damage Type: Sidewalks'), xlabel('Damage Length (m)'), ...
    ylabel('Magnitude of Damage'), yticks([1 2 3 4 5]), axis square,
subplot(1,2,2),
L = {'DA','DD','S'};
scatter(da_wall.DamageMetricsTilt_degreeFromVertical_, da_wall.MagnitudeOfDamage, 60, 'b'),
hold on
scatter(dd_wall.DamageMetricsTilt_degreeFromVertical_, dd_wall.MagnitudeOfDamage, 60, 'm'),
hold on
scatter(s_wall.DamageMetricsTilt_degreeFromVertical_, s_wall.MagnitudeOfDamage, 60, [0 0.8 0.7]), ...
    ylim([0.5 5.5]), legend(L), title('Damage Type: Walls'), ylabel('Magnitude Of Damage'), ...
    xlabel('Tilt (degrees from vertical)'), yticks([1 2 3 4 5]), xlim([-5 45]), axis square
hold off

%% plot histograms for each group for magnitude of damage, measured length, and measured tilt

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
    histogram(damage.MagnitudeOfDamage(strcmp(damage.Person_s_Collecting, grps{i})), 5)
    title([grps{i} ': Damage magnitude'])
end

%% plot CDU map unfiltered, and with min and mode filters

subplot(1,3,1), imagesc(cdu), colorbar, title('cdu unfiltered'), axis square
subplot(1,3,2), imagesc(cdu_avg), colorbar, title('cdu min filter'), axis square
subplot(1,3,3), imagesc(cdu_mode), colorbar, title('cdu mode filter'), axis square

%% compute means for the 4 most common CDU categories for wall tilt and sidewalk crack length
% walls
wallCDU4 = damage((damage.DamagedStructure == 4) & (damage.CDUmodenearest == 4), :);
wallCDU5 = damage((damage.DamagedStructure == 4) & (damage.CDUmodenearest == 5), :);
wallCDU6 = damage((damage.DamagedStructure == 4) & (damage.CDUmodenearest == 6), :);
wallCDU7 = damage((damage.DamagedStructure == 4) & (damage.CDUmodenearest == 7), :);
wallCDU4mean = nanmean(wallCDU4.DamageMetricsTilt_degreeFromVertical_);
wallCDU5mean = nanmean(wallCDU5.DamageMetricsTilt_degreeFromVertical_);
wallCDU6mean = nanmean(wallCDU6.DamageMetricsTilt_degreeFromVertical_);
wallCDU7mean = nanmean(wallCDU7.DamageMetricsTilt_degreeFromVertical_);
%sidewalks
walkCDU4 = damage((damage.DamagedStructure == 1) & (damage.CDUmodenearest == 4), :);
walkCDU5 = damage((damage.DamagedStructure == 1) & (damage.CDUmodenearest == 5), :);
walkCDU6 = damage((damage.DamagedStructure == 1) & (damage.CDUmodenearest == 6), :);
walkCDU7 = damage((damage.DamagedStructure == 1) & (damage.CDUmodenearest == 7), :);
walkCDU4mean = nanmean(walkCDU4.DamageMetricsLength_m_);
walkCDU5mean = nanmean(walkCDU5.DamageMetricsLength_m_);
walkCDU6mean = nanmean(walkCDU6.DamageMetricsLength_m_);
walkCDU7mean = nanmean(walkCDU7.DamageMetricsLength_m_);
%% CDU mode nearest - walls
figure
L = {'DA','DD','S'};
scatter(da_wall.CDUmodenearest, da_wall.DamageMetricsTilt_degreeFromVertical_, 'b'),
hold on
scatter(dd_wall.CDUmodenearest, dd_wall.DamageMetricsTilt_degreeFromVertical_, 'g'),
hold on
scatter(s_wall.CDUmodenearest, s_wall.DamageMetricsTilt_degreeFromVertical_, 'c'),
hold on, scatter(4,wallCDU4mean,100,'p','r'),scatter(5,wallCDU5mean,100,'p','r'),scatter(6,wallCDU6mean,100,'p','r'),scatter(7,wallCDU7mean,100,'p','r')
xlim([0.5 8.5]), legend(L), title('Walls: CDU score'), xlabel('CDU'), ...
    ylabel('Tilt(degrees from vertical'), xlim([0.5 8.5]), axis square
hold off
%% CDU mode nearest - sidewalks
figure
L = {'DA','DD','S'};
scatter(da_walk.CDUmodenearest, da_walk.DamageMetricsLength_m_, 'b'),
hold on 
scatter(dd_walk.CDUmodenearest, dd_walk.DamageMetricsLength_m_, 'g'),
hold on
scatter(s_walk.CDUmodenearest, s_walk.DamageMetricsLength_m_, 'c'),
hold on, scatter(4,walkCDU4mean,100,'p','r'),scatter(5,walkCDU5mean,100,'p','r'),scatter(6,walkCDU6mean,100,'p','r'),scatter(7,walkCDU7mean,100,'p','r')
xlim([0.5 8.5]), legend(L), title('Sidewalks: CDU score'), xlabel('CDU'), ...
    ylabel('Damage Length (m)'), xlim([0.5 8.5]), ylim([-1 15]), axis square
hold off

%% CDU: mag damage all vs cdu mode filter

L = {'DA','DD','S'};
scatter(daRows.CDUmodenearest, daRows.MagnitudeOfDamage, 80, 'b'),
hold on
scatter(ddRows.CDUmodenearest, ddRows.MagnitudeOfDamage, 80, '.', 'r'),
hold on
scatter(sRows.CDUmodenearest, sRows.MagnitudeOfDamage, 170, 'MarkerEdgeColor', [0 0.8 0.2]),
xlim([0.5 8.5]), legend(L), title('Magnitude: CDU mode filter'), xlabel('CDU'), ...
    ylabel('Mag Dam'), xlim([0.5 8.5]), ylim([0 6]), axis square
hold off

%% CDU: plot magnitude of damage vs CDU

subplot(1,3,1), gscatter(daRows.MagnitudeOfDamage, daRows.CDU, daRows.SiteType), ylabel('CDU score'), axis square
subplot(1,3,2), gscatter(sRows.MagnitudeOfDamage, sRows.CDU, sRows.SiteType), xlabel('Mag Dam'), axis square
subplot(1,3,3), gscatter(ddRows.MagnitudeOfDamage, ddRows.CDU, ddRows.SiteType), axis square

%% CDU: plot damage length (*all damage, not just sidewalks) vs CDU

subplot(1,3,1), gscatter(daRows.DamageMetricsLength_m_, daRows.CDU, daRows.SiteType), ylabel('CDU score'), axis square
subplot(1,3,2), gscatter(sRows.DamageMetricsLength_m_, sRows.CDU, sRows.SiteType), xlabel('Length (m), all damaged structure types combined'), axis square
subplot(1,3,3), gscatter(ddRows.DamageMetricsLength_m_, ddRows.CDU, ddRows.SiteType), axis square
