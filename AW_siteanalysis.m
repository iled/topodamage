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
damage_fn = 'TopoDataAnalysis_AW.csv';
opts = detectImportOptions(damage_fn);
damage = readtable(damage_fn, opts);
% rename columns
damage.Properties.VariableNames{'DamagedStructure_1_sidewalk_2_road_3_fence_4_house_wall_5'} = 'DamagedStructure';
damage.Properties.VariableNames{'MagnitudeOfDamage'} = 'DamageMagnitude';
damage.Properties.VariableNames{'MaterialType_1_concrete_2_brick_3_asphalt_4_stone_5_wood_'} = 'MaterialType';

%% ignore DA sites for now; they have issues with the coordinates
DAs = strcmp(damage.SiteType, 'DA_L') + strcmp(damage.SiteType, 'DA_H');
damage = damage(~DAs, :);

%% get topographic metric values from maps
damage.DrainageDensity = geotiffinterp('resources/drainage_density_fixed.tif', damage.Lat, damage.Long);
damage.DrainageArea = geotiffinterp('resources/drainage_area_mdf.tif', damage.Lat, damage.Long);
[damage.Slope, x, y] = geotiffinterp('resources/slope_gauss.tif', damage.Lat, damage.Long);
damage.WetnessIndex = geotiffinterp('resources/wetness_index.tif', damage.Lat, damage.Long);

%% PLOT: sample sites over DEM
imagesc(DEM)
hold on
%plot(x, y, 'r*')
gscatter(x, y, damage.SiteType, 'mgrbcy', '+', 8, 'on')
%legend('DD High', 'DD Low', 'Slope High', 'Slope Low', ...
        %'Location', 'northwestoutside')
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

%% Group stats: Damaged structure and magnitude per site
% grpstats requires the Statistics and Machine Learning Toolbox
grpstats(damage, 'SiteType', {'min', 'max', 'mean', @mode}, 'DataVars', ... 
    {'DamagedStructure', 'DamageMagnitude'})

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


%% load damage data *** Julio this is where I started using my version of the topo spreadsheet, and diverged from your siteanalysis.m file

damage_fn = 'TopoDataAnalysis_AW.csv';
%opts = detectImportOptions(damage_fn);
damage = readtable(damage_fn);
% rename columns
%damage.Properties.VariableNames{'DamagedStructure_1_sidewalk_2_road_3_fence_4_house_wall_5_'} = 'DamagedStructure';
%damage.Properties.VariableNames{'MagnitudeOfDamage'} = 'DamageMagnitude';
%damage.Properties.VariableNames{'MaterialType_1_concrete_2_brick_3_asphalt_4_stone_5_wood_'} = 'MaterialType';
%damage.Properties.VariableNames{'DrainageArea'} = 'DrainageArea';
%damage.Properties.VariableNames{'DrainageDensity'} = 'DrainageDensity';
%damage.Properties.VariableNames{'Slope'} = 'Slope';

%% calculating correlation coefficients for DD

DDL = damage(strcmp(damage.SiteType, {'DD_L'}), :);
DDL_R = corrcoef(DDL.MagnitudeOfDamage,DDL.DrainageArea);
DDH = damage(strcmp(damage.SiteType, {'DD_H'}), :);
DDH_R = corrcoef(DDH.MagnitudeOfDamage,DDH.DrainageArea);
DDLcorr2 = corr2(DDL.MagnitudeOfDamage,DDL.DrainageArea); % -0.0429
DDHcorr2 = corr2(DDH.MagnitudeOfDamage,DDH.DrainageArea); % -0.3606
DD = damage(strcmp(damage.SiteType, {'DD_H'}) | strcmp(damage.SiteType, {'DD_L'}), :);
DDcorr2 = corr2(DD.MagnitudeOfDamage,DD.DrainageArea); % 0.2341

%% corr coeffs for DA

DAL = damage(54:73,:);
DAH = damage(73:103,:);
DA = damage(54:103,:);
DALc2 = corr2(DAL.DamageMagnitude,DAL.DrainageArea); % 0.2641
DAHc2 = corr2(DAH.DamageMagnitude,DAH.DrainageArea); % -0.1125
DAc2 = corr2(DA.DamageMagnitude,DA.DrainageArea); % 0.1193
DALc2_slope = corr2(DAL.DamageMagnitude,DAL.Slope); % 0.4354
DAHc2_slope = corr2(DAH.DamageMagnitude,DAH.Slope); % -0.0407
DAc2_slope = corr2(DA.DamageMagnitude,DA.Slope); % 0.2843
%DALc2_dd = corr2(DAL.DamageMagnitude,DAL.DrainageDensity) % 
%DAHc2_dd = corr2(DAH.DamageMagnitude,DAH.DrainageDensity) % 
%DAc2_dd = corr2(DA.DamageMagnitude,DA.DrainageDensity) % 
DALc2_imp = corr2(DAL.DamageMagnitude,DAL.ImperviousCover); % -0.2310
DAHc2_imp = corr2(DAH.DamageMagnitude,DAH.ImperviousCover); % 0.2789
DAc2_imp = corr2(DA.DamageMagnitude,DA.ImperviousCover); % 0.3774
DALc2_WI = corr2(DAL.DamageMagnitude,DAL.WetnessIndex); % 0.2240
DAHc2_WI = corr2(DAH.DamageMagnitude,DAH.WetnessIndex); % -0.0309
DAc2_WI = corr2(DA.DamageMagnitude,DA.WetnessIndex); % 0.2091

%% more corr coeffs

DALc2_DD = corr2(DAL.MagnitudeOfDamage,DAL.DrainageDensity)
DAHc2_DD = corr2(DAH.MagnitudeOfDamage,DAH.DrainageDensity)
DAc2_DD = corr2(DA.MagnitudeOfDamage,DA.DrainageDensity)

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
slinf = slRows(1:28,:);
%sinf = sRows(1:43,:) | sRows(55:84,:)
sRows_noinf = damage(damage.WetnessIndex < Inf & strcmp(damage.SiteType, {'S_L'}) | strcmp(damage.SiteType, {'S_H'}), :);

% filter out NaNs for age of housing data
dah_nan = damage(~isnan(damage.AgeHousing) & strcmp(damage.SiteType, {'DA_H'}), :);
ddl_nan = damage(~isnan(damage.AgeHousing) & strcmp(damage.SiteType, {'DD_L'}), :);
sl_nan = damage(~isnan(damage.AgeHousing) & strcmp(damage.SiteType, {'S_L'}), :);
sh_nan = damage(~isnan(damage.AgeHousing) & strcmp(damage.SiteType, {'S_H'}), :);
%s_nan = damage(~isnan(damage.AgeHousing & strcmp(damage.SiteType, {'S_L'}) | strcmp(damage.SiteType, {'S_H'}), :))
s_nan = vertcat(sl_nan,sh_nan);
ddh_nan = damage(~isnan(damage.DrainageDensity) & strcmp(damage.SiteType, {'DD_H'}), :);

%% for DA group damage only, plot magnitude of damage against the interpolated value for each topo metric

subplot(1,4,1), gscatter(daRows.MagnitudeOfDamage, daRows.DrainageArea, daRows.SiteType), axis square
subplot(1,4,2), gscatter(daRows.MagnitudeOfDamage, daRows.Slope, daRows.SiteType), axis square
subplot(1,4,3), gscatter(daRows.MagnitudeOfDamage, daRows.WetnessIndex, daRows.SiteType), axis square
subplot(1,4,4), gscatter(daRows.MagnitudeOfDamage, daRows.ImperviousCover, daRows.SiteType), axis square
%subplot(1,5,5), gscatter(daRows.MagnitudeOfDamage, daRows.AgeHousing, daRows.SiteType), axis square

%% plot magnitude of damage against the respective isolated metric for each pair of sites

subplot(1,3,1), gscatter(daRows.MagnitudeOfDamage, daRows.DrainageArea, daRows.SiteType), xlabel('Magnitude Damage'), ylabel('Drainage Area'), axis square
subplot(1,3,2), gscatter(sRows.MagnitudeOfDamage, sRows.Slope, sRows.SiteType), xlabel('Magnitude Damage'), ylabel('Slope'), ...
    title('Magnitude Damage vs isolated topographic metric for each pair of sites'), axis square
subplot(1,3,3), gscatter(ddRows.MagnitudeOfDamage, ddRows.DrainageDensity, ddRows.SiteType), xlabel('Magnitude Damage'), ylabel('Drainage Density'), axis square

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

subplot(2,3,1), scatter(da_length.Low_1_High_2, da_length.DamageMetricsLength_m_), xlim([0.75 2.25]), ...
    ylabel('Damage Length (m)'), set(gca,'xticklabel',{[]}),xlabel('(L)    DA     (H)'), axis square
subplot(2,3,2), scatter(dd_length.Low_1_High_2, dd_length.DamageMetricsLength_m_), xlim([0.75 2.25]), set(gca,'xticklabel',{[]}),xlabel('(L)    DD     (H)'), axis square
subplot(2,3,3), scatter(s_length.Low_1_High_2, s_length.DamageMetricsLength_m_), xlim([0.75 2.25]), set(gca,'xticklabel',{[]}),xlabel('(L)    S     (H)'), axis square
subplot(2,3,4), scatter(datilt.Low_1_High_2, datilt.DamageMetricsTilt_degreeFromVertical_), xlim([0.75 2.25]), ylabel('Tilt (degrees from vertical)')
subplot(2,3,5), scatter(ddtilt.Low_1_High_2, ddtilt.DamageMetricsTilt_degreeFromVertical_), xlim([0.75 2.25]), ylim([0 45])
subplot(2,3,6), scatter(stilt.Low_1_High_2, stilt.DamageMetricsTilt_degreeFromVertical_), xlim([0.75 2.25]),

%% compare tilt for each pair of sites vs their respective topo metrics ***** I think axes should be switched...?*****

subplot(1,3,1), gscatter(datilt.DamageMetricsTilt_degreeFromVertical_, datilt.DrainageArea, datilt.SiteType), axis square
subplot(1,3,2), gscatter(ddtilt.DamageMetricsTilt_degreeFromVertical_, ddtilt.DrainageDensity, ddtilt.SiteType), axis square
subplot(1,3,3), gscatter(stilt.DamageMetricsTilt_degreeFromVertical_, stilt.Slope, stilt.SiteType), axis square

%% plot tilt for each pair of sites comparing H vs L

subplot(1,3,1), scatter(datilt.Low_1_High_2, datilt.DamageMetricsTilt_degreeFromVertical_), xlim([0.75 2.25]), ylabel('degrees from vertical'), set(gca,'xticklabel',{[]}),xlabel('(L)            DA             (H)'), axis square
subplot(1,3,2), scatter(ddtilt.Low_1_High_2, ddtilt.DamageMetricsTilt_degreeFromVertical_), xlim([0.75 2.25]), title('Tilt'), set(gca,'xticklabel',{[]}),xlabel('(L)            DD             (H)'), axis square
subplot(1,3,3), scatter(stilt.Low_1_High_2, stilt.DamageMetricsTilt_degreeFromVertical_), xlim([0.75 2.25]), set(gca,'xticklabel',{[]}),xlabel('(L)            S             (H)'), axis square

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

%% for sidewalk data only, plot damage magnitude vs damage length for DD, DA, S, with combined H + L sites

L = {'DA','DD','S'};
scatter(da_walk.MagnitudeOfDamage, da_walk.DamageMetricsLength_m_, 'b'),
hold on
scatter(dd_walk.MagnitudeOfDamage, dd_walk.DamageMetricsLength_m_, 'r'),
hold on
scatter(s_walk.MagnitudeOfDamage, s_walk.DamageMetricsLength_m_, 'c'),
ylim([-1 14]), xlim([0.5 5.5]), legend(L), title('Measured vs assigned damage: Sidewalks'), xlabel('Magnitude Of Damage'), ...
    ylabel('Damage Length (m)'), xticks([1 2 3 4 5])
hold off

%% for house-wall data only, plot damage magnitude vs damage length for DD, DA, S, with combined H + L sites

L = {'DA','DD','S'};
scatter(da_wall.MagnitudeOfDamage, da_wall.DamageMetricsTilt_degreeFromVertical_, 'b'),
hold on
scatter(dd_wall.MagnitudeOfDamage, dd_wall.DamageMetricsTilt_degreeFromVertical_, 'r'),
hold on
scatter(s_wall.MagnitudeOfDamage, s_wall.DamageMetricsTilt_degreeFromVertical_, 'c'), ...
    xlim([0.5 5.5]), legend(L), title('Measured vs assigned damage: Walls'), xlabel('Magnitude Of Damage'), ...
    ylabel('Tilt (degrees from vertical'), xticks([1 2 3 4 5])
hold off

%% scatter plot of magnitude damage for all groups combined, comparing L vs H sites

scatter(da_length.Low_1_High_2, da_length.DamageMetricsLength_m_, 'b'), xlim([0.75 2.25]), ...
    ylabel('Damage Length (m)'), set(gca,'xticklabel',{[]}),
hold on
scatter(dd_length.Low_1_High_2, dd_length.DamageMetricsLength_m_, 'r'),
hold on
scatter(s_length.Low_1_High_2, s_length.DamageMetricsLength_m_, 'c')
hold off
%scatter(dd.MagnitudeOfDamage, damage.SiteType),
%hold on
%scatter(ddRows.Low_1_High_2, ddRows.MagnitudeOfDamage, 'r'),
%hold on
%scatter(sRows.Low_1_High_2, sRows.MagnitudeOfDamage, 'c'),
%hold off

%% same as above but for H and L each group - bias results by knowing which site is which? ** looks like you can't tell, not enough mag 4&5 measured

scatter(dal_walk.MagnitudeOfDamage, dal_walk.DamageMetricsLength_m_, 'b'),
hold on
scatter(dah_walk.MagnitudeOfDamage, dah_walk.DamageMetricsLength_m_, 'b', '^'),
hold on
scatter(ddl_walk.MagnitudeOfDamage, ddl_walk.DamageMetricsLength_m_, 'r'),
hold on
scatter(ddh_walk.MagnitudeOfDamage, ddh_walk.DamageMetricsLength_m_, 'r', '^'),
hold on
scatter(sl_walk.MagnitudeOfDamage, sl_walk.DamageMetricsLength_m_, 'c'),
hold on
scatter(sh_walk.MagnitudeOfDamage, sh_walk.DamageMetricsLength_m_, 'c', '^')
ylim([-1 14]), xlim([0.5 5.5])
hold off

%% CDU: for sidewalks only, plot damage length vs CDU, all groups

L = {'DA','DD','S'};
scatter(da_walk.CDU, da_walk.DamageMetricsLength_m_, 'b'),
hold on
scatter(dd_walk.CDU, dd_walk.DamageMetricsLength_m_, 'r'),
hold on
scatter(s_walk.CDU, s_walk.DamageMetricsLength_m_, 'c'),
xlim([0.5 8.5]), legend(L), title('Measured damage vs CDU: Sidewalks'), xlabel('CDU'), ...
    ylabel('Damage Length (m)'), xlim([2.5 7.5]), ylim([-1 15]), axis square
hold off

%% CDU: for house-walls only, plot tilt vs CDU, all groups***cdu_avg (min) (linear)

L = {'DA','DD','S'};
scatter(datilt.CDU, datilt.DamageMetricsTilt_degreeFromVertical_, 'b'),
hold on
scatter(ddtilt.CDU, ddtilt.DamageMetricsTilt_degreeFromVertical_, 'r'),
hold on
scatter(stilt.CDU, stilt.DamageMetricsTilt_degreeFromVertical_, 'c'),
xlim([0.5 8.5]), legend(L), title('Measured damage vs CDU: Walls'), xlabel('CDU'), ...
    ylabel('Tilt(degrees from vertical'), xlim([2.5 7.5]), axis square
hold off

%% CDU mode - changed by julio, need to change back
L = {'DA','DD','S'};
scatter(DA_interp_CDUmode(tilted), datilt.DamageMetricsTilt_degreeFromVertical_, 'b'),
hold on
scatter(DD_interp_CDUmode(tilted), ddtilt.DamageMetricsTilt_degreeFromVertical_, 'r'),
scatter(S_interp_CDUmode(tilted), stilt.DamageMetricsTilt_degreeFromVertical_, 'c'),
%xlim([0.5 8.5]), legend(L), title('Measured damage vs CDU: Walls'), xlabel('CDU'), ...
%    ylabel('Tilt(degrees from vertical'), xlim([2.5 7.5]), axis square
hold off
shgf

%% CDU min nearest - tilt
L = {'DA','DD','S'};
scatter(datilt.CDUminnearest, datilt.DamageMetricsTilt_degreeFromVertical_, 'b'),
hold on
scatter(ddtilt.CDUminnearest, ddtilt.DamageMetricsTilt_degreeFromVertical_, 'r'),
hold on
scatter(stilt.CDUminnearest, stilt.DamageMetricsTilt_degreeFromVertical_, 'c'),
xlim([0.5 8.5]), legend(L), title('Walls: CDU min filter'), xlabel('CDU'), ...
    ylabel('Tilt (degrees from vertical)'), xlim([0.5 8.5]), axis square
hold off

%% CDU min nearest - sidewalks
L = {'DA','DD','S'};
scatter(da_walk.CDUminnearest, da_walk.DamageMetricsLength_m_, 'b'),
hold on
scatter(dd_walk.CDUminnearest, dd_walk.DamageMetricsLength_m_, 'r'),
hold on
scatter(s_walk.CDUminnearest, s_walk.DamageMetricsLength_m_, 'c'),
xlim([0.5 8.5]), legend(L), title('Sidewalks: CDU min filter'), xlabel('CDU'), ...
    ylabel('Damage Length (m)'), xlim([0.5 8.5]), ylim([-1 15]), axis square
hold off
%% CDU mode nearest - plot tilt
L = {'DA','DD','S'};
scatter(datilt.CDUmodenearest, datilt.DamageMetricsTilt_degreeFromVertical_, 'b'),
hold on
scatter(ddtilt.CDUmodenearest, ddtilt.DamageMetricsTilt_degreeFromVertical_, 'r'),
hold on
scatter(stilt.CDUmodenearest, stilt.DamageMetricsTilt_degreeFromVertical_, 'c'),
xlim([0.5 8.5]), legend(L), title('Walls: CDU mode filter'), xlabel('CDU'), ...
    ylabel('Tilt(degrees from vertical'), xlim([0.5 8.5]), axis square
hold off

%% CDU mode nearest - plot sidewalks
L = {'DA','DD','S'};
scatter(da_walk.CDUmodenearest, da_walk.DamageMetricsLength_m_, 'b'),
hold on
scatter(dd_walk.CDUmodenearest, dd_walk.DamageMetricsLength_m_, 'r'),
hold on
scatter(s_walk.CDUmodenearest, s_walk.DamageMetricsLength_m_, 'c'),
xlim([0.5 8.5]), legend(L), title('Sidewalks: CDU mode filter'), xlabel('CDU'), ...
    ylabel('Damage Length (m)'), xlim([0.5 8.5]), ylim([-1 15]), axis square
hold off

%% CDU: plot magnitude of damage vs CDU

subplot(1,3,1), gscatter(daRows.MagnitudeOfDamage, daRows.CDU, daRows.SiteType), axis square
subplot(1,3,2), gscatter(sRows.MagnitudeOfDamage, sRows.CDU, sRows.SiteType), axis square
subplot(1,3,3), gscatter(ddRows.MagnitudeOfDamage, ddRows.CDU, ddRows.SiteType), axis square

%% CDU: plot damage length (*all damage, not just sidewalks) vs CDU

subplot(1,3,1), gscatter(daRows.DamageMetricsLength_m_, daRows.CDU, daRows.SiteType), axis square
subplot(1,3,2), gscatter(sRows.DamageMetricsLength_m_, sRows.CDU, sRows.SiteType), axis square
subplot(1,3,3), gscatter(ddRows.DamageMetricsLength_m_, ddRows.CDU, ddRows.SiteType), axis square
