
% **** add this folder to your matlab path: geotiffinterp_documentation

%% get coordinates for each group's damage photos
DD_damage = xlsread('DD_damage.xlsx'); 
DA_damage = xlsread('DA_damage.xlsx'); 
S1_damage = xlsread('S1_damage.xlsx'); 
S2_damage = xlsread('S2_damage.xlsx'); 

%% fix refmat and georef for impervious_filtered
good = GRIDobj('tmp/drainage_area_mdf.tif');
bad = GRIDobj('resources/impervious_filtered.tif');
bad.refmat = good.refmat;
bad.georef = good.georef;
GRIDobj2geotiff(bad, 'resources/impervious_filtered.tif');
bad2 = GRIDobj('resources/housing_age_filtered.tif');
bad2.refmat = good.refmat;
bad2.georef = good.georef;
GRIDobj2geotiff(bad2, 'resources/housing_age_filtered.tif');

%% for Drainage Density group's damage measurement locations, get values of metrics from DEM: 
DD_interp_DD = geotiffinterp('resources/drainage_density_filtered.tif',DD_damage(:,1),DD_damage(:,2));
DD_interp_slope = geotiffinterp('resources/slope_filtered.tif',DD_damage(:,1),DD_damage(:,2));
DD_interp_DA = geotiffinterp('resources/drainage_area_mdf_filtered.tif',DD_damage(:,1),DD_damage(:,2));
DD_interp_WI = geotiffinterp('resources/wetness_index_filtered.tif',DD_damage(:,1),DD_damage(:,2));
DD_interp_imp = geotiffinterp('resources/impervious_filtered.tif',DD_damage(:,1),DD_damage(:,2));
DD_interp_DD_nearest = geotiffinterp('resources/drainage_density_filtered.tif',DD_damage(:,1),DD_damage(:,2),'nearest');
DD_interp_housing = geotiffinterp('resources/housing_age_filtered.tif',DD_damage(:,1),DD_damage(:,2),'nearest');
%% for Drainage Area group's damage measurement locations, get values of
%metrics from DEM:
DA_interp_DA = geotiffinterp('resources/drainage_area_mdf_filtered.tif',DA_damage(:,1),DA_damage(:,2));
DA_interp_DD = geotiffinterp('resources/drainage_density_filtered.tif',DA_damage(:,1),DA_damage(:,2));
DA_interp_slope = geotiffinterp('resources/slope_filtered.tif',DA_damage(:,1),DA_damage(:,2));
DA_interp_WI = geotiffinterp('resources/wetness_index_filtered.tif',DA_damage(:,1),DA_damage(:,2));
DA_interp_imp = geotiffinterp('resources/impervious_filtered.tif',DA_damage(:,1),DA_damage(:,2));
DA_interp_DD_nearest = geotiffinterp('resources/drainage_density_filtered.tif',DA_damage(:,1),DA_damage(:,2),'nearest');
DA_interp_housing = geotiffinterp('resources/housing_age_filtered.tif',DA_damage(:,1),DA_damage(:,2),'nearest');
%% for Slope1 group:
S1_interp_slope = geotiffinterp('resources/slope_filtered.tif',S1_damage(:,1),S1_damage(:,2));
S1_interp_DA = geotiffinterp('resources/drainage_area_mdf_filtered.tif',S1_damage(:,1),S1_damage(:,2));
S1_interp_DD = geotiffinterp('resources/drainage_density_filtered.tif',S1_damage(:,1),S1_damage(:,2));
S1_interp_WI = geotiffinterp('resources/wetness_index_filtered.tif',S1_damage(:,1),S1_damage(:,2));
S1_interp_imp = geotiffinterp('resources/impervious_filtered.tif',S1_damage(:,1),S1_damage(:,2));
S1_interp_DD_nearest = geotiffinterp('resources/drainage_density_filtered.tif',S1_damage(:,1),S1_damage(:,2),'nearest');
S1_interp_housing = geotiffinterp('resources/housing_age_filtered.tif',S1_damage(:,1),S1_damage(:,2),'nearest');
%% for Slope2 group:
S2_interp_slope = geotiffinterp('resources/slope_filtered.tif',S2_damage(:,1),S2_damage(:,2));
S2_interp_DA = geotiffinterp('resources/drainage_area_mdf_filtered.tif',S2_damage(:,1),S2_damage(:,2));
S2_interp_DD = geotiffinterp('resources/drainage_density_filtered.tif',S2_damage(:,1),S2_damage(:,2));
S2_interp_WI = geotiffinterp('resources/wetness_index_filtered.tif',S2_damage(:,1),S2_damage(:,2));
S2_interp_imp = geotiffinterp('resources/impervious_filtered.tif',S2_damage(:,1),S2_damage(:,2));
S2_interp_DD_nearest = geotiffinterp('resources/drainage_density_filtered.tif',S2_damage(:,1),S2_damage(:,2));
S2_interp_housing = geotiffinterp('resources/housing_age_filtered.tif',S2_damage(:,1),S2_damage(:,2),'nearest');
%%