
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
%% filter CDU: mean
cdu = GRIDobj('resources/cdu3final.tif');
cdu.Z(cdu.Z == 0) = NaN;
n_px2 = 3;
cdu_avg = cdu;
cdu_avg.Z = fillmissing(cdu_avg.Z, 'constant', 9999);
cdu_avg.Z = ordfilt2(cdu_avg.Z, 1, ones(n_px2, n_px2));
cdu_avg.Z(cdu_avg.Z == 9999) = NaN;
cdu_avg.Z(cdu_avg.Z == 0) = NaN;
 
good = GRIDobj('tmp/drainage_area_mdf.tif');
cdu_avg.refmat = good.refmat;
cdu_avg.georef = good.georef;
GRIDobj2geotiff(cdu_avg, 'resources/cdu_avg.tif');
%%
CDUnofilter = GRIDobj('resources/cdu3final.tif');
CDUnofilter.Z(CDUnofilter.Z == 0) = NaN;
CDUnofilter.refmat = good.refmat;
CDUnofilter.georef = good.georef;
GRIDobj2geotiff(CDUnofilter, 'resources/CDUnofilter.tif');
%% filter CDU: mode


%%
%CDU.Z=changem(CDU.Z,NaN,0);
CDU = GRIDobj('resources/cdu_avg.tif');
DD_interp_CDU = geotiffinterp('resources/cdu_avg.tif',DD_damage(:,1),DD_damage(:,2));
DA_interp_CDU = geotiffinterp('resources/cdu_avg.tif',DA_damage(:,1),DA_damage(:,2));
S1_interp_CDU = geotiffinterp('resources/cdu_avg.tif',S1_damage(:,1),S1_damage(:,2));
S2_interp_CDU = geotiffinterp('resources/cdu_avg.tif',S2_damage(:,1),S2_damage(:,2));
%%
good = GRIDobj('tmp/drainage_area_mdf.tif');
cdu_mode.refmat = good.refmat;
cdu_mode.georef = good.georef;

GRIDobj2geotiff(cdu_mode, 'resources/cdu_mode.tif');

%%
% using mode and linear interpolation:
DD_interp_CDUmode = geotiffinterp('resources/cdu_mode.tif',DD_damage(:,1),DD_damage(:,2));
DA_interp_CDUmode = geotiffinterp('resources/cdu_mode.tif',DA_damage(:,1),DA_damage(:,2));
S1_interp_CDUmode = geotiffinterp('resources/cdu_mode.tif',S1_damage(:,1),S1_damage(:,2));
S2_interp_CDUmode = geotiffinterp('resources/cdu_mode.tif',S2_damage(:,1),S2_damage(:,2));

%% using mode and nearest neighbor: 
DD_interp_CDUmodenearest = geotiffinterp('resources/cdu_mode.tif',DD_damage(:,1),DD_damage(:,2), 'nearest');
DA_interp_CDUmodenearest = geotiffinterp('resources/cdu_mode.tif',DA_damage(:,1),DA_damage(:,2), 'nearest');
S1_interp_CDUmodenearest = geotiffinterp('resources/cdu_mode.tif',S1_damage(:,1),S1_damage(:,2), 'nearest');
S2_interp_CDUmodenearest = geotiffinterp('resources/cdu_mode.tif',S2_damage(:,1),S2_damage(:,2), 'nearest');

%% CDU min filter, nearest (using cdu_avg
DD_interp_CDUminnearest = geotiffinterp('resources/cdu_avg.tif',DD_damage(:,1),DD_damage(:,2), 'nearest');
DA_interp_CDUminnearest = geotiffinterp('resources/cdu_avg.tif',DA_damage(:,1),DA_damage(:,2), 'nearest');
S1_interp_CDUminnearest = geotiffinterp('resources/cdu_avg.tif',S1_damage(:,1),S1_damage(:,2), 'nearest');
S2_interp_CDUminnearest = geotiffinterp('resources/cdu_avg.tif',S2_damage(:,1),S2_damage(:,2), 'nearest');

%% CDU no filter, nearest
DD_interp_CDUnofilter = geotiffinterp('resources/CDUnofilter.tif',DD_damage(:,1),DD_damage(:,2), 'nearest');
DA_interp_CDUnofilter = geotiffinterp('resources/CDUnofilter.tif',DA_damage(:,1),DA_damage(:,2), 'nearest');
S1_interp_CDUnofilter = geotiffinterp('resources/CDUnofilter.tif',S1_damage(:,1),S1_damage(:,2), 'nearest');
S2_interp_CDUnofilter = geotiffinterp('resources/CDUnofilter.tif',S2_damage(:,1),S2_damage(:,2), 'nearest');
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