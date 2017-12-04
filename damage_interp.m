
% **** add this folder to your matlab path: geotiffinterp_documentation

%get coordinates of our DD damage photos
DDdamage = xlsread('DDcoord.xlsx');

%get values for metrics computed from DEM:
DDinterp = geotiffinterp('tmp/drainage_density_fixed.tif',DDdamage(:,1),DDdamage(:,2), 'nearest');
Sinterp = geotiffinterp('tmp/slope_gauss.tif',DDdamage(:,1),DDdamage(:,2));
DAinterp = geotiffinterp('tmp/drainage_area_mdf.tif',DDdamage(:,1),DDdamage(:,2));
WIinterp = geotiffinterp('tmp/wetness_index.tif',DDdamage(:,1),DDdamage(:,2));