%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topographic Analysis - HW5
% Sam Mark, Arielle Woods, Julio Caineta
% Main file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generation of maps is performed within the script file generate_maps.m

% site selection is performed within the script file site_selection.m

%% load data
DEM = GRIDobj('resources/Clip_30mProject.tif');
% TODO: remove one
drainage_area = GRIDobj('resources/drainage_area_mdf.tif');
drainage_density = GRIDobj('resources/drainage_density_fixed.tif');
slope = GRIDobj('resources/slope_gauss.tif');
wetness_index = GRIDobj('resources/wetness_index.tif');
% TODO: load only the 6 selected sites
% 'locations' uses the same syntax as in site_selection
sites = load('site_selection.mat', 'locations');

%% c) Selected sites over a DEM
% TODO: merge points after loading the 6 selected
figure
imagesc(DEM), colorbar
hold on
plot(sites.locations{1}{1}, sites.locations{1}{2}, 'rx', 'MarkerSize', 8)
title('Site locations over the DEM')
shg

%% d) Selected sites over each topographic parameter
% TODO: exclude one
topo_params = {drainage_area, drainage_density, slope, wetness_index};
topo_title = {'Drainage area', 'Drainage density', 'Slope', 'Wetness index'};

for i = 1:3
    subplot(1, 3, i)
    imagesc(topo_params{i});
    hold on
    coord = sites.locations{i};
    plot(sites.locations{i}{1}, sites.locations{i}{2}, 'rx', 'MarkerSize', 8)
    title(topo_title{i});
    colorbar
    hold off
end