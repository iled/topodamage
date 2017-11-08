%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topographic Analysis - HW5
% Sam Mark, Arielle Woods, Julio Caineta
% Main file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generation of maps is performed within the script file generate_maps.m

% site selection is performed within the script file site_selection.m

%% load data
DEM = GRIDobj('resources/Clip_30mProject.tif');
% loading the averaged (mean filter) maps that are better for visualization
drainage_area = GRIDobj('resources/drainage_area_mdf_filtered.tif');
drainage_density = GRIDobj('resources/drainage_density_filtered.tif');
slope = GRIDobj('resources/slope_filtered.tif');
sites = load('coordinates.mat');

%% c) Selected sites over a DEM
figure
imagesc(DEM)
cz = colorbar;
cz.Label.String = 'Elevation [masl]';
colors = {'g', 'r', 'm'};
hold on
for i = 1:3
    plot(sites.coordinates{i}(:, 1), sites.coordinates{i}(:, 2), 'p', ...
        'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{i})
end
legend('Drainage area', 'Slope', 'Drainage density')
title('Site locations over the DEM')
shg

%% d) Selected sites over each topographic parameter
topo_params = {drainage_area, drainage_density, slope};
topo_title = {'Drainage area', 'Drainage density', 'Slope'};
bar_title = {'[m^2]', '[m^{-1}]', '[.]'};
colors = {'g', 'r', 'm'};

for i = 1:3
    subplot(1, 3, i)
    imagesc(topo_params{i});
    hold on
    plot(sites.coordinates{i}(:, 1), sites.coordinates{i}(:, 2), 'p', ...
        'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{i})
    title(topo_title{i});
    cz = colorbar;
    cz.Label.String = bar_title{i};
    hold off
end