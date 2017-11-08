function plot_site_selection(sub_images, locations, sub_titles, test_name, save)
% plot_site_selection(sub_images, locations, sub_titles, test_name, save)
% Plots the results of the given site selection test.
% sub_images is a GRIDobj or a cell array of GRIDobj's.
% locations is a cell array of locations, as given by one of the analysis
% functions.
% sub_titles is a cell array with the titles for each of the subplots.
% test_name (optional) is the name of the test, used to name the figure window
% and to save the figure in a file (optional). Default is an empty string.
% save (optional, {'on', 'off'}) indicates if the figure is to be saved in a png
% file. Default is 'off'.
% Julio Caineta Nov 7 2017

if nargin < 4
    test_name = '';
end
if nargin < 5
    save = 'off';
end
if strcmp(save, 'off')
    vis = 'on';
else
    vis = 'off';
end

f = figure('Visible', vis, 'NumberTitle', 'off', 'Name', test_name, ...
    'units','normalized','outerposition',[0 0 1 1]);

n_figs = numel(sub_titles);

for i = 1:n_figs
    if n_figs > 1
        subplot(2, 2, i)
    end
    % workaround to plot discrete scale bar for soils
    if numel(sub_images) == 1
        imagesc(sub_images)
        varname = 'sub_images.Z';
        name = sub_images.name;
    else
        imagesc(sub_images{i});
        varname = 'sub_images{i}.Z';
        name = sub_images{i}.name;
    end
    if sum(ismember('soil', name)) == 4
        colormap(jet(numel(unique(eval(varname)))))
    else
        colormap default
    end
    % end workaround for colorbar
    hold on
    coord = locations{i};
    gscatter(coord{1}, coord{2}, coord{5}, 'mgr', 'o', 8, 'off')
    gscatter(coord{3}, coord{4}, coord{6}, 'mgr', '+', 8, 'off')
    legend('low-low', 'high-low', 'low-med', 'high-med','low-high', 'high-high', ...
        'Location', 'northwestoutside')
    title(sub_titles{i});
    colorbar
    hold off
end

if strcmp(save, 'on')
    saveas(f, ['site_selection_' test_name], 'png');
end

end