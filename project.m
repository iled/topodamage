%% load DEM

DEM = GRIDobj('resources/Clip_30mProject.tif');
DEM.cellsize=30;

hsize=5;
sigma=0.7;
h = fspecial('gaussian', hsize, sigma);

DEM2=DEM;
%filter the new DEM
DEM2.Z=imfilter(DEM.Z, h);
%% plot filtered DEM vs unfiltered DEM
imageschs(DEM2);
shg

subplot(1,2,1), imagesc(DEM), title('unfiltered DEM')
subplot(1,2,2), imagesc(DEM2), title('filtered DEM')

%%

%fill sinks in DEM
DEMf = fillsinks(DEM2);
%flow direction object using MDF
FD = FLOWobj(DEMf,'Dinf');
%flow accumulation
A = flowacc(FD)*DEM2.cellsize^2;
%slope
S=gradient8(DEM2);
%wetness index
W=DEMf;
W.Z=log((A.Z/DEMf.cellsize)./S.Z);

%%
subplot(1,3,1), imagesc(log(A)), title('log(A)')
subplot(1,3,2), imagesc(S), title('slope')

% TODO not sure if this is correct, the points don't plot to the edges of the plot (like they in the lecture slide image) 

subplot(1,3,3), loglog(A.Z(:), S.Z(:), '.'), axis square
xlim([900-100, max(A.Z(ix))])

%% drainage basins
% Find basins that drain to DEM boundaries, using SDF (D8)
flow_direction_sdf = FLOWobj(DEMf);
% assign each pixel to its basin, according to some integer ID
drainage_basins = drainagebasins(flow_direction_sdf);
flow_acc_sdf = flowacc(flow_direction_sdf) * DEM2.cellsize ^ 2;
% index of most commont basin (the largest)
basin_id = mode(drainage_basins.Z(:));

%% plot basins over the DEM
imageschs(DEM, shufflelabel(drainage_basins))


%% plot only the drainage basins
imagesc(drainage_basins)

%%
%drainage_basins.Z(drainage_basins.Z~=13)=0; %This makes a grid of 0's for all basins except to that whose pixels have a value of 13 (i.e., basin 13)
%imagesc(drainage_basins),shg

%basin 13 is fully contained within the DEM

ix=find(drainage_basins.Z==13);
loglog(A.Z(ix), S.Z(ix), '.'),
xlim([900-100, max(A.Z(ix))])

%%
nbins=40;
S.Z(S.Z==0)=NaN;%make zero slopes NaNs - this to enable processing it in log space binning
B=bin(log10(A.Z(ix)), log10(S.Z(ix)), nbins);
%now plot the data for the basin and teh logged data on top of it
loglog(A.Z(ix), S.Z(ix), '.'),shg
hold on
loglog(10.^B(:,1), 10.^B(:,2), 'ok', 'markerfacecolor', 'y'),shg
xlim([900-100, max(A.Z(ix))]);
xlabel ('A [m^2]');
ylabel ('S []');

% our threshold is around 3000 m^2 based on binning the data
%% approximate the drainage density over the map
% using our threshold
area_threshold = 3000 * 10;
area_map = prod(DEM2.size) * DEM2.cellsize ^ 2;
channel_px = find(flow_acc_sdf > area_threshold);
drainage_density_approx = length(channel_px) * DEM2.cellsize ./ area_map;

%% plot approximate drainage density
% area threshold per unit area
area_threshold_unit = area_threshold / DEM2.cellsize ^ 2;
% create stream object defining the upslope area threshold for channel initiaiton (minarea)
stream = STREAMobj(FD, 'minarea', area_threshold_unit);
imageschs(DEM2);
hold on
plot(stream, 'k')

%% compute drainage densit
drainage_density = drainagedensity(stream, FD);
%% plot drainage density
imagesc(log(drainage_density * DEM2.cellsize)) % TODO: confirm that it should be expontentially distributed
colorbar