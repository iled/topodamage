

DEM = GRIDobj('resources/Clip_30mProject.tif');
DEM.cellsize=30;

hsize=5;
sigma=0.7;
h = fspecial('gaussian', hsize, sigma);

DEM2=DEM;
%filter the new DEM
DEM2.Z=imfilter(DEM.Z, h);
imageschs(DEM2);
shg

subplot(1,2,1), imagesc(DEM), title('unfiltered DEM')
subplot(1,2,2), imagesc(DEM2), title('filtered DEM')

%%

%fill sinks in DEM
DEMf = fillsinks(DEM2);
%flow direction object:
FD = FLOWobj(DEMf,'Dinf');
%flow accumulation
A = flowacc(FD)*DEM2.cellsize^2;
%slope
S=gradient8(DEM2);
%wetness index
W=DEMf;
W.Z=log((A.Z/DEMf.cellsize)./S.Z);

% left off at doing code for drainage density


