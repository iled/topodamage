function s = bin(x,y,nbin,overlap)

% s = bin(x,y,numwin)
%
% Bins the function y(x) using nbin bins of equal x-width (as opposed to 
% bins containing a constant number of data points, as is the case with 
% many other routines). modified by es to exclude nan
%
% Input arguments:
%    x        - a vector containing the values of the independent variable
%    y        - a vector containing the values of the dependent variable y
%               at the values in x. length(y) == length(x).
%    nbin     - number of bins of equal x-width
%    overlap  - a flag that determines whether the bins will be 
%               non-overlapping (overlap=0, default) or 50% overlap 
%               (overlap=1)
%
% Output arguments:
%    s        - matrix consisting of nbin rows and 7 columns: 
%               [x (center of bin), mean, standard deviation, standard 
%               error, n, max, min]

% Copyright 2005-2008 Taylor Perron

if nargin<4, overlap=0; end

% if x and y are row vectors, make them column vectors
x = x(:);
y = y(:);


%es avoid nans
x1=x(~isnan(x) & ~isnan(y));
y1=y(~isnan(x) & ~isnan(y));
x=x1;
y=y1;

% sort x and y by x
sorted = sortrows([x y],1);
x = sorted(:,1); y = sorted(:,2);

% find the extrema of x
xmin = x(1); xmax = x(end);
xrange = xmax - xmin;



% determine the window width
if overlap,
    w = 2*xrange/(nbin+1); % for windows with 50% overlap (smoother)
else
    w = xrange/nbin; % for windows with no overlap
end

% Allocate memory for the binned data
s = zeros(nbin,7);

% loop through the bins
for i=1:nbin
	
    % determine min and max x values of current window position
    if overlap
        xlo = xmin+(i-1)*(w/2); % for windows with 50% overlap
    else
        xlo = xmin+(i-1)*w; % for windows with no overlap
    end
    xhi = xlo+w;
	
    % find min and max indices of x vector corresponding to this range
	window = find((x >= xlo) & (x <= xhi));
    mini = min(window); maxi = max(window);
    
    % calculate mean, standard dev, standard error, and n of points that 
    % fall within this window, but watch out for windows with only one 
    % point:
    if isempty(window)
        %s(i,:) = [mean([xlo xhi]) 0 0 0 0 0 0];
        s(i,:) = [nanmean([xlo xhi]) NaN 0 0 0 0 0];
    else
        s(i,:) = [nanmean([xlo xhi]) nanmean(y(mini:maxi)) ...
                 nanstd(y(mini:maxi)) nanstd(y(mini:maxi))/sqrt(maxi-mini+1)...
                 maxi-mini+1 max(y(mini:maxi)) min(y(mini:maxi))];
    end
    
end
