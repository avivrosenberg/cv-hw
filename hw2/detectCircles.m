%% Detect Circles

function [ centers ] = detectCircles( im, radius, usegradient, varargin )
%%
%  Returns locations of circles with a given radius in an image.
%  Parameters:
% 
% * im - input image
% * radius - circle radius to look for
% * usegradient - whether or nor to use gradient direction on circle edges.
% 

%% Paramters Parsing
parser = inputParser;
addRequired(parser,'im', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x));
addRequired(parser,'radius', @isscalar);
addOptional(parser,'usegradient', true, @islogical);
addOptional(parser,'quantization',0.25,@isscalar);
addOptional(parser,'hough_thresh',0.75,@isscalar);
parse(parser, im, radius, usegradient, varargin{:});

%% Convert image to grayscale.
im = rgb2gray(im);

%% Create a smooth BW image of the edges:
% * Convert image to BW with a variying threshold.
% * Get edges using 'canny' method with high thresholds (since it's a BW image).
% * Accumulate edges.
edges = false(size(im));

for bwThresh = [.1, .2, .4, .8]
    bwCurrent = im2bw(im, bwThresh); 
    edgesCurrent = edge(bwCurrent, 'canny', [.8, .9]);
    edges = edges | edgesCurrent;
end

%% Hough Transform
% * Every pixel in the image corresponds to a circle in the hough plane.
% Each pixel in a hough-plane circle is a 'vote' for that pixel as a
% circle-center in the original image.
% * The hough plane accumulates these votes into bins.
% * Bins with high values correspond to circle-centers in the image plane.

q = parser.Results.quantization; % Quantization factor
%Initialize Hough Plane, which is madeup of bins that accumulate the 'votes'.
houghPlane = zeros(ceil(q * size(edges,1)), ceil(q * size(edges,2)));
r = round(q * radius);

% Create index-matrices for all indexes in the hough plane.
[hX, hY] = meshgrid(1:size(houghPlane,1), 1:size(houghPlane,2));
hX = hX'; hY = hY';

% Find nonzero pixels in the image.
[imX, imY] = find(edges);
n = length(imX);

fprintf(1,'\n Progress =      ');
for i=1:n
    fprintf(1,'\b\b\b\b\b%5.1f',(i/n) * 100);
    
    dist = round( sqrt( (hX - q * imX(i)).^2 + (hY - q * imY(i)).^2 ) );
    circle = dist == r;
    houghPlane = houghPlane + circle;
end

%% Find Circles
% Take maximal points in the Hough plane, and compute their index in the
% image plane.
% These indices are where circles exists in the image.

t = parser.Results.hough_thresh; % Threshold percent for number of votes needed for a circle

[cX, cY] = find(houghPlane >= t * max(max(houghPlane)));
cX = cX ./ q; cY = cY ./ q;
centers(:,1) = cX; centers(:,2) = cY;
%% DEBUG
figure; imshow(edges);
figure; imagesc(houghPlane);
figure; imshow(im); hold on; plot(cY, cX, 'r+', 'MarkerSize', 10); plot(cY, cX, 'ro', 'MarkerSize', 2*radius);
end

