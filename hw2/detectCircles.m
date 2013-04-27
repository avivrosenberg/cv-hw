%% Detect Circles

function [ centers ] = detectCircles( im, radius, varargin )
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
parser.addRequired('im', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x));
parser.addRequired('radius', @isscalar);
parser.addOptional('usegradient', true, @islogical);
parser.addParamValue('bw_thresh',[.1, .2, .4, .8], @(x) isrow(x) || iscolumn(x));
parser.addParamValue('quantization',0.25,@isscalar);
parser.addParamValue('hough_thresh',0.75,@isscalar);

parser.parse(im, radius, varargin{:});

%% Convert image to grayscale.
im = rgb2gray(im);

%% Create a smooth BW image of the edges:
% * Convert image to BW with a variying threshold.
% * Get edges using 'canny' method with high thresholds (since it's a BW image).
% * Accumulate edges.
edges = false(size(im));

for bwThresh = parser.Results.bw_thresh
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
houghRadius = round(q * radius);

% Create index-matrices for all indexes in the hough plane.
[hX, hY] = meshgrid(1:size(houghPlane,1), 1:size(houghPlane,2));
hX = hX'; hY = hY';

% Find nonzero pixels in the image of edges.
[imX, imY] = find(edges);
% Move nonzero indices to the hough space
imX = imX .* q; imY = imY .* q;
n = length(imX); % number of nonzero pixels

ticId = tic; fprintf(1,' Progress =      ');

% For each nonzero pixel at [imX(i), imY(i)]:
% calculate a matrix of distances from that pixel to the rest of the
% plane. In the places where the distance equals houghRadius, increment
% the houghPlane.
for i=1:n
    fprintf(1,'\b\b\b\b\b%5.1f',(i/n) * 100);

    dist = round( sqrt( (hX - imX(i)).^2 + (hY - imY(i)).^2 ) );
    circle = dist == houghRadius;
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

fprintf(1,'\n Done. Elapsed: %.5f [sec]\n', toc(ticId));
%% DEBUG
figure; imshow(edges);
figure; imagesc(houghPlane);
figure; imshow(im); hold on; plot(cY, cX, 'r+', 'MarkerSize', 10); plot(cY, cX, 'ro', 'MarkerSize', 2*radius);
end

