%% Detect Circles

function [ centers ] = detectCircles( im, radius, usegradient )
%%
%  Returns locations of circles with a given radius in an image.
%  Parameters:
% 
% * im - input image
% * radius - circle radius to look for
% * usegradient - whether or nor to use gradient direction on circle edges.
% 

%% Paramters Parsing
p = inputParser;
addRequired(p,'im', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x));
addRequired(p,'radius', @isscalar);
addOptional(p,'usegradient', true, @islogical);
parse(p, im, radius, usegradient);

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

q = 0.25; % Quantization factor
%Initialize Hough Plane, which is madeup of bins that accumulate the 'votes'.
houghPlane = zeros(q * size(edges,1), q * size(edges,2));
% Create index-matrices for all indexes in the hough plane.
[hX, hY] = meshgrid(1:size(houghPlane,1), 1:size(houghPlane,2));
% Find nonzero pixels in the image.
[imX, imY] = find(edges);
n = length(imX);

fprintf(1,'\n Progress =      ');

for i=1:n
    fprintf(1,'\b\b\b\b\b%5.1f',(i/n) * 100);
    
    dist = round( sqrt( (hX - q * imX(i)).^2 + (hY - q * imY(i)).^2 ) );
    circle = dist == radius;
    houghPlane = houghPlane + circle;
end

%% DEBUG
figure; imshow(edges);
figure; imagesc(houghPlane);
end

