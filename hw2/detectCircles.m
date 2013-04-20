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

%% DEBUG
figure; imshow(edges);
end

