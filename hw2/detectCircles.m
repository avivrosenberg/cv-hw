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

%% Convert image to grayscale.
im = rgb2gray(im);


%% Create a smooth BW image of the edges:
% * Convert image to BW with a variying threshold, and smooth it with a filter.
% * Get edges using 'canny' method with high thresholds (since it's a BW image)
% * Accumulate edges.
edges = false(size(im));
diskFilter3 = fspecial('disk', 3);
for bwThresh = [.1, .2, .4, .8]
    bwCurrent = im2bw(im, bwThresh);
    bwCurrent = imfilter(bwCurrent, diskFilter3);
    
    edgesCurrent = edge(bwCurrent, 'canny', [.8, .9]);
    edges = edges | edgesCurrent;
end

%% DEBUG
imshow(edges);
figure; imshow(edge(edges));
end

