function [ centers ] = detectCircles( im, radius, usegradient )
%DETECTCIRCLES Returns locations of circles with a given radius in an image
%   Detailed explanation goes here

im = rgb2gray(im);
diskFilter3 = fspecial('disk', 3);

edges = false(size(im));
for bwThresh = [.1, .2, .4, .8]
    bwCurrent = im2bw(im, bwThresh);
    bwCurrent = imfilter(bwCurrent, diskFilter3);
    
    edgesCurrent = edge(bwCurrent, 'canny', [.8, .9]);
    edges = edges | edgesCurrent;
end

imshow(edges);
figure; imshow(edge(edges));
end
