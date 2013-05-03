%% Detect Circles

function [ centers ] = detectCirclesMulti( im, varargin )
%%
%  Returns locations of circles in an image.
%  Parameters:
%
% * im - input image
%

%% Paramters Parsing
parser = inputParser;
parser.addRequired('im', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x));
parser.addParamValue('bw_thresh',[.1, .2, .4, .8], @(x) isrow(x) || iscolumn(x));
parser.addParamValue('quantization',0.25,@isscalar);
parser.addParamValue('hough_thresh',0.6,@isscalar);

parser.parse(im,varargin{:});

%% Convert image to grayscale.
im = rgb2gray(im);
[imSize1, imSize2] = size(im);
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


ticId = tic;

% radii to search for
radii = 3:0.25*min(size(im));

%Initialize Hough Plane, which is madeup of bins that accumulate the 'votes'.
houghPlane = zeros(ceil(q * imSize1), ceil(q * imSize2), length(radii));
[houghSize1, houghSize2, houghSize3] = size(houghPlane);

[gradmag, graddir] = imgradient(edges);
% Remove small gradients
gradmag(gradmag < 0.75 * max(max(gradmag))) = 0;

% Find indices in image corresponding to nonzero gradient
ind = find(gradmag);
% Convert linear indices to subscripts.
% switch y/x because in the image x is horizontal (like column).
[imY, imX] = ind2sub([imSize1, imSize2],ind);

nzGradDir = graddir(ind) * pi / 180;

slice = 0;
for radius = radii
    slice = slice+1;
    % Find indices where votes should be cast (in the image plane).
    votesX = imX + radius .* cos(nzGradDir);
    votesY = imY + radius .* sin(nzGradDir);
    
    % Move the above indices to the hough plane.
    votesX = round(q.*votesX);
    votesY = round(q.*votesY);
    
    % Use logical indexing to filter out illegal indices.
    logicalInd = (votesX > 0 & votesX <= houghSize2) & (votesY > 0 & votesY <= houghSize1);
    votesX = votesX(logicalInd);
    votesY = votesY(logicalInd);
    
    % Count number of votes for each index
    houghInd = sub2ind([houghSize1, houghSize2], votesY, votesX);
    uniqHoughInd = unique(houghInd);
    countHoughInd = hist(houghInd, uniqHoughInd)';
    
    % Cast the votes in the hough plane...
    houghSlice = zeros(houghSize1, houghSize2);
    houghSlice(uniqHoughInd) = houghSlice(uniqHoughInd) + countHoughInd;
    houghPlane(:,:,slice) = houghSlice;
end
%% Find Circles
% Take maximal points in the Hough plane, and compute their index in the
% image plane.
% These indices are where circles exists in the image.

% We're looking for peaks in the hough plane, so we'll make them 'stand out'.
houghPlane = houghPlane.^2;

t = parser.Results.hough_thresh; % Threshold percent for number of votes needed for a circle

cInd = find(houghPlane >= t * max(max(max(houghPlane))));
[cY, cX, cR] = ind2sub(size(houghPlane), cInd);
cX = cX ./ q; cY = cY ./ q;

centers(:,1) = cX;
centers(:,2) = cY;
centers(:,3) = cR + min(radii);

fprintf(1,'\n Done. Elapsed: %.5f [sec]\n', toc(ticId));

%% Output
if (size(centers,1) > 0)
    figure; imshow(im); hold on; plot(centers(:,1), centers(:,2), 'r+', 'MarkerSize', 10);
    
    for i=1:size(centers,1)
        plot(centers(i,1), centers(i,2), 'ro', 'MarkerSize', 2*centers(i,3));
    end   
end

end

