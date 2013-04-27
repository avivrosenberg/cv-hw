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
%Initialize Hough Plane, which is madeup of bins that accumulate the 'votes'.
houghPlane = zeros(ceil(q * imSize1), ceil(q * imSize2));
[houghSize1, houghSize2] = size(houghPlane);
houghRadius = round(q * radius);

% Create index-matrices for all indexes in the hough plane.
[hX, hY] = meshgrid(1:houghSize1, 1:houghSize2);
hX = hX'; hY = hY';

ticId = tic;

if (parser.Results.usegradient)
    [gradmag, graddir] = imgradient(edges);
    % Find indices in image corresponding to nonzero gradient
    ind = find(gradmag);
    [imX, imY] = ind2sub([imSize1, imSize2],ind); 
    
    n = length(ind); % number of nonzero pixels
    nzGradDir = graddir(ind) * pi / 180;
    
    % Find indices where votes should be cast (in the image plane).
    votesX1 = imX + houghRadius .* sin(nzGradDir);
    votesY1 = imY + houghRadius .* cos(nzGradDir);
    %votesX2 = imX - houghRadius .* sin(nzGradDir);
    %votesY2 = imY - houghRadius .* cos(nzGradDir);
    
    % Move the above indices to the hough plane.
    votesX1 = round(q.*votesX1); %votesX2 = round(q.*votesX2);
    votesY1 = round(q.*votesY1); %votesY2 = round(q.*votesY2);

    % Use logical indexing to filter out illegal indices.
    logicalInd1 = (votesX1 > 0 & votesX1 <= houghSize1) & (votesY1 > 0 & votesY1 <= houghSize2);
    %logicalInd2 = (votesX2 > 0 & votesX2 <= houghSize1) & (votesY2 > 0 & votesY2 <= houghSize2);
    votesX1 = votesX1(logicalInd1);
    votesY1 = votesY1(logicalInd1);
    %votesX2 = votesX2(logicalInd2);
    %votesY2 = votesY2(logicalInd2);
    
    % Convert the vote indices to houghPlane (quantize and round).
    houghInd1 = sub2ind(size(houghPlane), votesX1, votesY1);
    %houghInd2 = sub2ind(size(houghPlane), votesX2, votesY2);
    
    % Count number of votes for each index
    uniqHoughInd1 = unique(houghInd1);
    countHoughInd1 = hist(houghInd1, uniqHoughInd1)';
    %uniqHoughInd2 = unique(houghInd2);
    %countHoughInd2 = hist(houghInd2, uniqHoughInd2)';
    
    % Cast the votes in the hough plane...
    houghPlane(uniqHoughInd1) = houghPlane(uniqHoughInd1) + countHoughInd1;
    %houghPlane(uniqHoughInd2) = houghPlane(uniqHoughInd2) + countHoughInd2;
else
    % Find nonzero pixels in the image of edges.
    [imX, imY] = find(edges);
    % Move nonzero indices to the hough space
    imX = imX .* q; imY = imY .* q;
    n = length(imX); % number of nonzero pixels
    
    fprintf(1,' Progress =      ');
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

