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
houghRadius = q * radius;

ticId = tic;

if (parser.Results.usegradient)
    [gradmag, graddir] = imgradient(edges);
    % Remove small gradients
    gradmag(gradmag < 0.75 * max(max(gradmag))) = 0;
    % Find indices in image corresponding to nonzero gradient
    ind = find(gradmag);
    % Convert linear indices to subscripts.
    % switch y/x because in the image x is horizontal (like column).
    [imY, imX] = ind2sub([imSize1, imSize2],ind);
    
    nzGradDir = graddir(ind) * pi / 180;
    
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
    
    % Convert the vote indices to houghPlane (quantize and round).
    houghInd = sub2ind(size(houghPlane), votesY, votesX);
    
    % Count number of votes for each index
    uniqHoughInd = unique(houghInd);
    countHoughInd = hist(houghInd, uniqHoughInd)';
    
    % Cast the votes in the hough plane...
    houghPlane(uniqHoughInd) = houghPlane(uniqHoughInd) + countHoughInd;
else
    % Create index-matrices for all indexes in the hough plane.
    % Switch y/x because in the image x is horizontal (like column).
    [hY, hX] = meshgrid(1:houghSize1, 1:houghSize2);
    hX = hX'; hY = hY';
    
    % Find nonzero pixels in the image of edges.
    [imY, imX] = find(edges);
    % Move nonzero indices to the hough space
    imX = imX .* q; imY = imY .* q;
    n = length(imX); % number of nonzero pixels
    rhRadius = round(houghRadius);
    
    fprintf(1,' Progress =      ');
    % For each nonzero pixel at [imX(i), imY(i)]:
    % calculate a matrix of distances from that pixel to the rest of the
    % plane. In the places where the distance equals houghRadius, increment
    % the houghPlane.
    for i=1:n
        fprintf(1,'\b\b\b\b\b%5.1f',(i/n) * 100);
        
        dist = round( sqrt( (hX - imX(i)).^2 + (hY - imY(i)).^2 ) );
        circle = dist == rhRadius;
        houghPlane = houghPlane + circle;
    end
end

%% Find Circles
% Take maximal points in the Hough plane, and compute their index in the
% image plane.
% These indices are where circles exists in the image.

% We're looking for peaks in the hough plane, so we'll make them 'stand out'.
houghPlane = imfilter(houghPlane, fspecial('average',2));
houghPlane = houghPlane.^1.25;
%[houghPlane, ~] = imgradient(houghPlane);


t = parser.Results.hough_thresh; % Threshold percent for number of votes needed for a circle

[cY, cX] = find(houghPlane >= t * max(max(houghPlane)));
cX = cX ./ q; cY = cY ./ q;
centers(:,1) = cX; centers(:,2) = cY;

fprintf(1,'\n Done. Elapsed: %.5f [sec]\n', toc(ticId));
%% DEBUG
figure; imshow(edges);
figure; imagesc(houghPlane); title('hough plane');
if (parser.Results.usegradient)
    figure; imagesc(gradmag); title('gradmag');
    figure; imagesc(graddir); title('graddir');
end
figure; imshow(im); hold on; plot(cX, cY, 'r+', 'MarkerSize', 10); plot(cX, cY, 'ro', 'MarkerSize', 2*radius);
end

