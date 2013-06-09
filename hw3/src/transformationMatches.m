function [tform, inliers, totalError] = transformationMatches(matches1, matches2, varargin)
%TRANSFORMATIONMATCHES Computes affine or homography transformation using RANSAC.
%   The function recieves two sets of coordinates in matches1 and matches2.
%   These sould be matrices sized 2xN.
%   
%   The function returns the best transformation found for which
%   tform * matches1 = matches2, with minimal square error (total error for
%   the transformation on all of the points is returned).
%
%   It also returns the column-indices of inliers (as a logical array), so
%   that the inliers are matches1(:,inliers) and matches2(:,inliers).

%% Parameter parsing
%
parser = inputParser;
parser.addRequired('matches1', @(x) ndims(x) == 2 && size(x,1) == 2);
parser.addRequired('matches2', @(x) ndims(x) == 2 && size(x,1) == 2);
parser.addParamValue('TransformType', 'affine', @(x) strcmp(x,'affine') || strcmp(x,'homography'));
parser.addParamValue('MaxIterations', 1000, @isscalar);
parser.addParamValue('MinSquareDistThresh', 0.25, @isscalar);
parser.addParamValue('MinInliersPerTransform', 5, @isscalar);

parser.parse(matches1, matches2, varargin{:});
maxIter = parser.Results.MaxIterations;
distThresh = parser.Results.MinSquareDistThresh;
minInliers = parser.Results.MinInliersPerTransform;

[~, nPoints] = size(matches1);
[~, n2] = size(matches2);

if (nPoints ~= n2)
    error('matches1 and matches2 must contain same number of elements.');
end

% Set number of random points to select each iteration.
if (strcmp(parser.Results.TransformType, 'affine'))
    nRand = 3;
    computeTransform = @affineTransform;
else
    nRand = 4;
    computeTransform = @homographyTransform;
end

%% Initialize
%
H = zeros(3,3);
minDist = Inf;
bestInliers = zeros(1,nPoints);

% create homogeneous coordinates 
matches1 = [matches1; ones(1,nPoints)];
matches2 = [matches2; ones(1,nPoints)];

% shift and scale points so they are evenly distanced from origin
[matches1, T1] = normalise2dpts(matches1);
[matches2, T2] = normalise2dpts(matches2);
%% RANSAC
%
for i = 1:maxIter
    % select random points
    randInd = randperm(nPoints);
    pts1 = matches1(:,randInd(1:nRand));
    pts2 = matches2(:,randInd(1:nRand));
    
    % compute transformation matrix
    currH = computeTransform(pts1, pts2);
    
    % check the transform
    [currDist, inliers] = calcTransformDist(matches1, matches2, currH, distThresh);
    numInliers = nnz(inliers);
    
    if (numInliers > minInliers && currDist < minDist)
        H = currH;
        minDist = currDist;
        bestInliers = inliers;
    end
end

% save outputs
inliers = bestInliers;
totalError = minDist;

% denormalize
tform = T2\H*T1;

end

%% Helper functions
%

% Error function for a a given transformation.
function [error, inliers] = calcTransformDist(x1, x2, H, th)
        
    % transform x1 with H (and normalize homogeneous coordinates)
    x2New = hnormalise( H * x1 );
    
    % check distance of new x2 from original x2
    diff2 = (x2New - x2).^2;
    error = sqrt(sum(diff2,1));
    
    inliers = error < th;
    error = sum(error);
end

% Credit: The following functions were authored by
% % Peter Kovesi
% % http://www.csse.uwa.edu.au/~pk

% Normalize homogeneous coordinates so that w = 1
function nx = hnormalise(x)
    [rows,~] = size(x);
    nx = x;

    % Find the indices of the points that are not at infinity
    finiteind = find(abs(x(rows,:)) > eps);

    for r = 1:rows-1
        nx(r,finiteind) = x(r,finiteind)./x(rows,finiteind);
    end
    nx(rows,finiteind) = 1;
end

% Shift and scale points so that they're centered at origin with a 
% mean distance of sqrt(2).
function [newpts, T] = normalise2dpts(pts)
    if size(pts,1) ~= 3
        error('pts must be 3xN');
    end

    % Find the indices of the points that are not at infinity
    finiteind = find(abs(pts(3,:)) > eps);
    
    % For the finite points ensure homogeneous coords have scale of 1
    pts(1,finiteind) = pts(1,finiteind)./pts(3,finiteind);
    pts(2,finiteind) = pts(2,finiteind)./pts(3,finiteind);
    pts(3,finiteind) = 1;
    
    c = mean(pts(1:2,finiteind)')';            % Centroid of finite points
    newp(1,finiteind) = pts(1,finiteind)-c(1); % Shift origin to centroid.
    newp(2,finiteind) = pts(2,finiteind)-c(2);
    
    dist = sqrt(newp(1,finiteind).^2 + newp(2,finiteind).^2);
    meandist = mean(dist(:));
    
    scale = sqrt(2)/meandist;
    
    T = [scale   0   -scale*c(1)
         0     scale -scale*c(2)
         0       0      1      ];
    
    newpts = T*pts;
end