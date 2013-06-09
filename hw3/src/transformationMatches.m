function [tform, inliers, error] = transformationMatches(matches1, matches2, varargin)
%TRANSFORMATIONMATCHES Summary of this function goes here
%   Detailed explanation goes here

%% Parameter parsing
%
parser = inputParser;
parser.addRequired('matches1', @(x) ndims(x) == 2 && size(x,1) == 2);
parser.addRequired('matches2', @(x) ndims(x) == 2 && size(x,1) == 2);
parser.addParamValue('TransformType', 'affine', @(x) strcmp(x,'affine') || strcmp(x,'homography'));
parser.addParamValue('MaxIterations', 1000, @isscalar);
parser.addParamValue('MinSquareDistThresh', 0.7, @isscalar);
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

% normalize points
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
error = minDist;

% denormalize
tform = T2\H*T1;

end

%% Helper functions
%
function [error, inliers] = calcTransformDist(x1, x2, H, th)
        
    % transform x2 with H and check distance from x1
    diff2 = (x1 - H * x2).^2;
    error = sqrt(sum(diff2,1));
    
    inliers = error < th;
    error = sum(error);
end