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
parser.addParamValue('Iterations', 1000, @isscalar);
parser.addParamValue('MaxSquareDistThresh', 0.25, @isscalar);
parser.addParamValue('MinInliersPerTransform', 5, @isscalar);

parser.parse(matches1, matches2, varargin{:});
maxIter = parser.Results.Iterations;
distThresh = parser.Results.MaxSquareDistThresh;
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
useInliers = false;

% create homogeneous coordinates 
matches1 = [matches1; ones(1,nPoints)];
matches2 = [matches2; ones(1,nPoints)];

% shift and scale points so they are evenly distanced from origin
[matches1, T1] = normalise2dpts(matches1);
[matches2, T2] = normalise2dpts(matches2);
%% RANSAC
%
for i = 1:maxIter
    
    if (useInliers)
        % select inliers from previous iteration
        pts1 = matches1(:,bestInliers);
        pts2 = matches2(:,bestInliers);
    else
        % select random points
        randInd = randperm(nPoints);
        pts1 = matches1(:,randInd(1:nRand));
        pts2 = matches2(:,randInd(1:nRand));
    end
    
    % compute transformation matrix from pts1 and pts2
    currH = computeTransform(pts1, pts2);
    
    % check the transform
    [currDist, inliers] = calcTransformDist(matches1, matches2, currH, distThresh);
    numInliers = nnz(inliers);
    
    if (numInliers > minInliers && currDist < minDist)
        H = currH;
        minDist = currDist;
        bestInliers = inliers;
        % next iteration, use inliers from this iteration
        useInliers = true;
    else
        % next iteration, use random points
        useInliers = false;
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
