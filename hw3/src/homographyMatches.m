function [ tform, inliers, error ] = homographyMatches( matches1, matches2, varargin )
%HOMOGRAPHYMATCHES given matching points from two images, returns affine
%transform and inliers.

[ tform, inliers, error ] = ...
    transformationMatches( matches1, matches2, 'TransformType' , 'homography', varargin{:});

end

