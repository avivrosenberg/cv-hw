function [ tform, inliers, error ] = affineMatches( matches1, matches2, varargin )
%AFFINEMATCHES given matching points from two images, returns affine transform and inliers

[ tform, inliers, error ] = ...
    transformationMatches( matches1, matches2, 'TransformType' , 'affine', varargin{:});

end

