function [ im_out ] = extractTemplate( im_template, im_test, varargin )
%EXTRACTTEMPLATE Find a given template image in a test image
%   Detailed explanation goes here

%% Input parsing
%

featFunc = @(im) sift(im);

%% Features calculation
% 
% Extract frames (location, size, orientation of keypoints) and the
% descriptors for each keypoint.

[templFrames, templDescr] = featFunc( im_template );
[testFrames, testDescr] = featFunc( im_test );

%% Features matching
%
% Match descriptors using Lowe's method.
% Use RANSAC to find an affine transform between the matching point, and
% discard all matches not fitting the transform.

matches = findMatches(templDescr, testDescr, 'thresh', 0.8);

% extract image coordinates of the matches.
pts1 = templFrames(1:2,matches(1,:));
pts2 = testFrames (1:2,matches(2,:));

% find affine transform from matches 
[~, inliers, ~] = transformationMatches(pts1, pts2, 'TransformType', 'affine');
% use transform's inliers to filter out bad matches
matches = matches(:, inliers);

pts1 = templFrames(1:2,matches(1,:));
pts2 = testFrames (1:2,matches(2,:));

%% DEBUG
%
%

plotMatches(im_template,im_test,templFrames(1:2,:),testFrames(1:2,:),matches,'Placement','horz');
end

