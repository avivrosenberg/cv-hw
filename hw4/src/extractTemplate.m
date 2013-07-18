function [ im_out ] = extractTemplate( im_tmpl, im_test, varargin )
%EXTRACTTEMPLATE Find a given template image in a test image
%   Detailed explanation goes here

%% Input parsing
%

featFunc = @(im) sift(im);

%% Features calculation
% 
% Extract frames (location, size, orientation of keypoints) and the
% descriptors for each keypoint.

[kp_tmpl, desc_tmpl] = featFunc( im_tmpl );
[kp_test, desc_test] = featFunc( im_test );

%% Features matching
%
% Match descriptors using Lowe's method.
% Use RANSAC to find an affine transform between the matching point, and
% discard all matches not fitting the transform.

matches = findMatches(desc_tmpl, desc_test, 'thresh', 0.8);

% Use voting for image center to determine which of the matches are on the object
[votes, vote_results] = keypointVotes(im_tmpl, kp_tmpl(:, matches(1,:)), kp_test(:, matches(2,:)), ...
                                      'clustertype','ward','maxclusters',7);

% Get new matches
matches = matches(:, vote_results);

% extract image coordinates of the matches.
pts1 = kp_tmpl(1:2,matches(1,:));
pts2 = kp_test (1:2,matches(2,:));

% find affine transform from matches 
[~, inliers, ~] = transformationMatches(pts1, pts2, 'TransformType', 'affine');

% use transform's inliers to filter out bad matches
matches = matches(:, inliers);

%% DEBUG

% Show votes
figure;
imshow(im_test); hold on; scatter(votes(1,vote_results),votes(2,vote_results),'r+');
hold off;

figure;
plotMatches(im_tmpl,im_test,kp_tmpl(1:2,:),kp_test(1:2,:),matches,'Placement','horz');
end

