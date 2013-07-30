function [ template_points_in_test, tform, varargout] = findTemplate( im_tmpl, im_test, varargin )
%FINDTEMPLATE Find a given template image in a test image.
%   Given a template image (IM_TMPL), and a test image (IM_TEST) the
%   function tries to find the template in the test.
%
%   The template image is assumed to contain a single object, and the test
%   image could be anything.
%
%   The function returns coordinates of points in the test image that have
%   been found be be in an instance of the object in the test image.
%   Also, a transformation matrix (TFORM) is returns such that IM_TMP *
%   TFORM should give the entire object's location in IM_TEST.

%% Paramters Parsing

if (nargout ~= 0 && nargout ~= 2 && nargout ~= 4)
    error('Wrong number of output parameters: Should be 0, 2 or 4.');
end

parser = inputParser;
parser.addRequired('im_tmpl', @(x) ~isempty(x) );
parser.addRequired('im_test', @(x) ~isempty(x) );
parser.addOptional('kp_tmpl', [], @(x) ismatrix(x) && (isempty(x) || (ndims(x) == 2 && size(x,1) == 4)));
parser.addOptional('desc_tmpl', [], @(x) ismatrix(x) && (isempty(x) || (ndims(x) == 2)));
parser.addParamValue('clustertype', 'kmeans', @(x) strcmp(x,'ward') || strcmp(x,'kmeans'));
parser.addParamValue('maxclusters', 7, @isscalar);
parser.addParamValue('matchthresh', 0.9, @isscalar);
parser.addParamValue('transformtype', 'affine', @(x) strcmp(x,'affine') || strcmp(x,'homography'));
parser.addParamValue('descriptor', 'sift', @(x) strcmp(x,'sift') || strcmp(x,'ssift'));

parser.parse(im_tmpl, im_test, varargin{:});

if (ndims(im_tmpl) ~= ndims(im_test))
    error('Template and Image should have same dementionality.');
end

calc_tmpl_features = false;
if (isempty(parser.Results.kp_tmpl) || isempty(parser.Results.desc_tmpl))
    calc_tmpl_features = true;
else if(size(parser.Results.kp_tmpl,2) ~= size(parser.Results.desc_tmpl,2))
        error('Number of template keypoints and descriptors supplied does not match');
     end
end

if (strcmpi(parser.Results.transformtype, 'affine'))
    fitFunc = @(pts1,pts2) transformationMatches(pts1, pts2, 'TransformType', 'affine');
else
    fitFunc = @(pts1,pts2) ransacfithomography(pts1, pts2, .05);
end

if (strcmpi(parser.Results.descriptor, 'sift'))
    featFunc = @(im) sift(im);
    if (ndims(im_tmpl) == 3 && ndims(im_test) == 3)
        im_tmpl = rgb2gray(im_tmpl);
        im_test = rgb2gray(im_test);
    end
else
    featFunc = @(im) ssift(im);
end


%% Features calculation
% 
% Extract frames (location, size, orientation of keypoints) and the
% descriptors for each keypoint.
[kp_test, desc_test] = featFunc( im_test );

% Might not need to do it for the template image.
if (calc_tmpl_features)
    [kp_tmpl, desc_tmpl] = featFunc( im_tmpl );
else
    kp_tmpl = parser.Results.kp_tmpl;
    desc_tmpl = parser.Results.desc_tmpl;
end


%% Descripor matching
%
% Match descriptors using Lowe's method.
% Use RANSAC to find an affine transform between the matching point, and
% discard all matches not fitting the transform.
% matches will contain indices for kp_tmpl in row 1, and indices for kp_test in row2.
matches = findMatches(desc_tmpl, desc_test, 'thresh', parser.Results.matchthresh);

% Extract matching keypoints from both images (including scale and
% orientation data).
matching_kp_tmpl = kp_tmpl(:, matches(1,:));
matching_kp_test = kp_test(:, matches(2,:));

%% Voting for template object's center
%
% Use voting for image center to determine which of the matches are on the object.
% votes contains the calculated center location of the object in the test
% image, as calculated from each keypoint in matching_kp_test.
[votes, vote_results] = keypointVotes(im_tmpl, ...
                                      matching_kp_tmpl, matching_kp_test, ...
                                      'clustertype',parser.Results.clustertype, ...
                                      'maxclusters',parser.Results.maxclusters);
                                  
% Take only the votes (object center locations) which are in the largest
% cluster of center locations.
majority_votes = votes(:, vote_results);
                                  
% Get new matches: Take only matching keypoints which voted with the
% majority on the location of the object.
matches = matches(:, vote_results);

%% Filter matching keypoint based on image transform and it's inliers
%
% Extract image coordinates of the new matching keypoints.
pts1 = kp_tmpl(1:2,matches(1,:));
pts2 = kp_test (1:2,matches(2,:));

% Find image transform from keypoint locations in both images.
[tform, inliers] = fitFunc(pts1, pts2);
if (strcmpi(parser.Results.transformtype, 'affine'))
    tform = affine2d(tform');
else
    tform = projective2d(tform');
end
    
% use transform's inliers to filter matches even further.
% matches will now contain the keypoint indicies of keypoints which votes
% with the majority on the object center AND ALSO conform with the affine transform.
if (nnz(inliers) > 0)
    matches = matches(:, inliers);
    
    % extract keypoint locations in the test image using the new matches.
    % these keypoint location are assumed to all be inside the desired object,
    % and can therefor be use for segmentation of the object.
    template_points_in_test = kp_test (1:2,matches(2,:));
    template_points_in_test = unique(template_points_in_test','rows')';
else
    % images can't be matched...
    template_points_in_test = [];
    tform = [];
end


%% DEBUG

% only plot debug info if no outputs
if (nargout == 4)
    varargout{1} = kp_tmpl;
    varargout{2} = desc_tmpl;
    return;
else
    if (nargout > 0)
        return;
    end
end

% Show votes
figure('Name','Keypoint Votes: Computed center locations');
imshow(im_test); hold on;
scatter(majority_votes(1,:),majority_votes(2,:),'r+');
scatter(votes(1,~vote_results),votes(2,~vote_results),'yx');
scatter(majority_votes(1,inliers),majority_votes(2,inliers),'go');
scatter(template_points_in_test(1,:),template_points_in_test(2,:),'mp');
legend('Majority votes', 'Non-Majority Votes',['Majority votes of ' parser.Results.transformtype ' keypoints'], 'Actual Keypoints');
hold off;

figure;
plotmatches(im_tmpl, im_test, kp_tmpl(1:2,:), kp_test(1:2,:), matches);
end

