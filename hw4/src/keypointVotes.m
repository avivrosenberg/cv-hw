function [ votes, voting_result ] = keypointVotes( im_tmpl, kp_tmpl, kp_test, varargin)
%KEYPOINTVOTES Determines which keypoints from the test image belong to the
%object, by using a voting system.
%   im_tmpl - teplate image showing the desired object
%   fr_tmpl - 4xK matrix of keypoints in the template image
%   im_test - test image we would like to find the object in
%   fr_test - 4xK matrix of keyopints in the test image
%
%   The given keypoints are assumed to be MATCHING.
%
%   The function returns:
%       votes - 2xK matrix of the x,y location each keypoint in the test
%       image voted for as the test image's center.
%
%   This function is based on:
%       A. Suga, K. Fukuda, T. Takiguchi, Y. Ariki
%       "Object Recognition and Segmentation Using SIFT and Graph Cuts"
%       2008

%% Paramters Parsing
parser = inputParser;
parser.addRequired('im_tmpl', @ismatrix);
parser.addRequired('kp_tmpl', @ismatrix);
parser.addRequired('kp_test', @ismatrix);
parser.addParamValue('clustertype', 'ward', @(x) strcmp(x,'ward') || strcmp(x,'kmeans'));
parser.addParamValue('maxclusters', 5, @isscalar);

parser.parse(im_tmpl, kp_tmpl, kp_test, varargin{:});

if (size(kp_tmpl) ~= size(kp_test))
    error('Keypoint matrices must have same dimentions');
end

% if 10% of the keypoints is less than maxclusters, we can't initialize
% kmeans with 'cluster' mode.
if (size(kp_tmpl,2) * 0.1 <= parser.Results.maxclusters)
    kmeansinit = 'sample'; else kmeansinit = 'cluster';
end
%% Extract keypoint data

% frame x/y values are ZERO based, so add 1
x_tmpl = kp_tmpl(1,:) + 1;
y_tmpl = kp_tmpl(2,:) + 1;
x_test = kp_test(1,:) + 1;
y_test = kp_test(2,:) + 1;

scale_tmpl = kp_tmpl(3,:);
scale_test = kp_test(3,:);

theta_tmpl = kp_tmpl(4,:);
theta_test = kp_test(4,:);

%% Template image center

% in images, y corresonds to ROW and x to COL.
[n1, m1] = size(im_tmpl);
yC = floor(n1 * 0.5);
xC = floor(m1 * 0.5);

% Offsets of template keypoints from center
dx = x_tmpl - xC;
dy = y_tmpl - yC;

%% Votes Computation

% absolute angle of each template keypoint LOCATION (not orientation) relative image center
theta = atan(dy ./ dx);

% absolute distance of each template keypoint LOCATION from the image center
r = sqrt(dx.^2 + dy.^2);

% ratios of scale between test and template image keypoints
scale_ratio = scale_test ./ scale_tmpl;

% delta of keypoint ORIENTATIONS between test and template
dtheta = theta_tmpl - theta_test;

votesX = x_test + scale_ratio .* r .* cos(theta + dtheta);
votesY = y_test + scale_ratio .* r .* sin(theta + dtheta);
votes = [votesX; votesY];

%% Cluster votes

switch (parser.Results.clustertype)
    case 'ward'
        clusters_idx = clusterdata(votes','linkage','ward', 'maxclust', parser.Results.maxclusters);
        
    case 'kmeans'
        clusters_idx = kmeans(votes', parser.Results.maxclusters, ...
                              'emptyaction','drop', 'start', kmeansinit, 'replicates', 2);
end

% actual voting: select the cluster with the most members
largest_cluster = mode(clusters_idx);

% the voting results is a logical indexing of the keypoints that belong to
% the cluster that received the most votes.
voting_result = clusters_idx == largest_cluster;
end

