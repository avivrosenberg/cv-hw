function [ keypoints, descriptors ] = ssift( im_color, varargin )
%SSIFT Self-Similarity SIFT descriptor
%   This function computes the SSIFT keypoints and descriptors for a given
%   image.
%   The SSIFT descriptors uses SIFT to find scale invariant keypoints, but
%   then calculates the Self-Similarity descriptor at these keypoints while
%   accounting for scale (the Self-Similarity descriptor itself does not
%   account for scale).

%% Parameter Parsing
parser = inputParser;
parser.addRequired('im_color', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x) && ~isempty(x));
parser.addParamValue('scaleclusters',5, @isscalar);
parser.addParamValue('patch_size',5, @isscalar);
parser.addParamValue('desc_rad',40, @isscalar);
parser.addParamValue('nrad',3, @isscalar);
parser.addParamValue('nang',12, @isscalar);
parser.addParamValue('var_noise',2700, @isscalar);
parser.addParamValue('saliency_thresh',1.0, @isscalar);
parser.addParamValue('homogeneity_thresh',1.0, @isscalar);
parser.addParamValue('snn_thresh',1.0, @isscalar);

parser.parse(im_color, varargin{:});

% create parameters struct for self similarity descriptor
params = struct;
params.patch_size = parser.Results.patch_size;
params.desc_rad = parser.Results.desc_rad;
params.nrad = parser.Results.nrad;
params.nang = parser.Results.nang;
params.var_noise = parser.Results.var_noise;
params.saliency_thresh = parser.Results.saliency_thresh;
params.homogeneity_thresh = parser.Results.homogeneity_thresh;
params.snn_thresh = parser.Results.snn_thresh;

%%
[rows, cols, ~] = size(im_color);

%% SIFT

if (ndims(im_color) == 3), im = rgb2gray(im_color);
else im = im_color; end

[ sift_kp, sift_desc ] = sift( im );

sift_kp_xy = round(sift_kp(1:2,:)) + 1; % Add one because sift keypoints are zero-based


%% Self Similarity

im = double(im_color);

% cluster sift keypoints according to their sigma (proportional to radius
% of region around the keypoint).
[sigma_idx, sigma] = kmeans(sift_kp(3,:)', parser.Results.scaleclusters, 'EmptyAction', 'drop');
sigma = sigma';
sigma_idx = sigma_idx';

% scale each keypoint xy coordinate and replace sigma with clustered values
% these will be used as the keypoints for self similarity
ss_kp = [
    round(sift_kp_xy(1,:) ./ sigma(sigma_idx));
    round(sift_kp_xy(2,:) ./ sigma(sigma_idx));
    sigma(sigma_idx);
    ];

% remove keypoints which are located in the margin of their scaled image
% (the self similarity descriptor cannot be calculated for these).
margin = params.desc_rad + (params.patch_size - 1)./2;

scaled_rows = floor(rows ./ ss_kp(3,:));
scaled_cols = floor(cols ./ ss_kp(3,:));

bad_kp_x_idx = (ss_kp(1,:) < margin + 1) | (ss_kp(1,:) > scaled_cols - margin);
bad_kp_y_idx = (ss_kp(2,:) < margin + 1) | (ss_kp(2,:) > scaled_rows - margin);
good_kp_idx = ~(bad_kp_x_idx | bad_kp_y_idx);

% throw away bad keypoints and their descriptors
sift_kp = sift_kp(:, good_kp_idx);
sift_desc = sift_desc(:, good_kp_idx);
ss_kp = ss_kp(:, good_kp_idx);
sigma_idx = sigma_idx(:, good_kp_idx);

% scale the input image K times
im_scaled = cell(1, length(sigma));
for i=1:length(sigma)
    if (isnan(sigma(i))), continue; end
    im_scaled{i} = imresize(im, 1/sigma(i));
end

% iterate over keypoints, and calculate self-similarity
nkp = size(sift_kp, 2);
ss_desc = zeros(params.nrad * params.nang, nkp);
for k=1:nkp
    x = ss_kp(1, k);
    y = ss_kp(2, k);
    
    [ss_desc_xy, ~, ~, ~, ~] = mexCalcSsdescs(im_scaled{sigma_idx(k)}, params,  [ x y ]);
    ss_desc(:, k) = ss_desc_xy;
end

%% Results

% take keypoint locations in the original scale
keypoints = sift_kp;

% normalize all descriptors to [0,1]
%descriptors = [mat2gray(sift_desc); mat2gray(ss_desc)];
descriptors = mat2gray(ss_desc);

end

