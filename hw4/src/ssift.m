function [ keypoints, descriptors ] = ssift( im_color, varargin )
%SSIFT Summary of this function goes here
%   Detailed explanation goes here

%% Parameter Parsing
parser = inputParser;
parser.addRequired('im_color', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x) && ~isempty(x));
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
[rows, cols, ncolors] = size(im_color);

%% SIFT

if (ndims(im_color) == 3), im = rgb2gray(im_color);
else im = im_color; end

[ sift_kp, sift_desc ] = sift( im );

sift_kp_xy = round(sift_kp(1:2,:)) + 1; % Add one because sift keypoints are zero-based


%% Self Similarity

im = double(im_color);

% calculate the margin in which SSIM descriptors can't be calculated
margin = params.desc_rad + (params.patch_size - 1)./2;
xmargins = [margin + 1, cols - margin];
ymargins = [margin + 1, rows - margin];

% keypoints inside the margin need to be removed.
bad_kp_x_idx = (sift_kp_xy(1,:) < xmargins(1)) | (sift_kp_xy(1,:) > xmargins(2));
bad_kp_y_idx = (sift_kp_xy(2,:) < ymargins(1)) | (sift_kp_xy(2,:) > ymargins(2));
good_kp_idx = ~bad_kp_x_idx & ~bad_kp_y_idx;

sift_kp = sift_kp(:, good_kp_idx);
sift_desc = sift_desc(:, good_kp_idx);
sift_kp_xy = sift_kp_xy(:, good_kp_idx);

% iterate over keypoints, and calculate self-similarity
nkp = size(sift_kp, 2);
ss_desc = zeros(params.nrad * params.nang, nkp);
for k=1:nkp
    x = sift_kp_xy(1, k);
    y = sift_kp_xy(2, k);
    
    [ss_desc_xy, ~, ~, ~, ~] = mexCalcSsdescs(im, params,  [ x y ]);
    ss_desc(:, k) = ss_desc_xy;
end

% normalize
%ssim_desc = mat2gray(ssim_desc) * max(sift_desc(:));

%% 


keypoints = sift_kp;
%descriptors = [sift_desc; ss_desc];
descriptors = ss_desc;
%weights = [2*ones(1, size(sift_desc,1)), 1*ones(1, ncolors*nfilters)];
%[~, pca_desc, ~] = pca([sift_desc; filt_desc]', 'NumComponents', 200, 'VariableWeights', weights);
%descriptors = pca_desc';
end

