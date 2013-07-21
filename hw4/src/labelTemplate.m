function [ labels ] = labelTemplate( tmpl_im, im, seeds )
%LABELTEMPLATE Apply GraphCut segmentation to an image mask.
%   Seeds should be a matrix containing a mask for pixels of some image.
%   Value in 'seeds' should be as follows:
%       o '0' in pixels considered background;
%       o '1' in pixels which might belong to the object
%       o '2' in pixels which are known to belong to the object
%
%   The returned value, 'labels' is a logical mask the same size as 'seeds'.
%   it contains '1' at pixels which have been found to belong to the
%   object, and '0' elsewhere.

[nrow, ncol] = size(seeds);

%% Preprocessing
%

tmpl_im = im2double(tmpl_im);
im = im2double(im);

%seeds = imdilate(seeds,strel('disk',2,4));

% cluster the template image colors into k regions
k = 5;
data = reshape(tmpl_im, size(tmpl_im,1)*size(tmpl_im,2),3);
[~, means] = kmeans(data, k, 'maxiter',200, 'emptyaction', 'drop');

[ii,jj] = find( seeds == 1 );


%% Set up the cost matrix for each edge (Boundary properties/Pairwise potentials)
%

%[vC, hC] = pairwisePotential(im);
[vC, hC] = spatialCues(im);

maxPairwise = 1   + max(max(vC + hC));

%% Set up the cost matrix for each pixel (Region properties/Unary potentials)
% We'll define label '1' as the object, and label '2' as the background.
% The cost should be a propability of the pixel belonging to the object or
% background.

% propability that a pixel in the 'unknown' region actually does belong to
% the object.
pr = zeros(size(seeds));
for kk=1:length(ii)
    rgb = reshape( im(ii(kk),jj(kk),:), 1, 3);
    rgb = repmat(rgb, k, 1);
    diff2 = sum((means - rgb).^2,2);
    
    pr(ii(kk), jj(kk)) = min(diff2);
end
% Normalize ot [0,1] and inverse, so high pr means high probablity of being
% an object pixel.
pr = (pr-min(pr(:))) ./ (max(pr(:)-min(pr(:))));
pr(ii,jj) = 1 - pr(ii,jj);
% Add a small value to pr to prevent having zero probability
pr = pr + exp(-maxPairwise);
pr(pr >= 1) = 1 - exp(-maxPairwise);

lambda = 1.0;

% Costs for labling these pixels as 'object':
% should be low for actual object pixels
obj_cost = zeros(nrow, ncol);
obj_cost(seeds == 2) = 0; % known object pixels - zero cost
obj_cost(seeds == 0) = maxPairwise; % known background pixels - high cost
obj_cost(seeds == 1) = -log(pr(seeds==1)) * lambda; % unknown, use probablity. Low pr should be high cost.

% Costs for labling these pixels as 'background':
% should be high for actual object pixels
bkg_cost = zeros(nrow, ncol);
bkg_cost(seeds == 2) = maxPairwise; % known pbject pixels - high cost
bkg_cost(seeds == 0) = 0; % known background pixels - zero cost
bkg_cost(seeds == 1) = -log(1 - pr(seeds==1)) * lambda; % unknown, use probablity. High pr should be high cost.


% data_cost(r,c,l) equals the cost for assigning label l to pixel at (r,c).
data_cost(:,:,1) = obj_cost;
data_cost(:,:,2) = bkg_cost;

%% Initial labels
% label all object seeds with '0' (the object's label).

init_labels = ones(nrow, ncol);
init_labels(seeds == 2) = 0;

%% Applying GC

smoothness_cost = 2.0 * (ones(2)-eye(2));

tic;[gch] = GraphCut('open', data_cost, smoothness_cost, vC, hC);toc;
%tic;[gch] = GraphCut('set', gch, init_labels);toc;
tic;[gch, labels] = GraphCut('swap', gch);toc;
GraphCut('close', gch);
figure;imagesc(labels);

end

function p = unaryPotential(seeds, type)
if (strcmpi(type,'obj'))
    pr = 1 * (seeds == 2) + 0.9 * (seeds == 1) + 1e-9 * (seeds == 0);
else
    pr = 1 * (seeds == 0) + 0.1 * (seeds == 1) + 1e-9 * (seeds == 2);
end
p = -log(pr);
end

function [vert, horz] = pairwisePotential(im)
    
    [n, m] = size(im);
    
    vert_diff = im - [im(2:end,:); zeros(1, m)];
    horz_diff = im - [im(:,2:end)  zeros(n, 1)];
    
    vert = exp(- vert_diff.^2 ./ 2 ./ 1.0^2 );
    horz = exp(- horz_diff.^2 ./ 2 ./ 1.0^2 );
end

function [vC, hC] = spatialCues(im)
g = fspecial('gauss', [13 13], sqrt(13));
dy = fspecial('sobel');
vf = conv2(g, dy, 'valid');
sz = size(im);

vC = zeros(sz(1:2));
hC = vC;

for b=1:size(im,3)
    vC = max(vC, abs(imfilter(im(:,:,b), vf, 'symmetric')));
    hC = max(hC, abs(imfilter(im(:,:,b), vf', 'symmetric')));
end

vC = exp(-vC ./ 2 ./ 3^2);
hC = exp(-hC ./ 2 ./ 3^2);
end
