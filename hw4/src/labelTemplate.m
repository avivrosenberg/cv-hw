function [ labels ] = labelTemplate( tmpl_im, im, tform, seeds )
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

tmpl_im = rgb2gray( im2double(tmpl_im) );
im = rgb2gray( im2double(im) );


%% Set up the cost matrix for each edge (Boundary properties/Pairwise potentials)
%

[vC, hC] = pairwisePotential(im);
maxPairwise = 1   + max(max(vC + hC));

%% Set up the cost matrix for each pixel (Region properties/Unary potentials)
% We'll define label '1' as the object, and label '2' as the background.
% The cost should be a propability of the pixel belonging to the object or
% background.

% warp template into image coordinates
tmpl_warped = imwarp(tmpl_im, tform);
[x0_im, y0_im] = tform.transformPointsForward(1,1);
x_im = max( ceil(x0_im),1 ) : min( (x0_im+size(tmpl_warped,2)), size(im,2) );
y_im = max( ceil(y0_im),1 ) : min( (y0_im+size(tmpl_warped,1)), size(im,1) );
% get the region in the image that should contain the object
im_obj_region = im(y_im,x_im,:);

% crop the warped template to the same size as the image region
% do it such that it's center will remain in the center
yxdiff = size(tmpl_warped) - size(im_obj_region);
yxdiff = [floor(yxdiff./2); ceil(yxdiff./2)];
tmpl_warped = tmpl_warped((1+yxdiff(1,1)):(end-yxdiff(2,1)), (1+yxdiff(1,2)):(end-yxdiff(2,2)));

% Compute gradients to get more 'important' regions
[gm1,gd1] = imgradient(tmpl_warped);
[gm2,gd2] = imgradient(im_obj_region);

% Spead out the energy since the images are probably not exactly aligned
gfilt = fspecial('gauss', [13 13], sqrt(13));
gm1 = imfilter(gm1, gfilt);
gm2 = imfilter(gm2, gfilt);

% Compute diff
diff = sqrt((gm1 - gm2).^2 + (gd1 - gd2).^2);
diff = mat2gray(diff);

% ignore places where the gradients were low, they give us little info even
% if the diff was small.
diff(gm1<.01 | gm2<.01) = 1;

% propability that a pixel in the 'unknown' region actually does belong to
% the object.
pr = ones(size(seeds));
pr(y_im, x_im) = diff;

% High pr means high probablity of being an object pixel.
pr = 1 - pr;
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

smoothness_cost = 0.22 * (ones(2)-eye(2));

tic;[gch] = GraphCut('open', data_cost, smoothness_cost, vC, hC);toc;
tic;[gch] = GraphCut('set', gch, init_labels);toc;
tic;[gch, labels] = GraphCut('expand', gch, 200);toc;
GraphCut('close', gch);

%% Close holes to smooth out the label image

labels = 1-labels;
labels = imclose(labels, strel('disk', 5, 4));
labels = imfill(labels,'holes');

labels = 1-labels;
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

function [vC, hC] = pairwisePotential(im)   
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
