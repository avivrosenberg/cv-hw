function [ outIm ] = alignImages( im1, im2, H )
%ALIGNIMAGES Aligns two images, given the transformatino between them.
%   This functions will "blend" the images by simply taking the RG channels
%   from im2, and the B channel from im1.
%   H must be a transformation that transforms im1 into im2's coordinates,
%   so that H * im1 will be blended into im2.

%% Parameter parsing
%
parser = inputParser;
parser.addRequired('im1', @(x) ndims(x) == 3 && size(x,3) == 3);
parser.addRequired('im2', @(x) ndims(x) == 3 && size(x,3) == 3);
parser.addRequired('H',   @(x) all(size(x) == [3,3]));

%% Align
% construct the blue channel by first transforming im1 into im2's
% coordinate system with the given transformation.

% initialize output
outIm = zeros(size(im2));

% take Red and Green channels from im2.
outIm(:,:,1:2) = im2(:,:,1:2);

% take the blue channels
blue1 = im1(:,:,3);
blue2 = outIm(:,:,3);

% transform im1's blue channel with the transformation H
if (all(H(3,:) == [0,0,1]))
    tform = affine2d(H');
else
    tform = projective2d(H');
end
blue1t = imwarp(blue1, tform, 'cubic');

% calculate top-left coordinate of im1 trasformed to im2
x2y2 = hnormalise(H * [1;1;1]);
x2 = floor(x2y2(1));
y2 = floor(x2y2(2));

% place the transformed image into the correct im2 coordinates
[t1,t2] = size(blue1t);
x2vec = x2 + (1:t1);
y2vec = y2 + (1:t2);
blue2(x2vec, y2vec) = blue1t;

% clip to im2
blue2 = blue2(1:size(im2,1), 1:size(im2,2));

outIm(:,:,3) = blue2;
outIm = uint8(outIm);

end
