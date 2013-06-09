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

% first we need to get all the coordinates in im1, so we can transform them
% into im2 coordinates.
[m1, n1, ~] = size(im1);
[X1,Y1] = meshgrid(1:m1,1:n1);
im1Coords = [ X1(:)'; Y1(:)'; ones(1,m1*n1) ]; % 3xMN matrix of all the coordinates in im1.

% transform into im2, then normalize the homogeneous part (so that w = 1)
im2Coords = H * im1Coords; % 3xMN matrix of coordinates in im2 reachable from im1.
im2Coords = hnormalise(im2Coords);

% use the transformed grid to transform the blue channel from im1
blue1 = im1(:,:,3);
blue2 = outIm(:,:,3);

for k = 1:(m1*n1)
    val = blue1(im1Coords(1,k), im1Coords(2,k));
    if (val == 0), continue; end;
    
    x2 = im2Coords(1,k);
    y2 = im2Coords(2,k);
    
    dx2 = x2 - floor(x2);
    dy2 = y2 - floor(y2);
    
    blue2(floor(x2),floor(y2)) = blue2(floor(x2),floor(y2)) + val * dx2 * dy2;
    blue2(ceil(x2), floor(y2)) = blue2(ceil(x2), floor(y2)) + val * (1-dx2) * dy2;
    blue2(floor(x2), ceil(y2)) = blue2(floor(x2), ceil(y2)) + val * dx2 * (1-dy2);
    blue2(ceil(x2),  ceil(y2)) = blue2(ceil(x2),  ceil(y2)) + val * (1-dx2) * (1-dy2);
end

%blue2 = blue2 - min(blue2(:));
%blue2 = 255 * blue2 / max(blue2(:));

outIm(:,:,3) = blue2;
outIm = uint8(outIm);

end

%% Helper functions

% Credit: The following functions were authored by
% % Peter Kovesi
% % http://www.csse.uwa.edu.au/~pk

% Normalize homogeneous coordinates so that w = 1
function nx = hnormalise(x)
    [rows,~] = size(x);
    nx = x;

    % Find the indices of the points that are not at infinity
    finiteind = find(abs(x(rows,:)) > eps);

    for r = 1:rows-1
        nx(r,finiteind) = x(r,finiteind)./x(rows,finiteind);
    end
    nx(rows,finiteind) = 1;
end