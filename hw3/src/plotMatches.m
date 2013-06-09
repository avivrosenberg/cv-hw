function h = plotMatches(im1, im2, kp1, kp2, matches, varargin)
%% PLOTMATCHES  Plot keypoint matches
%   H=PLOTMATCHES(IM1,IM2,KP1,KP2,MATCHES) plots the two images IM1 and IM2
%   and lines connecting the frames (keypoints) KP1 and KP2 as specified
%   by MATCHES.
%
%   KP1 and KP2 specify keypoint coordinates, one keypoint per column.
%   Their size should therefore be 2xN, where the first row contains
%   x-values and the second contains y-values.
%
%   MATCHES specifies a set of matches, one per column. The two
%   elementes of each column are two indexes in the sets KP1 and KP2
%   respectively.
%
%   The images IM1 and IM2 might be either both grayscale or both color.
%
%   Optional parameter: 'Placement'.
%       Sets the placement of the images relative each other in the output.
%       Values can be either 'horz' (default) or 'vert'.
%
%   This function returns a handles to the line plot created.

%% Inner functions
%
    function out = imToScaledDouble(in)
        if (isa(in,'integer'))
            in = double(in);
        else
            if (~isa(in,'double')) 
                error('Unsupported image type');
            end
        end
        in = in - min(in(:));
        out = in / max(in(:));
    end
%% Paramters Parsing
%
parser = inputParser;
parser.addRequired('im1', @(x) ndims(x) == 2 || ndims(x) == 3);
parser.addRequired('im2', @(x) ndims(x) == 2 || ndims(x) == 3);
parser.addRequired('kp1', @(x) ndims(x) == 2 && size(x,1) == 2);
parser.addRequired('kp2', @(x) ndims(x) == 2 && size(x,1) == 2);
parser.addRequired('matches', @(x) ndims(x) == 2 && size(x,1) == 2);
parser.addOptional('Placement', 'horz', @(x) strcmp(x,'horz') || strcmp(x,'vert'));

[m1,n1,k1] = size(im1);
[m2,n2,k2] = size(im2);

if (k1 ~= k2)
    error('im1 and im2 must both be grayscale or both color.');
end

parser.parse(im1, im2, kp1, kp2, matches, varargin{:});
placement = parser.Results.Placement;


im1 = imToScaledDouble(im1);
im2 = imToScaledDouble(im2);

%% Calculate output dimentions
%

% 'horz' - im1 will be left of im2, so we need to sum the columns
if (strcmp(placement, 'horz'))
    n3 = n1 + n2;
    m3 = max(m1,m2);
    offsetRow = 0;
    offsetCol = n1;
else %'vert' - im1 will be above im2, so we need to sum the rows
    m3 = m1 + m2;
    n3 = max(n1,n2);
    offsetRow = m1;
    offsetCol = 0;
end

%% Create output image
%
outIm = zeros(m3,n3,k1);

% Copy im1 and im2 into ouput image
outIm(offsetRow + (1:m2), offsetCol + (1:n2), :) = im2;
outIm(1:m1,1:n1,:) = outIm(1:m1,1:n1, :) + im1;

% Average pixel values on image boundaries
outIm(1:min(m1,m2), 1:min(n1,n2), :) = 0.5 * outIm(1:min(m1,m2),1:min(n1,n2),:);

%% Calcualte matches
%
numMatches = size(matches, 2);
padding = NaN * ones(1, numMatches);
%padding = [];

% Take keypoint XY values from indices of 'kpX' specifeid in 'matches'.
% take all x-coordinates (row 1) from kp1 and kp2. Offset so that they are
% correctly placed in the output image.

x = [ kp1(1, matches(1,:));             % matches on im1
    kp2(1, matches(2,:)) + offsetCol; % matches on im2
    padding ];

% do the same for y-coordinates.
y = [ kp1(2, matches(1,:));             % matches on im1
    kp2(2, matches(2,:)) + offsetRow; % matches on im2
    padding ];
%% Draw image and match lines
%
axes('Position', [0 0 1 1]);
imagesc(outIm); colormap gray; hold on;
axis off; axis image;

h = line(x(:)', y(:)', 'Marker','+', 'MarkerSize',7, 'Color','g');
end

