function [ out_im, levels ] = kHSVquantizer(im, varargin)
%KHSVQUANTIZER Quantized image colors in the HSV plane with k-means

%% Paramters Parsing
parser = inputParser;
parser.addRequired('im', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x) && size(x,3) == 3);
parser.addParamValue('k',16, @isscalar);
parser.parse(im, varargin{:});

%% Convert to HSV
im_hsv = rgb2hsv(im);

%% Convert image to observations matrix
[imSz1,imSz2,imSz3] = size(im_hsv);
% Reshape HSV so the first column is Hue, second is Saturation, third is Value.
im_hsv = reshape(im_hsv, imSz1 * imSz2, imSz3);
% Take only first column, which contains the Hue values, as observations.
observations = double(im_hsv(:,1));

%% Use k-means
[idx, c] = kmeans(observations, parser.Results.k,...
                  'start','cluster',...
                  'emptyaction','drop');
              
% The rows of c now contain Hue values of cluster centers.
% We'll take the rows of c, in the order specified in idx, and put them
% into the first column of the hsv image (still reshaped).
im_hsv(:,1) = c(idx,:);

% Reshape it back and convert to RGB.
out_im = reshape(im_hsv, imSz1, imSz2, imSz3);
out_im = uint8(round(hsv2rgb(out_im) * 255));
levels = c(unique(idx));


end

