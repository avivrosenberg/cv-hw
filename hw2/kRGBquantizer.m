function [ out_im ] = kRGBquantizer( im , varargin)
%KRGBQUANTIZER Quantized image colors with k-means

%% Paramters Parsing
parser = inputParser;
parser.addRequired('im', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x) && size(x,3) == 3);
parser.addParamValue('k',16, @isscalar);
parser.parse(im, varargin{:});

%% Convert image to observations matrix
[imSz1,imSz2,imSz3] = size(im);
observations = reshape(im, imSz1 * imSz2, imSz3);
observations = double(observations);
%% Use k-means
[idx, c] = kmeans(observations, parser.Results.k,...
                  'start','cluster',...
                  'emptyaction','drop');
              
% The rows of c now contain RBG values of cluster centers.
% We'll round them to get RBG values.
c = uint8(round(c));

% We'll take the rows of c, in the order specified in idx, and put them
% into the output image rows (the output image is shaped as observations).
out_im = uint8(zeros(size(observations)));
out_im(:,:) = c(idx,:);

out_im = reshape(out_im, imSz1, imSz2, imSz3);
end

