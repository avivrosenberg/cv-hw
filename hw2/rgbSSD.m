function [ ssd ] = rgbSSD( im1, im2 )
%RBGSSD Compute Sum of Squared Difference between two RGB images.

%% Paramters Parsing
parser = inputParser;
parser.addRequired('im1', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x) && size(x,3) == 3);
parser.addRequired('im2', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x) && size(x,3) == 3);
parser.parse(im1, im2);
if (size(im1) ~= size(im2))
    error('Images sizes must match');
end

%% Calc SSD
diff = im1(:) - im2(:);
ssd = sum(diff.^2);

end

