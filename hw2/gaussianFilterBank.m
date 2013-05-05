function [ filterBank ] = gaussianFilterBank( angles, sizes, varargin )
%GAUSSIANFILTERBANK Create a filter bank of gaussian derivative filters of
%multiple orientations and scales.

%% Paramters Parsing
parser = inputParser;
parser.addRequired('angles', @(x) isrow(x) || iscolumn(x));
parser.addRequired('sizes', @(x) isrow(x) || iscolumn(x));
parser.addOptional('sigma', 1, @isscalar);
parser.parse(angles, sizes, varargin{:});

%% Create Filter Bank
filterBank = cell(length(sizes), length(angles));
bankRow = 0;

ddx = fspecial('sobel')';
ddy = fspecial('sobel');
sigma = parser.Results.sigma;

for sz=sizes
    bankRow = bankRow + 1; bankCol = 0;
    gx = conv2(fspecial('gaussian', sz, sigma), ddx, 'same');
    gy = conv2(fspecial('gaussian', sz, sigma), ddy, 'same');
    
    for ang=angles
        bankCol = bankCol + 1;
        theta = ang * pi / 180;
        filter = gx .* cos(theta) + gy .* sin(theta);
        
        filterBank{bankRow, bankCol} = filter;
    end
end

end
