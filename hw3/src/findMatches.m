function [ matches ] = findMatches( desc1, desc2, varargin )
%% MATCHES=findMatches(DESC1, DESC2) matches the two sets of SIFT
%    descriptors DESCR1 and DESCR2. DESCR1 and DESCR2 should be matrices with
%    the same number of rows, since each column in them is a descriptor.
%
%    The function uses the same algorithm suggested by D. Lowe to
%    reject matches that are too ambiguous.
%    
%    Output is a 2xK matrix, where K is the number of matches found.
%    The first row contains indices of matches from desc1, and the second
%    from desc2.
% 
%    findMatches(DESC1, DESC2, THRESH) accepts a specified
%    threshold THRESH. A descriptor D1 is matched to a descriptor D2
%    only if the distance d(D1,D2) multiplied by THRESH is not greather
%    than the distance of D1 to all other descriptors. The default
%    value of THRESH is 0.8.

%% Paramters Parsing
parser = inputParser;
parser.addRequired('desc1', @(x) ndims(x) == 2);
parser.addRequired('desc2', @(x) ndims(x) == 2);
parser.addParamValue('thresh', 0.8, @isscalar);

parser.parse(desc1, desc2, varargin{:});

[m1,n1] = size(desc1);
[m2,n2] = size(desc2);

if(m1 ~= m2)
  error('desc1 and desc2 must have same number of rows.') ;
end

thresh = parser.Results.thresh;
%% Compute distances for each descriptor pair
%  The element d(i,j) contains the distance between column i from descA
%  and column j from descB, where descA is the one with less columns.
%  This means that d should always have less rows than columns.

d = pdist2(desc1', desc2');
flipOutput = false;

% If desc2 has more columns than desc1, transpose d so it will have less
% rows than columns.
if (n1 > n2)
    d = d';
    flipOutput = true;
end

[m3,~] = size(d);

%% Find matches
%  Sort each row i of d to obtain the closest and second-closest distances
%  from descriptor i in desc1 to each descriptor j in desc2.

% prealocate matches with it's maximal possible size.
matches = zeros(2, m3);
currMatch = 0;

% Loop over rows of d (remember that m3 = min(m1,m2))
for i = 1:m3
    [~, ind] = sort(d(i,:));
    
    % compute the ratio closest/second_closest
    ratio = d(i,ind(1)) / d(i,ind(2));
    
    % discard matches where ratio is over the threshold
    if (ratio > thresh), continue; end
    
    % keep all other matches
    currMatch = currMatch + 1;
    matches(:,currMatch) = [i; ind(1)];
end

%% Prepare output

% Handle case of no matches found
if (currMatch == 0)
    matches = [];
    return;
end

% Take only the actual matches, discard empty columns.
matches = matches(:, 1:currMatch);

% Make sure that row 1 in matches corresponds to indices of desc1,
% and row 2 to indices of desc2.
if (flipOutput)
    tmp = matches(1,:);
    matches(1,:) = matches(2,:);
    matches(2,:) = tmp;
end

end
