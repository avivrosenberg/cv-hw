function H = homographyTransform(x1, x2)
%HOMOGRAPHYTRANSFORM Calculate the homography transform that transforms points x1
%into points x2.
%   x1 and x2 must be 3xN matrices with homogeneous coordinates.

[m1,numPts] = size(x1);
[m2,    n2] = size(x2);

if (m1 ~= m2 || m1 ~= 3)
    error('both x1 and x2 must have three rows (homogeneous coordinates)');
else if (numPts ~= n2)
        error('x1 and x2 must have same number of elements');
    end
end

A = zeros(3*numPts,9);

O = [0 0 0];
for n = 1:numPts
    X = x1(:,n)';
    x = x2(1,n); y = x2(2,n); w = x2(3,n);
    A(3*n-2,:) = [  O  -w*X  y*X];
    A(3*n-1,:) = [ w*X   O  -x*X];
    A(3*n  ,:) = [-y*X  x*X   O ];
end

% apply SVD
[~,~,V] = svd(A,0);

% Extract homography
H = reshape(V(:,9),3,3)';

end


