function H = affineTransform(x1, x2)
%AFFINETRANSFORM Calculate the affine transform that transforms points x1
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

A = zeros(2 * numPts, 6);
b = zeros(2 * numPts, 1);

for n = 1:numPts
    A(2*n-1,:) = [ x1(:,n)' 0 0 0 ];
    A(2*n  ,:) = [ 0 0 0 x1(:,n)' ];
    
    b(2*n-1,1) = x2(1,n);
    b(2*n  ,1) = x2(2,n);
end

x = A\b;
H = [ reshape(x,3,2)'; 0 0 1 ]; 

end


