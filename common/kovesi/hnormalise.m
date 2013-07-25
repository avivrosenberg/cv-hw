function nx = hnormalise(x)
%HNORMALISE % Normalize homogeneous coordinates so that w = 1

    [rows,~] = size(x);
    nx = x;

    % Find the indices of the points that are not at infinity
    finiteind = find(abs(x(rows,:)) > eps);

    for r = 1:rows-1
        nx(r,finiteind) = x(r,finiteind)./x(rows,finiteind);
    end
    nx(rows,finiteind) = 1;

end

