function xInsert = insertPt(x1,x2)

%   INSERTPT uses the MATLAB's built-in cubic spline interpolation for
%   insertion of a point between x1 and x2;

    temp = [x1;x2];
    pp = spline(linspace(0,1,2),temp');
    yy = ppval(pp,linspace(0,1,3));

    xInsert = yy(:,2)';

end