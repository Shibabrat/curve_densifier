function flag = hobsonVariant(traj,alpha,delta,dalpha) 

% Checks to insert points according to Hobson-Variant Criterion (Mancho et
% al., 2003) This combines the distance criteria using parameter tau and
% local radii of curvature using parameter alpha and delta

traj1 = traj(1,:);
traj2 = traj(2,:);
traj3 = traj(3,:);

x = [traj1;traj2;traj3];
 
% distance check
dj = norm(traj1-traj2,2);
if (dj > delta)
    flag = 1;
else
    flag = 0;  % points are within minimum distance and of the order of integration errror
    return
end

% angle check
xbar = x(2,:)  + ((x(2,:) - x(3,:))/(abs((x(2,:) - x(3,:))))).*(abs(x(2,:) - x(1,:)));

alphaj = 2*asin((xbar-x(1,:))/(2*(x(2,:)-x(1,:))));

productj = alphaj*dj;

if ((alphaj > alpha) || (productj > dalpha))
	flag = 1;
else
	flag = 0;
end

return
