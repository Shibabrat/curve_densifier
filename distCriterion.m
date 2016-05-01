function flag = distCriterion(traj1,traj2,tau) 

%Implements the distance criterion for adjacent trajectories 
%x_j and x_{j+1}- flag = 1:insert point

dj = norm(traj1-traj2,2);
if (dj >= tau)
    flag = 1;
else
    flag = 0;
end

return