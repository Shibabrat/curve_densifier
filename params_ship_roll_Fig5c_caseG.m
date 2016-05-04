% ship roll model, Fig 5c, Case G
% xDHT = [-0.8575;   -0.0300];
% eSt = [-0.6003; 0.7998];
% eUt = [0.6037; 0.7972];

[xDHT,eUt, eSt] = get_hyp_traj_lin_subspace([-0.8561; 0]);
% xDHT = [0.8765;  -2e-4];
% eSt = [-0.6007; 0.7995];
% eUt = [0.6037; 0.7972];
% 
omegaN = 0.62;
omegaE = 0.527;
omegaBar = omegaE/omegaN;
TBar = (2*pi)/omegaBar;
t0 = 0; tf = 3*TBar;


flag = 1;          %1:forward time,-1:Backward time
maniFlag = 1;      %1: unstable 
branchMani = 1;     %1: in the same direction of eigenvector, -1: 180 rotation of eigenvector
et = eUt;           %eigenvector
dT = flag*((tf-t0)/3.0);


