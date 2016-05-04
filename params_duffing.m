% Test using Duffing equation's HT
epsilon = 0.01;
t0 = 0; tf = 6*pi;
xDHT = -epsilon*0.5*[sin(t0); cos(t0)];
eUt = (1/sqrt(2))*[1; 1];
eSt = (1/sqrt(2))*[-1; 1];


flag = 1;           %1:forward time,-1:Backward time
maniFlag = 1;       %1: unstable 
branchMani = -1;    %1: in the same direction of eigenvector, -1: 180 rotation of eigenvector
et = eUt;           %eigenvector
dT = flag*((tf-t0)/3.0);



