% RUN_CURVE_DENSIFIER Implement point insertion and interpolation for
% densifying curves. Point insertion can be done based on distance
% criterion or Hobson variant criterion (see Mancho, A. et al., 2003). The
% interpolation is done using cubic spline function in MATLAB.



clear all;close all;clc

% diary;

%parameters and global variable definition 
% global ALPHA BETA y0 EPS ODEFUNC epsilon
% ALPHA = 0.1;
% EPS = 0.5;
% BETA = 1/sqrt(1 + 4*ALPHA);
% y0 = (1+sqrt(1+4*ALPHA))/2;
% % ODEFUNC = 'kinematicModelHurricane';
% ODEFUNC = 'duffing'

% Parameters for insertion and interpolation of points
% tau = 0.001;		%distance tolerance for points on the manifolds
tau = 2e-2;
delta = 1e-6;		
alpha = 3e-1;       %local curvature tolerance for a segment 
dalpha = 1e-1;      %product of angle and distance tolerance

Delta = 4e-3;		%Initial separation of points of the linear manifold
lambdaMax = 1e-1;        %length of the eigen subspace at the starting time of the computation of the unstable manifold.

% t0 = 0; tf = 2*pi;
%Distinguished hyperbolic trajectory
% xDHT = -EPS*0.5*[0.0; 0.0]; %======> this looks suspicious

%Unstable and stable directions
%eUt = (1/sqrt(2))*[1; 1];
% eUt = [0.6796;    0.7336];
% eSt = (1/sqrt(2))*[-1; 1];

% ship roll model, Fig 5c, Case G
% xDHT = [-0.8575;   -0.0300];
% eSt = [-0.6003; 0.7998];
% eUt = [0.6037; 0.7972];

xDHT = [0.8765;  -2e-4];
eSt = [-0.6007; 0.7995];
eUt = [0.6037; 0.7972];
omegaN = 0.62;
omegaE = 0.527;
omegaBar = omegaE/omegaN;
TBar = (2*pi)/omegaBar;
t0 = 0; tf = 3*TBar;


% Test using Duffing equation's HT
% epsilon = 0.01;
% t0 = 0; tf = 6*pi;
% xDHT = -epsilon*0.5*[sin(t0); cos(t0)];
% eUt = (1/sqrt(2))*[1; 1];
% eSt = (1/sqrt(2))*[-1; 1];


flag = 1;		%1:forward time,-1:Backward time
maniFlag = 1;   %1: unstable 
et = eUt;		%eigendirection
dT = flag*((tf-t0)/3.0);

%Initialization
manifold = [];
trajectories = [];

% tf = flag*tf;
% dT = flag*dT;

% Setting the time slice vector
if maniFlag == 1,
    tVec	=	t0:dT:tf;
elseif maniFlag == -1,
    tVec	=	tf:dT:t0;
end

paramVec = [tf; tau; Delta; lambdaMax];
lambda = [-lambdaMax:Delta:lambdaMax];
% tVec	=	[t0:dT:tf];
manifold=cell(length(tVec),1);

%Generating initial manifold 
for i = 1:length(lambda)
	maniOld(i,:) = xDHT+lambda(i)*et;
	linMani(i,:) = xDHT+lambda(i)*et;
end

% maniOld(1:1000,:) = [-0.88*ones(1000,1), linspace(-0.88,0.88,1000)'];
% linMan(1:1000,:) = [-0.88*ones(1000,1), linspace(-0.88,0.88,1000)'];

k = 1;
manifold{k} = linMani;
tic;
%Obtaining the manifolds for all the time steps, this is also the simple
%method (Hobson, 1993)
% while (k < (length(tVec))),
% % 	for i = 1:size(maniOld,1)	
% % 		[tTemp, xTemp] = ode45(ODEFUNC, [tVec(k) tVec(k+1)],manifold{k}(i,:));
% %         trajectories(i,:) = xTemp(end,:);
% %     end
%     [xTemp, ~] = mex_integration([tVec(k) tVec(k+1)], ...
%         manifold{k}, maniFlag);
%     trajectories = xTemp;
% 
% 	manifold{k+1} = trajectories;
% 	maniOld = trajectories;
% 	trajectories =  [];
% 	k = k + 1;
% end
% manifold{k} = maniOld;

% for ii = 1:length(tVec)
%    plot(manifold{ii}(:,1),manifold{ii}(:,2),'r.');hold on
%    waitforbuttonpress
% end

% save('manifold.mat', 'manifold');

%%Insertion and interpolation of points
iterate = 1;
while (iterate < length(tVec)),

    % Just evolving the segment for the next time interval
    tic;
    [xTemp, ~] = mex_integration([tVec(iterate) tVec(iterate+1)], ...
            manifold{iterate},maniFlag);
    manifold{iterate+1} = xTemp;
    toc
    
	finished = false;
    oldTemp = [];
    newTemp = [];
	while ~finished,	% until all points pass the insertion criteria
		oldTemp = manifold{iterate};
		newTemp = manifold{iterate+1};
		flagVec = zeros(size(oldTemp,1),1);
		flagVec(1) = 0;
		oldInsert = [];
        
        tic;
        % Step through the old manifold points and check the distance
        for k=1:size(newTemp,1)-2	
			flagVec(k+1) = distCriterion(newTemp(k,:),newTemp(k+1,:),tau);
%             flagVec(k+1) = hobsonVariant(newTemp(k:k+2,:),alpha,delta,tau); 
			if (flagVec(k+1) == 1)	%Insert the point in the old manifold
				oldInsert = [oldInsert; insertPt(oldTemp(k,:),oldTemp(k+1,:))];
            end
        end
		toc
                
        
% 		for k =1:size(oldInsert,1)	%Iterate the inserted points and save in the newInsert 
% 			[tTemp, xTemp] = ode45(ODEFUNC, [tVec(iterate) tVec(iterate+1)], oldInsert(k,:));
% 			newInsert(k,:) = xTemp(end,:);
% 		end
 		
        tic;
        [xTemp, ~] = mex_integration([tVec(iterate) tVec(iterate+1)], ...
                oldInsert,maniFlag);
		newInsert = xTemp;
        toc
        
		%plot(newTemp(:,1),newTemp(:,2),'.g',newInsert(:,1),newInsert(:,2),'xg');hold on;
		
		%%Arranging points in ordered sequence along the new manifold
		newMani = [];
		oldMani = [];
		newMani = [newMani; newTemp(1,:)];
		oldMani = [oldMani; oldTemp(1,:)];
		j = 1;			%Counter for the inserted points
        
        tic;
		for k=1:(length(flagVec)-1)
			
			if (flagVec(k+1) == 0) 	%No point was inserted between k and k+1 point
				newMani = [newMani; newTemp(k+1,:)]; 
				oldMani = [oldMani; oldTemp(k+1,:)];
			else
				newMani = [newMani; newInsert(j,:); newTemp(k+1,:)];
				oldMani = [oldMani; oldInsert(j,:); oldTemp(k+1,:)];
				j = j + 1;
			end
				
		end
		toc
        
		if (length(find(flagVec ==1)) == 0)
			finished = true;
		else
			finished = false;	
		end
	
		manifold{iterate} = oldMani;
		manifold{iterate+1} = newMani;
		size(newMani,1)
		%plot(newMani(:,1),newMani(:,2),'or');
	
    end

    iterate = iterate + 1
%     trajectories = [];
    
    %Evolve the new manifold to the next step, this is the 
%     for i = 1:size(newMani,1)
%         [tTemp, xTemp] = ode45(ODEFUNC, [tVec(iterate) tVec(iterate+1)], newMani(i,:));
%         trajectories(i,:) = xTemp;
%     end

    
%     
%     manifold{iterate} = trajectories;
    
end

runtime = toc;

% save;

diary off

