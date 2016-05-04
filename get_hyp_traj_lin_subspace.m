function [xDHT,eUt, eSt] = get_hyp_traj_lin_subspace(xDHTGuess)

% GET_HYP_TRAJ_LIN_SUBSPACE Evolves a line segment near a guess saddle
% point to compute the equilibrium point in time dependent flow.

format long e

% xDHT = [0; 0];
% xDHTGuess = [-0.8561; 0];
% initialize the vector in random direction
et = rand(2,1);
et = et/norm(et);

segmentFor = cell(2,1);
segmentBack = cell(2,1);

H = 9.8362;
HBar = (0.73*pi)*(H/221.94);
epsilon = HBar;
omegaN = 0.62;
omegaE = 0.527;
omegaBar = omegaE/omegaN;
TBar = (2*pi)/omegaBar;
tI = 0; tF = 0.05*TBar;

% epsilon = 0.01;
% tI = 0; tF = 0.1*2*pi;

numPtsInitSeg = 5e3;
delta = 1e-3;
iterMax = 20;

lambda = linspace(-epsilon,epsilon,numPtsInitSeg+1);

iterate = 1;
eUt = et;
eSt = et;
xDHT = xDHTGuess;
while iterate < iterMax
    % generate initial segement 
    for i = 1:length(lambda)
        segmentFor{1}(i,:) = xDHT+lambda(i)*eUt;
        segmentBack{1}(i,:) = xDHT+lambda(i)*eSt;
    end

    % evolve the initial segment forward for unstable subspace
    tic;
    [xTemp, ~] = mex_integration([tI tF], ...
            segmentFor{1},1);
    segmentFor{2} = xTemp;
    toc;
    
    % evolve the initial segment backward for stable subspace
    tic;
    [xTemp, ~] = mex_integration([tF tI], ...
            segmentBack{1},-1);
    segmentBack{2} = xTemp;
    toc;
    
    % compute intersection point of the forward and backward segments 
    % Assumption is flow time in small enough to keep the segment linear
    numPtsSegment = size(xTemp,1);
    
    out = lineSegmentIntersect([segmentBack{2}(1,:),segmentBack{2}(end,:)], ...
        [segmentFor{2}(1,:),segmentFor{2}(end,:)]);

    % Initialize the next segment
    % update xDHT and the eigendirection
    xDHT = [out.intMatrixX; out.intMatrixY]

    if (out.intAdjacencyMatrix == 0), % evolved lines didn't intersect
        xDHT = xDHTGuess;
    end
    
    eUt = [segmentFor{2}(1,:)' - segmentFor{2}(end,:)'];
    eUt = eUt/norm(eUt)
    eSt = [segmentBack{2}(1,:)' - segmentBack{2}(end,:)'];
    eSt = eSt/norm(eSt)

    
    % Visual test for the algorithm
%     plot(segmentFor{2}(:,1),segmentFor{2}(:,2),'.r')
%     hold on
%     plot(segmentBack{2}(:,1),segmentBack{2}(:,2),'.b');
%     waitforbuttonpress
%     hold off
    
    iterate = iterate + 1
end


end


