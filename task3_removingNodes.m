%task 3 - remove nodes.
clear all
clc

R = 0.6;                    % taken from task 2, approx 80% gets infected with small world network
N = 20;                     % number of nodes
probInfect = 1/(1 + 1/R);   % prob to remove edge
f = (0.05:0.05:0.2);        % fraction of removed nodes
nbrIterations = length(f);
maxSizes_smallWorld = zeros(2, nbrIterations);
maxSizes_Preferential = zeros(2, nbrIterations);

% Small world parameters
neighbours = 2;         % nbr of closest neighbours connected to
p_smallWorld = 0.4;     % prob of connecting shortcut.
removeProb = 0;
% Preferential network parameters
nInitial = 5;           % initial nbr of Nodes
nFinal = N;             % final nbr of nodes
m = 1;                  % number of added edges for each added node
p_preferential = 0.5;                % prob of edge between initial nodes

% create small world
% [x_smallWorld, y_smallWorld, A_smallWorld, XY_smallWorld] = CreateSmallWorldNetwork(N, p_smallWorld, neighbours, removeProb);
%create preferential
% [x_preferential, y_preferential, A_preferential, XY_preferential] = CreatePreferentialNetwork(nInitial, nFinal, p_preferential, m);
% figure(1)
% gplot(A_smallWorld, XY_smallWorld', '*-')
% figure(2)
% gplot(A_preferential, XY_preferential', '*-')


% for i_step = 1:nbrIterations
for i_step = 1:1
    nbrRemovedNodes = round(f(i_step)*N);
    
    %create networks
    % small world
    tic
    [x_sw1, y_sw1, A_sw1, XY_sw] = CreateSmallWorldNetwork(N, p_smallWorld, neighbours, removeProb);
    A_sw2 = A_sw1;      % A_sw1 for random removal, A_sw2 for degree removal
    toc
    tic
    % preferential - SOMETIMES IT GETS STUCK!!!
    [x_p1, y_p1, A_p1, XY_p] = CreatePreferentialNetwork(nInitial, nFinal, p_preferential, m);
    A_p2 = A_p1;
    toc
    
    %remove fraction of nodes
    %remove random Nodes
    removeNodes = randperm(N, nbrRemovedNodes);
    A_sw1(removeNodes,:) = 0;
    A_sw1(:, removeNodes) = 0;
    A_p1(removeNodes,:) = 0;
    A_p1(:, removeNodes) = 0;
        
    %remove nodes with highest degree
    degree_sw = full(sum(A_sw2, 2));
    degree_p = full(sum(A_p2, 2));
    removeNodes_sw = max(degree_sw, nbrRemovedNodes);   NOT WORKING!!! START HERE
    removeNodes_p = max(degree_p, nbrRemovedNodes);
    
    
    %remove edges that pass the recovery check
%     indexToRemove = rand(1,length(x)) < 1-probInfect;
%     xToRemove = [x(indexToRemove), y(indexToRemove)];
%     yToRemove = [y(indexToRemove), x(indexToRemove)];
%     ind = sub2ind([N,N], xToRemove, yToRemove);
%     B(ind) = 0;
%     
%     %Calculate biggest cluster
%     % maxCluster = GetBiggestCluster(A);
%     allPaths = GetAllPathLengths(B, N);
%     
%     GetGroups(N, allPaths);
    
end