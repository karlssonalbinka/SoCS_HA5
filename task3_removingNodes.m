%task 3 - remove nodes.
clear all
clc

R = 0.6;                    % taken from task 2, approx 80% gets infected with small world network
N = 40;                     % number of nodes
probInfect = 1/(1 + 1/R);   % prob to remove edge
f = (0.05:0.05:0.5);        % fraction of removed nodes
nbrIterations = length(f);
nbrAvg = 10;

%initiate variables
maxSizes_smallWorld = zeros(2, nbrIterations);
maxSizes_Preferential = zeros(2, nbrIterations);
maxGroup_sw1 = zeros(1,nbrIterations);
maxGroup_sw2 = zeros(1,nbrIterations);
maxGroup_p1 = zeros(1,nbrIterations);
maxGroup_p2 = zeros(1,nbrIterations);
sw1_avg = zeros(1, nbrAvg);
sw2_avg = zeros(1, nbrAvg);
p1_avg = zeros(1, nbrAvg);
p2_avg = zeros(1, nbrAvg);

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


for i_step = 1:nbrIterations
    for avg_i = 1:nbrAvg
    nbrRemovedNodes = round(f(i_step)*N);
    index_1 = zeros(1, nbrRemovedNodes);
    index_2 = zeros(1, nbrRemovedNodes);
    
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
%     removeNodes = randperm(N, nbrRemovedNodes);   Not working in matlab
%     2011...
    removeNodes = ceil(rand(1, nbrRemovedNodes)*N);     %SWITCH TO THE OTHER ONE IN SCHOOL
    A_sw1(removeNodes,:) = 0;
    A_sw1(:, removeNodes) = 0;
    A_p1(removeNodes,:) = 0;
    A_p1(:, removeNodes) = 0;
        
    %remove nodes with highest degree
    degree_sw = full(sum(A_sw2, 2));
    degree_p = full(sum(A_p2, 2));
    sw_sorted = sort(degree_sw, 'descend');
    p_sorted = sort(degree_p, 'descend');
    for i = 1: nbrRemovedNodes
        index_1(i) = find(degree_sw == sw_sorted(i), 1);
        index_2(i) = find(degree_p == p_sorted(i), 1);
    end
    A_sw2(index_1, :) = 0;
    A_sw2(:, index_1) = 0;
    A_p2(index_2, :) = 0;
    A_p2(:, index_2) = 0;
    
    
    
    %remove edges that pass the recovery check
    A_sw1 = edgesToRecover(A_sw1, probInfect);
    A_sw2 = edgesToRecover(A_sw2, probInfect);
    A_p1 = edgesToRecover(A_p1, probInfect);
    A_p2 = edgesToRecover(A_p2, probInfect);
    
%     %Calculate biggest Group
    % maxCluster = GetBiggestCluster(A);
    allPaths_sw1 = GetAllPathLengths(A_sw1, N);
    allPaths_sw2 = GetAllPathLengths(A_sw2, N); 
    allPaths_p1 = GetAllPathLengths(A_p1, N); 
    allPaths_p2 = GetAllPathLengths(A_p2, N); 
    sw1_avg(i_avg) = GetGroups(N, allPaths_sw1);
    sw2_avg(i_avg) = GetGroups(N, allPaths_sw2);
    p1_avg(i_avg) = GetGroups(N, allPaths_p1);
    p2_avg(i_avg) = GetGroups(N, allPaths_p2);
    end
    maxGroup_sw1(i_step) = sum(sw1_avg)/nbrAvg;
    maxGroup_sw2(i_step) = sum(sw2_avg)/nbrAvg;
    maxGroup_p1(i_step) = sum(p1_avg)/nbrAvg;
    maxGroup_p2(i_step) = sum(p2_avg)/nbrAvg;
end

hold off
plot(f, maxGroup_sw1, 'b')
hold on
plot(f, maxGroup_p1, 'r')
plot(f, maxGroup_sw2, 'b--')
plot(f, maxGroup_p2, 'r--')
title('max group size as function over number of removed nodes (randomly or highest node removed)');
legend('small world - random','preferential - random', 'small world - highest deg','preferential - highest degree');
xlabel(['fraction nodes removed - (N = ' num2str(N) ', R = ' num2str(R) ')'])