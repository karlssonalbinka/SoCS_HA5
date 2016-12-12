%task 2 the SIR MODEL
clear all
clc

% Create small world network
N = 100;                 % number of nodes
neighbours = 2;         % nbr of closest neighbours connected to
p = 0.4;                % prob of connecting shortcut.
removeProb = 0;
[x, y, A, coord_XY] = CreateSmallWorldNetwork(N, p, neighbours, removeProb);
old_A = A;

R_all = 0.01:0.02:1.5;
averageIterations = 40;
nbrIterations = length(R_all);
maxSizes = zeros(1, nbrIterations);
%iterate over infection rate
for i = 1:nbrIterations
    i
    probInfect = 1/(1 + 1/R_all(i));
    maxAverage = zeros(1, averageIterations);
    for i_average = 1:averageIterations
        B = old_A;
        %remove edges that pass the recovery check
        for j = 1:N
            for k = 1:j
                removeEdge = rand(1) < 1-probInfect;
                if( removeEdge )
                    B(j, k) = 0;
                end
            end
        end
        B = B+B';
        B = B == 2;
        
        %Calculate biggest cluster
        % maxCluster = GetBiggestCluster(A);
        allPaths = GetAllPathLengths(B, N);
        
        maxAverage(i_average) = GetGroups(N, allPaths);
    end
    maxSizes(i) = sum(maxAverage)/averageIterations;
end
maxSizes = maxSizes/N;

plot(R_all, maxSizes)
title('largest group as function of R')
xlabel(['amll world network, N = ' num2str(N) ' p = ' num2str(p)]); 

%plot initial network and network with removed edges.
% hold off
% figure(1)
% gplot(old_A, coord_XY', 'r')
% hold on
% gplot(A, coord_XY', '*-b')
% title('small world with shortcuts')
% xlabel(['N = ' num2str(N) ', p = ' num2str(p)])

%% test get kargest cluster
clear all
clc
XY = [1 1 2 3 3; 0 2 1 0 2]; 
B = [0 1 0 0 0;...  %node 1 - lower left corner
    1 0 0 0 0; ...  %node 2 - upper left corner
    0 0 0 1 1;...   %node 3 - center node
    0 0 1 0 1;...   %node 4 - upper right corner
    0 0 1 1 0];     %node 5 - lower right corner
gplot(B, XY', '-*')
axis([0 4 -1 3])
N = length(B);

allPaths = GetAllPathLengths(B, N);

finalGroups = zeros(1,N);       % keep track of the groups we find
maxGroup = 0;                   % if a groupsize of N/2 is found, we can't find any bigger - stop while loop
halfNbrNodes = N/2;
while( sum(sum(finalGroups, 1)) < N && maxGroup <= halfNbrNodes )
    
    startOfGroup = find(sum(finalGroups,1) == 0, 1);            % get first node in new group
    group = allPaths(startOfGroup,:) == 1;                      % keep track of all nodes that are in this group
    nodesToCheck = allPaths(startOfGroup,:) == 1;               % what new nodes should we check where they leed to?
    while( sum(nodesToCheck) > 0 && sum(group(1,:)) < N)
        newNodesToCheck = 0;
        
        pathToNodes = sum(allPaths(nodesToCheck,:) == 1, 1);    % new nodes we can reach from the ones under investigation
        pathToNodes = pathToNodes > 0;
        newNodes = (pathToNodes ~= group).*(group~=1);          % have we found any new nodes that haven't already ben investigated
        
        if( sum(newNodes) > 0 )                                 % if new nodes - update nodes to check
            group = group + newNodes;
            newNodesToCheck = newNodes;
        end
        
        nodesToCheck = newNodesToCheck == 1;
    end
    
    if( sum(group) > maxGroup )                                 % update the max size of groups
        maxGroup = sum(group);
    end
    finalGroups = [finalGroups; group];                         % store the group configuration found
end