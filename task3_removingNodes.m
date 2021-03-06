%task 3 - remove nodes.
clear all
clc

R = 4;                    % taken from task 2, 100% gets infected with small world network
N = 100;                     % number of nodes
probInfect = 1/(1 + 1/R);   % prob to remove edge
f = (0.05:0.025:0.8);        % fraction of removed nodes
nbrIterations = length(f);
nbrAvg = 5;

%initiate variables
maxSizes_smallWorld     = zeros(2, nbrIterations);
maxSizes_Preferential   = zeros(2, nbrIterations);
maxGroup_sw1            = zeros(1,nbrIterations);
maxGroup_sw2            = zeros(1,nbrIterations);
maxGroup_p1             = zeros(1,nbrIterations);
maxGroup_p2             = zeros(1,nbrIterations);
sw1_avg                 = zeros(1, nbrAvg);
sw2_avg                 = zeros(1, nbrAvg);
p1_avg                  = zeros(1, nbrAvg);
p2_avg                  = zeros(1, nbrAvg);

% Random network parameter
p = 0.3;                % prob of edge between initial nodes

% Preferential network parameters
nInitial = 5;           % initial nbr of Nodes
nFinal = N;             % final nbr of nodes
m = 3;                  % number of added edges for each added node
p_preferential = 0.5;   % prob of edge between initial nodes

for i_step = 1:nbrIterations
    disp([num2str(i_step) ' out of ' num2str(nbrIterations)]);
  nbrRemovedNodes = round(f(i_step)*N);
  index_1 = zeros(1, nbrRemovedNodes);
  index_2 = zeros(1, nbrRemovedNodes);
    for avg_i = 1:nbrAvg
        
        %create networks
        % Random (obs, denoted sw1 and sw2)
        [x_sw1, y_sw1, A_sw1, XY_sw] = CreateRandomNetwork(N, p);
        A_sw2 = A_sw1;      % A_sw1 for random removal, A_sw2 for degree removal
        x_sw2 = x_sw1;
        y_sw2 = y_sw1;
        % preferential
        [x_p1, y_p1, A_p1, XY_p] = CreatePreferentialNetwork(nInitial, nFinal, p_preferential, m);
        A_p2 = A_p1;
        x_p2 = x_p1;
        y_p2 = y_p1;
        
        %remove fraction of nodes
        %remove random Nodes
        removeNodes = randperm(N, nbrRemovedNodes);   %Not working in matlab 2011...
%             removeNodes = ceil(rand(1, nbrRemovedNodes)*N);     %SWITCH TO THE OTHER ONE IN SCHOOL
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
            index_sw = find(degree_sw == sw_sorted(i));
            index_p  = find(degree_p == p_sorted(i));
            A_sw2(index_sw,:) = 0;
            A_sw2(:,index_sw) = 0;
            A_p2(index_p,:) = 0;
            A_p2(:,index_p) = 0;
        end
        
        %remove edges that pass the recovery check
        A_sw1 = edgesToRecover(A_sw1, probInfect);
        A_sw2 = edgesToRecover(A_sw2, probInfect);
        A_p1 = edgesToRecover(A_p1, probInfect);
        A_p2 = edgesToRecover(A_p2, probInfect);
        
        %     %Calculate biggest Group
        allPaths_sw1 = GetAllPathLengths(A_sw1, N);
        allPaths_sw2 = GetAllPathLengths(A_sw2, N);
        allPaths_p1  = GetAllPathLengths(A_p1, N);
        allPaths_p2  = GetAllPathLengths(A_p2, N);
        sw1_avg(avg_i) = GetGroups(N, allPaths_sw1);
        sw2_avg(avg_i) = GetGroups(N, allPaths_sw2);
        p1_avg(avg_i)  = GetGroups(N, allPaths_p1);
        p2_avg(avg_i)  = GetGroups(N, allPaths_p2);
    end
    maxGroup_sw1(i_step) = sum(sw1_avg)/(nbrAvg*N);
    maxGroup_sw2(i_step) = sum(sw2_avg)/(nbrAvg*N);
    maxGroup_p1(i_step) = sum(p1_avg)/(nbrAvg*N);
    maxGroup_p2(i_step) = sum(p2_avg)/(nbrAvg*N);
end

figure(1)
hold off
plot(f, maxGroup_sw1, 'b')
hold on
plot(f, maxGroup_sw2, 'r')
title('random network');
legend('random network - random', 'random network - highest deg');
xlabel(['fraction nodes removed - (N = ' num2str(N) ', R = ' num2str(R) ')'])

hold off
figure(2)
plot(f, maxGroup_p1, 'b')
hold on
plot(f, maxGroup_p2, 'r')
title('Preferential (randomly or highest node removed)');
legend('preferential - random','preferential - highest degree');
xlabel(['fraction nodes removed - (N = ' num2str(N) ', R = ' num2str(R) ')'])