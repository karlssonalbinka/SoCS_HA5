% Task 4
clear all
clc

R = 4;                    % taken from task 2, 100% gets infected with small world network
N = 200;                     % number of nodes
probInfect = 1/(1 + 1/R);   % prob to remove edge
f = (0.05:0.025:0.8);        % fraction of removed nodes
nbrIterations = length(f);
nbrAvg = 1;

% Random network parameter
p = 0.5;                % prob of edge between initial nodes

% Preferential network parameters
nInitial = 5;           % initial nbr of Nodes
nFinal = N;             % final nbr of nodes
m = 3;                  % number of added edges for each added node
p_preferential = 0.5;   % prob of edge between initial nodes

%initiate variables
maxGroup_sw1            = zeros(1,nbrIterations);
maxGroup_sw2            = zeros(1,nbrIterations);
maxGroup_sw3            = zeros(1,nbrIterations);
maxGroup_p1             = zeros(1,nbrIterations);
maxGroup_p2             = zeros(1,nbrIterations);
maxGroup_p3             = zeros(1,nbrIterations);
sw1_avg                 = zeros(1, nbrAvg);
sw2_avg                 = zeros(1, nbrAvg);
sw3_avg                 = zeros(1, nbrAvg);
p1_avg                  = zeros(1, nbrAvg);
p2_avg                  = zeros(1, nbrAvg);
p3_avg                  = zeros(1, nbrAvg);


for i_step = 1:nbrIterations
    disp([num2str(i_step) ' out of ' num2str(nbrIterations)]);
  nbrRemovedNodes = round(f(i_step)*N);
  index_1 = zeros(1, nbrRemovedNodes);
  index_2 = zeros(1, nbrRemovedNodes);
  nbrNodes = N - nbrRemovedNodes;
    for avg_i = 1:nbrAvg
        
        %create networks
        % Random (obs, denoted sw1 and sw2)
        [x_sw1, y_sw1, A_sw1, XY_sw] = CreateRandomNetwork(N, p);
        A_sw2 = A_sw1;      % A_sw1 for random removal, A_sw2 for degree removal
        A_sw3 = A_sw1;      % for vaccinating through neighbour method

        % preferential
        [x_p1, y_p1, A_p1, XY_p] = CreatePreferentialNetwork(nInitial, nFinal, p_preferential, m);
        A_p2 = A_p1;        % A_p2 for degree removal
        A_p3 = A_p1;        % for vaccinating through neighbour method
        
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
        
        %remove nodes that are neighbours to a randomed node
        gotRemoved_sw3 = 0;
        gotRemoved_p3 = 0;
        for i = 1:nbrRemovedNodes
            %get neighbours
            randomNode = ceil(rand(1)*(N-i));
            neighbours = find(A_sw3(randomNode,:) == 1);
            test = 0;
            while( isempty(neighbours)  && test < 5);
                test = test + 1;
                randomNode = ceil(rand(1)*(N-i));
                neighbours = find(A_sw3(randomNode,:) == 1);
            end
            if( ~isempty(neighbours) )
                gotRemoved_sw3 = gotRemoved_sw3 +1;
                removeNode = neighbours( ceil(rand(1)*length(neighbours)) );
                A_sw3(removeNode,:) = [];
                A_sw3(:,removeNode) = [];
            end
            %get neighbours
            neighbours = find(A_p3(randomNode,:) == 1);
            test = 0;
            while( isempty(neighbours)  && test < 5);
                test = test + 1;
                randomNode = ceil(rand(1)*(N-i));
                neighbours = find(A_p3(randomNode,:) == 1);
            end
            if( ~isempty(neighbours) )
                gotRemoved_p3 = gotRemoved_p3 +1;
                removeNode = neighbours( ceil(rand(1)*length(neighbours)) );
                A_p3(removeNode,:) = [];
                A_p3(:,removeNode) = [];
            end

        end
        
        %remove edges that pass the recovery check
        A_sw1 = edgesToRecover(A_sw1, probInfect);
        A_sw2 = edgesToRecover(A_sw2, probInfect);
        A_sw3 = edgesToRecover(A_sw3, probInfect);
        A_p1 = edgesToRecover(A_p1, probInfect);
        A_p2 = edgesToRecover(A_p2, probInfect);
        A_p3 = edgesToRecover(A_p3, probInfect);
        
        %     %Calculate biggest Group
        allPaths_sw1 = GetAllPathLengths(A_sw1, N);
        allPaths_sw2 = GetAllPathLengths(A_sw2, N);
        allPaths_sw3 = GetAllPathLengths(A_sw3, N-gotRemoved_sw3);
        allPaths_p1  = GetAllPathLengths(A_p1, N);
        allPaths_p2  = GetAllPathLengths(A_p2, N);
        allPaths_p3  = GetAllPathLengths(A_p3, N-gotRemoved_p3);
        sw1_avg(avg_i) = GetGroups(N, allPaths_sw1);
        sw2_avg(avg_i) = GetGroups(N, allPaths_sw2);
        sw3_avg(avg_i) = GetGroups( N-gotRemoved_sw3, allPaths_sw3);
        p1_avg(avg_i)  = GetGroups(N, allPaths_p1);
        p2_avg(avg_i)  = GetGroups(N, allPaths_p2);
        p3_avg(avg_i)  = GetGroups(N-gotRemoved_p3 , allPaths_p3);
    end
    maxGroup_sw1(i_step) = sum(sw1_avg)/(nbrAvg*N);
    maxGroup_sw2(i_step) = sum(sw2_avg)/(nbrAvg*N);
    maxGroup_sw3(i_step) = sum(sw3_avg)/(nbrAvg*N);
    maxGroup_p1(i_step) = sum(p1_avg)/(nbrAvg*N);
    maxGroup_p2(i_step) = sum(p2_avg)/(nbrAvg*N);
    maxGroup_p3(i_step) = sum(p3_avg)/(nbrAvg*N);
end

figure(1)
hold off
plot(f, maxGroup_sw1, 'b')
hold on
plot(f, maxGroup_sw2, 'r')
plot(f, maxGroup_sw3, 'k')
title('random network');
legend('vaccinate random', 'vaccinate highest deg', 'vaccinate neighbour');
xlabel(['fraction nodes removed - (N = ' num2str(N) ', R = ' num2str(R) ')'])

figure(2)
hold off
plot(f, maxGroup_p1, 'b')
hold on
plot(f, maxGroup_p2, 'r')
plot(f, maxGroup_p3, 'k')
title('Preferential');
legend('vaccinate random','vaccinate highest degree', 'vaccinate neighbour');
xlabel(['fraction nodes removed - (N = ' num2str(N) ', R = ' num2str(R) ')'])
