% HW 5 - task 1

clc
clear all
A = load('communityExample.mat');
A = A.A;

k = sum(A,2);                   %nbr edges of each node

%########################## CALC B and Q once ############3
N = length(A);
allNodes = 1:N;
m = sum(sum(A))/2;                %number of edges

%creade modular matrix B
k_ij = zeros(N,N);
for i = 1:N
    k_ij(i,:) = k(i)*k;
end
B = A - k_ij/(2*m);

[eigVec, eigVal] = eig(B); 
[max_i, max_j] = find(eigVal == max(max(eigVal)));       %is it always the last element that is the largest eigenvalue?
maxEigVec = eigVec(:,max_i);

%assign group
s = ones(N,1);
s(maxEigVec < 0) = -1;

%Calculate Modularity
Q = s'*B*s/(4*m)
%#####################################################
%#####################################################
%####################WHILE LOOP TO GO THROUGH ALL SUB GROUPS###############
nbrGroups = 2;
groups = [allNodes(s>0); allNodes(s<0)];             %initializing full group
test = 0;
B_tmp = zeros(N,N);

while (nbrGroups > 0 && test < 4)
    test = test + 1;
    test
    
    newNbrGroups = 0;
    offset = 1;
    [edge, ~] = size(groups);
    for i_groups = 1:nbrGroups
        currentGroup = groups( edge-i_groups+1, groups(edge-i_groups+1,:)>0);
        
        N_tmp = length(currentGroup);
%         m = sum(sum( A(currentGroup, :) ))/2;                %number of edges on the sub group
        
        %create modular matrix B for this group
        B_tmp = zeros(N_tmp, N_tmp);
        for i = 1:N_tmp
            for j = 1:N_tmp
                B_tmp(i,j) = B(currentGroup(i),currentGroup(j)) - (currentGroup(i)==currentGroup(j))*sum(B(currentGroup(i),currentGroup));
            end
        end
        
        [eigVec, eigVal] = eig(B_tmp);
%         [max_i, max_j] = size(eigVal);       %is it always the last element that is the largest eigenvalue?
        maxEigVec = eigVec(:,end);
        
        %assign group
        s = ones(N_tmp,1);
        s(maxEigVec < 0) = -1;
        
        %Calculate Modularity
        Q_tmp = s'*B_tmp*s/(4*m)
        
        if(Q_tmp > exp(-7))
            newNbrGroups = newNbrGroups + 2;
            Q = Q + Q_tmp
%             g1 = zeros(1, length(groups));
%             tmp = false(1,N);
%             tmp(currentGroup) = true;
%             indexOfGroup = tmp'.*(s>0) == 1;
%             g1(1:sum(indexOfGroup)) = allNodes(indexOfGroup);
%             
%             g2 = zeros(1, length(groups));
%             tmp = false(1,N);
%             tmp(currentGroup) = true;
%             indexOfGroup = tmp'.*(s>0) == 1;
%             g2(1:sum(indexOfGroup)) = allNodes(indexOfGroup);
            g1 = zeros(1, length(groups));
            g2 = g1;
            g1(1:sum(s>0)) = currentGroup(s>0);
            g2(1:sum(s<0)) = currentGroup(s<0);
            groups = [groups; g1; g2];
        end
    end
    nbrGroups = newNbrGroups;

end
%%
%plot the mess
%create positions depending on group a node is assigned
interval = 5;
sepperation = 2*interval;
nbrGroup1 = sum(groups(end,:) > 0);
nbrGroup2 = sum(groups(end-1,:) > 0);
nbrGroup3 = N - nbrGroup1 - nbrGroup2;
x = zeros(1,N);
y = zeros(1,N);
x(groups(end,1:nbrGroup1)) = rand(1,nbrGroup1)*interval + sepperation/2;
x(groups(end-1,1:nbrGroup2)) = rand(1,nbrGroup2)*interval - sepperation/2;
x(groups(end-2,1:nbrGroup3)) = rand(1,nbrGroup3)*interval;
y(groups(end,1:nbrGroup1)) = rand(1,nbrGroup1)*interval;
y(groups(end-1,1:nbrGroup2)) = rand(1,nbrGroup2)*interval;
y(groups(end-2,1:nbrGroup3)) = rand(1,nbrGroup3)*interval + sepperation/2;

hold off
gplot(A, [x;y]')
hold on
plot(x(groups(end,1:nbrGroup1)), y(groups(end,1:nbrGroup1)), 'or', 'MarkerFaceColor', 'r')
plot(x(groups(end-1,1:nbrGroup2)), y(groups(end-1,1:nbrGroup2)), 'og', 'MarkerFaceColor', 'g')
plot(x(groups(end-2,1:nbrGroup3)), y(groups(end-2,1:nbrGroup3)), 'og', 'MarkerFaceColor', 'b')
title('task 1, grouping. Q=0.549677678601963')