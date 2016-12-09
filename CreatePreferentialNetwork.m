function [ x, y, A, XY ] = CreatePreferentialNetwork(nInitial, nFinal, p, m )
% Create a preferential network

% nInitial = initial nbr nodes
% nFinal = final number of nodes
% p = chance of edges between initial nodes
% m = number of new edges connected to each new node

nbrTimesteps = floor(nFinal-nInitial);

% create circle coordinates
angleStep = 2*pi/nFinal;
k = 1:nFinal;
XY = [cos(angleStep*k); sin(angleStep*k)];

% initialize configuration
A = sparse(nFinal,nFinal);
x = 1:nInitial;
y = 1:nInitial;
for i = 1:nInitial
    randCoord = (rand(1,i) < p);
    coordX(1:i) = x(i);
    coordY = y(randCoord);
    A(coordX', coordY') = 1;
end
A = A+A';   %make symmetric

for i_timestep = 1:nbrTimesteps
    i_timestep

    degree = full(sum(A,2));      %calculate degree
    totalDegree = sum(degree);
    degree = degree/totalDegree;

    %create added degrees for randoming out later on.
    cumulatedDegree = zeros(1,nInitial+(i_timestep-1)*m);
    for i = 1:nInitial+i_timestep;
        cumulatedDegree(i) = sum(degree(1:i));
    end

    %random what nodes to connect the new node to
    newNodes = zeros(1,m);
    newNode = 0;
    for i = 1:m
        while find(newNode == newNodes,1)
            newNode = find(cumulatedDegree > rand(1,1), 1);
        end
        newNodes(i) = newNode;
    end
    A(i_timestep + nInitial, newNodes) = 1;
    A(newNodes, i_timestep + nInitial) = 1;
end

end

