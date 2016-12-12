function [x, y, A, coordinates ] = CreateRandomNetwork( N, p )
% Create a network with random edges
%   N = nbr of nodes
%   p = prob of edges between nodes


k = 1:N;
randConnections = rand(N,N);
randConnections = triu(randConnections);
randConnections = randConnections + randConnections'; %make symmetric
[x,y] = find(randConnections < p);
z = ones(1,length(x));
A = sparse(x,y,z, N, N);



%Create circular points
angleStep = 2*pi/N;
coordinates = [cos(angleStep*k); sin(angleStep*k)];

end

