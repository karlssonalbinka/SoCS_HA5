function [ x, y, A, XY ] = CreateSmallWorldNetwork( n, p, c, removeProb)
%Create a small world network (so far for task 2)
% N = nbr nodes
% p = prob of creating new shortcut
% c = nbr edges to neighbours to start with
% removeProb = probabillity to remove edges before adding shortcuts

nbrEdges = n*c;

% create circle coordinates
angleStep = 2*pi/n;
k = 1:n;
XY = [cos(angleStep*k); sin(angleStep*k)];

x = zeros(1, nbrEdges);
y = zeros(1, nbrEdges);
z = ones(1, nbrEdges);
for i = 1:n
    x(i*c-1:i*c) = 1+i:i+ c;
    y(i*c-1:i*c) = i;
end
%remove some random edges;
randNbr = rand(1,nbrEdges);
z(randNbr < removeProb) = 0;
x = mod(x-1,n) + 1;

%make shortcuts
%for each edge maybe place a new edge at random nodes
randNbr = rand(1,nbrEdges);
nbrNewEdges = sum(randNbr < p);
newEdgesXY = ceil(rand(2,nbrNewEdges)*n);
x = [x, newEdgesXY(1,:)];
y = [y, newEdgesXY(2,:)];
z = [z, ones(1,nbrNewEdges)];
A = sparse(x', y', z');
A(A+A'>0) = 1;


end

