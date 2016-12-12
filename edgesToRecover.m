function [ A ] = edgesToRecover( A, probInfect )
%See what edges should be removed thanks to nodes recovering
% A = adjacency matrix
% R = infection/recovery rate

N = length(A);
[x, y] = find(A == 1);

indexToRemove = rand(1,length(x)) < 1-probInfect;
xToRemove = [x(indexToRemove), y(indexToRemove)];
yToRemove = [y(indexToRemove), x(indexToRemove)];
ind = sub2ind([N,N], xToRemove, yToRemove);
A(ind) = 0;


end

