function [ allPaths ] = GetAllPathLengths( A, maxIteration )
%OUT - Path: shortest paths for each node
%   Detailed explanation goes here
N = length(A);
accessUpperDiag = triu(true(N,N), 1);
allPaths = zeros(N,N);
i = 0;
C = 1;
while ( ~isempty(find(allPaths(accessUpperDiag) == 0, 1)) && i < maxIteration)
    i = i + 1;
    C = C*A;
    i_notTaken = allPaths == 0;
    i_existingPaths = C > 0;
    i_toUpdate = i_notTaken.*i_existingPaths == 1;
    allPaths(i_toUpdate) = i;
end

diag = 1:N;
diag = sub2ind([N,N], diag,diag);
allPaths(diag) = 0;

% nodeLongestPath = max(allPaths, [], 2);

return

end

