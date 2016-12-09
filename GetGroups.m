function [ maxSize, finalGroups ] = GetGroups(N, allPaths  )
%Get the sepperation of groups.
% - Also calculate the largest group found

% N = nbr of nodes
% allPaths = all paths that nodes can take

finalGroups = zeros(1,N);       % keep track of the groups we find
maxGroup = 0;                   % if a groupsize of N/2 is found, we can't find any bigger - stop while loop
halfNbrNodes = N/2;
% outerWhile = 0;
while( sum(sum(finalGroups, 1)) < N && maxGroup <= halfNbrNodes )
%     outerWhile = outerWhile + 1
    
    startOfGroup = find(sum(finalGroups,1) == 0, 1);            % get first node in new group
    group = allPaths(startOfGroup,:) == 1;                      % keep track of all nodes that are in this group
    group(startOfGroup) = true;
    nodesToCheck = allPaths(startOfGroup,:) == 1;               % what new nodes should we check where they leed to?
%     innerWhile = 0;
    while( sum(nodesToCheck) > 0 && sum(group(1,:)) < N)
%         innerWhile = innerWhile + 1
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

maxSize = max(sum(finalGroups, 2));

end

