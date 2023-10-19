function remove_ind = check_conn_batch(A_curr, subgraph_list)

subgraph_size = size(subgraph_list, 2);

remove_ind = false(size(subgraph_list, 1), 1);

for ii = 1:size(subgraph_list, 1)
    
    % get subgraph adjacency matrix
    sg_A = A_curr(subgraph_list(ii,:), subgraph_list(ii,:));
 
    % use the Dulmage-Mendelsohn decomposition to find the number of connected components in
    % the subgraph
    [~, ~, r] = dmperm(sg_A | speye(subgraph_size));

    % if there exist more than one connected component (i.e. length(r) > 2) then remove that
    % subgraph
    if length(r) > 2
        remove_ind(ii) = true;
    end

end
end

