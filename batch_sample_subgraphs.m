function subgraph_list_new = batch_sample_subgraphs(subgraph_list, batch_size, samples_per_frame, ...
    A_prev, A_curr, inview_ind, weights)
subgraph_size = size(subgraph_list, 2);

if ~isempty(weights)
    % if 'weights' is not empty compute the cumulative probability distribution for initially sampling a
    % specific cell
    weights_norm = cumsum(weights) / sum(weights);
end

% calculate the number of subgraph samples to create. If 'samples_per_frame' is less than
% 'samples_per_frame + batch_size' only sample the remainder number of subgraphs
remainder_subgraphs = samples_per_frame - size(subgraph_list, 1);
if remainder_subgraphs > batch_size
    remainder_subgraphs = batch_size;
end

% variables used during sampling
u_neighbors_init = sparse(false(size(A_prev, 1), 1));
group_logical_not_init = true(size(A_prev, 1), 1);

subgraph_list_new = zeros(0, subgraph_size);
while size(subgraph_list_new, 1) < remainder_subgraphs
    num_samples = size(subgraph_list_new, 1) + 1;
    subgraph_list_new = [subgraph_list_new; zeros(remainder_subgraphs - num_samples + 1, subgraph_size)];     

    sample = num_samples;
    while sample <= remainder_subgraphs
        if isempty(weights)
            % if weights is empty sample using a uniform distribution
            group0 = inview_ind(randi(length(inview_ind)));
        else
            % if weights is not empty sample the initial cell weighted by 'weights' (the standard error)
            group0 = inview_ind(find(rand <= weights_norm, 1));
        end

        % sample a subgraph
        group = sample_subgraph(group0, A_prev, A_curr, subgraph_size, u_neighbors_init, group_logical_not_init);

        % if the size of the subgraph matches 'subgraph_size' add it to the list of new subgraphs
        if size(group, 2) == subgraph_size
            subgraph_list_new(sample,:) = group;
            sample = sample + 1;
        end
    end

    % remove subgraphs which are non unique
    subgraph_list_new = unique(subgraph_list_new, 'rows');

    % exclude subgraphs which are already present in 'subgraph_list'
    exclude_ind = ismember(subgraph_list_new, subgraph_list, 'rows');
    subgraph_list_new(exclude_ind,:) = [];
end
end

function group = sample_subgraph(group0, A_prev, A_curr, subgraph_size, u_neighbors, group_logical_not)

group = group0;

group_logical_not(group) = false;

while size(group, 2) < subgraph_size
    % neighbors of the group in the previous frame
    g_neighbors_prev = u_neighbors | A_prev(:,group(end));

    % neighbors of the group in the current frame
    g_neighbors_curr = u_neighbors | A_curr(:,group(end));
    
    % neighbors of the group in both frames (excluding neighbors which are already in the
    % group)
    u_neighbors = g_neighbors_prev & g_neighbors_curr & group_logical_not;
    
    % if the number of neighbors is zero, then we exit (group cannot be used)
    if nnz(u_neighbors) == 0
        break;
    end

    % get indices of the group's neighbors
    temp  = find(u_neighbors);

    % add a random neighbor to the group
    group(end+1) = temp(randi(length(temp)));
    group_logical_not(group(end)) = false;
end

group = sort(group); % sort the group by the cell index
end