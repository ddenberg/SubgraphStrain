function [incF_list, exclude_group] = batch_inc_deformation_gradient(A_prev, cell_lineage_centroids, subgraph_list, t)

incF_list = zeros(3, 3, size(subgraph_list, 1));
exclude_group = false(size(subgraph_list, 1), 1);
% cell_neighbors = zeros(size(subgraph_list));
for ii = 1:size(subgraph_list, 1)
    
    % get subgraph adjacency matrix
    sg_A = A_prev(subgraph_list(ii,:), subgraph_list(ii,:));

    % get edges of subgraph
    [I,J] = find(triu(sg_A));
    edges = [I,J];
    
    % get the centroids of cells in the subgraph at the previous and current time point
    x_prev = cell_lineage_centroids(:,subgraph_list(ii,:),t-1).';
    x_curr = cell_lineage_centroids(:,subgraph_list(ii,:),t).';
    
    % create the lambda and lambda0 arrays
    lambda0 = x_prev(edges(:,1),:) - x_prev(edges(:,2),:);
    lambda = x_curr(edges(:,1),:) - x_curr(edges(:,2),:);

    % compute the incremental deformation gradient
    incF_temp = (lambda.' * lambda0) / (lambda0.' * lambda0);
    
    if any(isnan(incF_temp), 'all')
        % exclude subgraph if system is ill-conditioned (likely a co-planar set of cells)
        exclude_group(ii) = true;
    end

    incF_list(:,:,ii) = incF_temp;
end
end

