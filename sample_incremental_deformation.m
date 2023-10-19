function [incF_cells, cell_frequency, max_stde] = sample_incremental_deformation(G, cell_lineage_centroids, inview_cells, lineage, ...
                                                     final_frame, subgraph_size, samples_per_frame, batch_size, dropout)
% num_frames = length(G);
num_cells = size(inview_cells, 1);

% create output array for incremental deformation gradient
incF_cells = zeros(3, 3, num_cells, final_frame);
% Initialize to identity matrix
incF_cells(1,1,:,:) = 1;
incF_cells(2,2,:,:) = 1;
incF_cells(3,3,:,:) = 1;

% output array for cell frequency (# subgraphs include a cell per frame)
cell_frequency = zeros(num_cells, final_frame);
inview_ind = find(inview_cells(:,1));
for ii = 1:length(inview_ind)
    cell_frequency(inview_ind(ii),1) = 1; % initialize first frame index to 1
end

% output array for the maximum standard error per cell per frame
max_stde = nan(num_cells, final_frame);
max_stde(:,1) = 0;

% loop over all frame pairs (each frame pair is [t-1, t]) Frame pair [0, 1] is skipped because
% frame 0 does not exist
for t = 2:final_frame
    count_str = '';
    fprintf('Sampling Frame %d/%d... ', t, final_frame);

    % get indices of cells which are inview for frames [t-1, t]
    inview_ind = find(all(inview_cells(:,t-1:t), 2));

    % get indices of divisions which occur on frame t
    div_idx = find(lineage.div_time == t);

    % get graphs from frame t and t-1
    G_curr = G{t};
    G_prev = G{t-1}; 

    % if there are divisions in this frame, copy the history of the incremental deformations of the
    % mother cell to the daughter cells
    if ~isempty(div_idx)
        for ii = 1:length(div_idx)
            for jj = 1:t-1
                intF_temp = incF_cells(:,:,lineage.mother(div_idx(ii)),jj);
                incF_cells(:,:,lineage.daughters(div_idx(ii),1),jj) = intF_temp;
                incF_cells(:,:,lineage.daughters(div_idx(ii),2),jj) = intF_temp;
            end
            
            % append indices of daughter cells to the list of inview cells
            inview_ind = union(inview_ind, lineage.daughters(div_idx(ii),:));
            
            % add duaghter cells to the previous frame graph such that they are connected to
            % the same neighbors as the mother. These "ficticious" cells are added so that
            % subgraphs can be sampled which include the division event
            mother_edges = neighbors(G_prev, lineage.mother(div_idx(ii)));
            new_edges = [[lineage.daughters(div_idx(ii),1) * ones(length(mother_edges), 1); ...
                            lineage.daughters(div_idx(ii),2) * ones(length(mother_edges), 1)], [mother_edges; mother_edges]];
            G_prev = addedge(G_prev, new_edges(:,1), new_edges(:,2)); 
        end
    end
    % extract adjacency matrices from the graph objects and convert them to sparse logical arrays
    A_curr = logical(adjacency(G_curr));
    A_prev = logical(adjacency(G_prev));
    
    % create mean (incF_nodes_M) and standard deviation (incF_nodes_S) arrays to store the mean
    % and standard deviation of the cells' incremental deformation gradient for this frame pair
    incF_nodes_M = zeros(3, 3, num_cells);
    incF_nodes_S = zeros(3, 3, num_cells);

    % 'cell_frequency_temp' is used to hold the intermediate sample numbers for each cell until
    % sampling is completed.
    cell_frequency_temp = zeros(num_cells, 1);
    weights = [];
    if t == 2
        % if we are at the beginning of the movie there are no prior subgraph samples so the
        % subgraph list is initialized to an empty array
        subgraph_list = zeros(0, subgraph_size);
    else
        % if we are after the beginning of the movie there may be prior subgraphs which are
        % still connected and do not need to be re-sampled. This can drastically speed up the process 

        % remove a random fraction of subgraphs (with probability 'dropout') 
        % also remove subgraphs which contain cells which are no longer in view
        remove_ind = ~all(ismember(subgraph_list, inview_ind), 2) | rand(size(subgraph_list, 1), 1) < dropout;
        subgraph_list(remove_ind,:) = [];
        
        % check if each subgraph is connected and remove the ones which are not
        remove_ind = check_conn_batch(A_curr, subgraph_list);
        subgraph_list(remove_ind,:) = [];
        
        % compute the incremental deformation gradient for each subgraph. 'exclude_group'
        % indicates which subgraphs are ill-conditioned and are likely co-planar. These are
        % removed
        [incF_list, exclude_group] = batch_inc_deformation_gradient(A_prev, cell_lineage_centroids, subgraph_list, t);
        incF_list = incF_list(:,:,~exclude_group);
        subgraph_list = subgraph_list(~exclude_group,:);

        % update the cells' incremental deformation means, standard deviations and cell
        % frequencies for this frame pair
        [incF_nodes_M, incF_nodes_S, cell_frequency_temp] = accumulate_batch(incF_nodes_M, incF_nodes_S, incF_list, ...
            cell_frequency_temp, subgraph_list);
        
        % if there exists at least 2 subgraph samples for every cell, compute the standard error
        if min(cell_frequency_temp(inview_ind)) >= 2
            weights = squeeze(max(incF_nodes_S(:,:,inview_ind), [], [1,2])) ./ (cell_frequency_temp(inview_ind) - 1);
            weights = sqrt(weights ./ cell_frequency_temp(inview_ind));
        else
            weights = [];
        end
    end
    
    % while the size of the current list of subgraphs is less than the number of subgraphs per
    % frame we desire, add subgraphs in batches to the full list of subgraph
    while size(subgraph_list, 1) < samples_per_frame
        % sample a batch of subgraphs which are distinct from the previously sampled set of
        % subgraphs
        subgraph_list_new = batch_sample_subgraphs(subgraph_list, batch_size, samples_per_frame, ...
            A_prev, A_curr, inview_ind, weights);
        
        % compute the incremental deformation gradients for each subgraph
        [incF_list, exclude_group] = batch_inc_deformation_gradient(A_prev, cell_lineage_centroids, subgraph_list_new, t);
        incF_list = incF_list(:,:,~exclude_group);
        subgraph_list_new = subgraph_list_new(~exclude_group,:);
        
        % update the cells' incremental deformation means, standard deviations and cell
        % frequencies for this frame pair
        [incF_nodes_M, incF_nodes_S, cell_frequency_temp] = accumulate_batch(incF_nodes_M, incF_nodes_S, incF_list, ...
            cell_frequency_temp, subgraph_list_new);
        subgraph_list = cat(1, subgraph_list, subgraph_list_new);
        
        % if there exists at least 2 subgraph samples for every cell, compute the standard error
        if min(cell_frequency_temp(inview_ind)) >= 2
            weights = squeeze(max(incF_nodes_S(:,:,inview_ind), [], [1,2])) ./ (cell_frequency_temp(inview_ind) - 1);
            weights = sqrt(weights ./ cell_frequency_temp(inview_ind));
        end
        
        % update progress bar
        fprintf(repmat('\b', 1, length(count_str)-1))
        count_str = [num2str(100 * size(subgraph_list, 1) / samples_per_frame, '%.0f'), '%%'];
        fprintf(count_str);
    end
    
    % append intermediate quantities to the output variables
    if ~isempty(weights)
        max_stde(inview_ind,t) = weights;
    else
        max_stde(inview_ind,t) = nan(length(inview_ind), 1);
    end
    cell_frequency(:,t) = cell_frequency_temp;
    incF_cells(:,:,inview_ind,t) = incF_nodes_M(:,:,inview_ind);
    
    % update progress bar
    fprintf(repmat('\b', 1, length(count_str)-1));
    fprintf('Done.\n');
end

end

