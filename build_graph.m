function [G, inview_cells] = build_graph(C, E)
G = cell(length(C), 1);

num_cells = max(cellfun(@length, C));

inview_cells = false(max(cellfun(@length, C)), length(C));
for ii = 1:length(C)
    edges = {E{ii}.cells}.';
    edges = cell2mat(edges);
    
    G{ii} = graph(edges(:,1), edges(:,2), [], num_cells);

    %Check to make sure all nodes have neighbors
    is_neighbored_cell = degree(G{ii}) > 2;
    inview_cells(is_neighbored_cell,ii) = true;  
end

end

