function [incF_nodes_M, incF_nodes_V, node_frequency] = accumulate_batch(incF_nodes_M, incF_nodes_V, incF_list, ...
    node_frequency, subgraph_list)

for ii = 1:size(subgraph_list, 1)
    node_frequency(subgraph_list(ii,:)) = node_frequency(subgraph_list(ii,:)) + 1;
    for jj = 1:size(subgraph_list, 2)
        if node_frequency(subgraph_list(ii,jj)) == 1
            incF_nodes_M(:,:,subgraph_list(ii,jj)) = incF_list(:,:,ii);
        else
            temp = incF_nodes_M(:,:,subgraph_list(ii,jj));
            incF_nodes_M(:,:,subgraph_list(ii,jj)) = temp + (incF_list(:,:,ii) - temp) / node_frequency(subgraph_list(ii,jj));
            incF_nodes_V(:,:,subgraph_list(ii,jj)) = incF_nodes_V(:,:,subgraph_list(ii,jj)) + (incF_list(:,:,ii) - temp) .* ...
                (incF_list(:,:,ii) - incF_nodes_M(:,:,subgraph_list(ii,jj)));
        end
    end
end

end