function V = vertices_to_array(V)
for ii = 1:length(V)
    temp = {V{ii}.coords}';
    temp(cellfun(@isempty, temp)) = {[nan, nan, nan]};
    temp = cell2mat(temp);
    V{ii} = temp;
end
end

