function export_vtk(writeFolder, fileFormat, C, V, inview_cells, cell_frequency, F_cells, Strain_cells, Stretch_cells, ...
    Eigenvector_cells, max_stde, cell_ids, frame_list)

% Write File for cells

if ~strcmp(fileFormat, 'ASCII') && ~strcmp(fileFormat, 'BINARY')
    error('Please enter a valid file format');
end

for frame = frame_list
    fprintf('Writing vtk %d/%d... ', frame, length(frame_list));
    
    filename = sprintf('%s/cells%04d.vtk', writeFolder, frame);

    fID = fopen(filename, 'w');

    fprintf(fID, '# vtk DataFile Version 3.0\n');
    fprintf(fID, '%s\n', 'Cell Mesh');
%     fprintf(fID, 'ASCII\n\n');
    fprintf(fID, [fileFormat, '\n\n']);
    fprintf(fID, 'DATASET POLYDATA\n' );

    fprintf(fID, 'POINTS %d double\n', length(V{frame}));

%     vertices = {V{frame}.coords}';
%     is_used_idx = ~cellfun(@isempty, vertices);
%     vertices(~is_used_idx) = {[nan nan nan]};
%     vertices = cell2mat(vertices);
    vertices = V{frame};
    if strcmp(fileFormat, 'ASCII')
        fprintf(fID, '%.12f\t%.12f\t%.12f\n', vertices'); %Write coordinates of nodes
    else
        fwrite(fID, vertices.', 'double' , 'b');
    end
    

    cell_verts = {C{frame}.vertices}.';
    cell_size = cellfun(@length, cell_verts);
    
%     output_cells = cell_size > 0;
    output_cells = inview_cells(:,frame);

    % Only output cells which have cell_size > 0
    fprintf(fID, '\nPOLYGONS\t%d\t%d\n', sum(output_cells), sum(cell_size(output_cells)) + sum(output_cells));

    for m = 1:length(C{frame})
        if output_cells(m)
            if strcmp(fileFormat, 'ASCII')
                format_spec = '%d';
                for ii = 1:cell_size(m)
                    format_spec = [format_spec, ' %d'];
                end
                format_spec = [format_spec, '\n'];

                fprintf(fID, format_spec, [cell_size(m), cell_verts{m} - 1]); %Write out element connectivity data
            else
                fwrite(fID, [cell_size(m), cell_verts{m} - 1], 'int', 'b');
            end
        end
    end
    %Note VTK uses a 0 based start point system vs MATLAB's 1 starting point

    fprintf(fID, '\nCELL_DATA\t%d\n', sum(output_cells));
    fprintf(fID, 'FIELD Label 16\n');
    
    %% Writing Number of samples per cell

    fprintf(fID, 'Samples_per_cell 1\t%d\tint\n', sum(output_cells));
    if strcmp(fileFormat, 'ASCII')
        fprintf(fID, '%d\n', cell_frequency(output_cells,frame).');
    else
        fwrite(fID, cell_frequency(output_cells,frame).', 'int', 'b');
    end
    
    %% Writing Cell IDs
    
    fprintf(fID, 'Cell_IDs 1\t%d\tint\n', sum(output_cells));
    if strcmp(fileFormat, 'ASCII')
        fprintf(fID, '%d\n', cell_ids(output_cells).');
    else
        fwrite(fID, cell_ids(output_cells).', 'int', 'b');
    end
    
    %% Writing Max_Stde
    
    fprintf(fID, 'Max_STDe 1\t%d\tdouble\n', sum(output_cells));
    temp = max_stde(output_cells,frame);
    if strcmp(fileFormat, 'ASCII')
        fprintf(fID, '%f\n', temp.');
    else
        fwrite(fID, temp.', 'double', 'b');
    end
    
    %% Writing Deformation Gradient Norm    

    fprintf(fID, 'F_norm 1\t%d\tdouble\n', sum(output_cells));
    F_norm = zeros(sum(output_cells), 1);
    count = 1;
    for m = 1:length(C{frame})
        if output_cells(m)
            F_norm(count) = norm(F_cells(:,:,m,frame));
            count = count + 1;
        end
    end
    if strcmp(fileFormat, 'ASCII')
        fprintf(fID, '%f\n', F_norm.');
    else
        fwrite(fID, F_norm.', 'double', 'b');
    end
    
    %% Writing Deformation Gradient
    
    for ii = 1:3
        fprintf(fID, 'F_%d_XYZ 3\t%d\tdouble\n', ii, sum(output_cells));
        temp = squeeze(F_cells(ii,:,output_cells,frame));
%             temp(isnan(temp)) = 0;

        if strcmp(fileFormat, 'ASCII')
            fprintf(fID, '%f\n', temp);
        else
            fwrite(fID, temp, 'double', 'b');
        end
    end
    
    %% Writing Strain
    
    for ii = 1:3
        fprintf(fID, 'Strain_%d_XYZ 3\t%d\tdouble\n', ii, sum(output_cells));
        temp = squeeze(Strain_cells(ii,:,output_cells,frame));
%             temp(isnan(temp)) = 0;

        if strcmp(fileFormat, 'ASCII')
            fprintf(fID, '%f\n', temp);
        else
            fwrite(fID, temp, 'double', 'b');
        end
    end
    
    %% Writing Stretches
    
    fprintf(fID, 'Stretch_1 1\t%d\tdouble\n', sum(output_cells)); 
    temp = squeeze(Stretch_cells(1,output_cells,frame));
    if strcmp(fileFormat, 'ASCII')
        fprintf(fID, '%f\n', temp);
    else
        fwrite(fID, temp, 'double', 'b');
    end
    
    fprintf(fID, 'Stretch_2 1\t%d\tdouble\n', sum(output_cells)); 
    temp = squeeze(Stretch_cells(2,output_cells,frame));
    if strcmp(fileFormat, 'ASCII')
        fprintf(fID, '%f\n', temp);
    else
        fwrite(fID, temp, 'double', 'b');
    end
    
    fprintf(fID, 'Stretch_3 1\t%d\tdouble\n', sum(output_cells)); 
    temp = squeeze(Stretch_cells(3,output_cells,frame));
    if strcmp(fileFormat, 'ASCII')
        fprintf(fID, '%f\n', temp);
    else
        fwrite(fID, temp, 'double', 'b');
    end
    
    %% Writing Eigenvectors
    
    for ii = 1:3
        fprintf(fID, 'Eigenvector_%d 3\t%d\tdouble\n', ii, sum(output_cells)); 
        temp = squeeze(Eigenvector_cells(:,ii,output_cells,frame));
%             temp(isnan(temp)) = 0;

        if strcmp(fileFormat, 'ASCII')
            fprintf(fID, '%f\n', temp);
        else
            fwrite(fID, temp, 'double', 'b');
        end
    end

    
	fclose(fID);
     
	fprintf('Done!\n');
end

end

