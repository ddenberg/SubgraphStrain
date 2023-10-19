function F_cells = batch_integrate_inc_def_gradient(incF_cells)
num_frames = size(incF_cells, 4);

F_cells = incF_cells;
for t = 2:num_frames
    F_cells(:,:,:,t) = pagemtimes(incF_cells(:,:,:,t), F_cells(:,:,:,t-1));
end
end