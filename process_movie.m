
data_dir = './data';
addpath(data_dir);

% load segmentation data
[C, E, V] = DG_load_segmented_embryo(fullfile(data_dir, 'Img_1620 (intercalations)/Mesh')); % raw mesh
load(fullfile(data_dir, 'full_lineage_preprocess.mat')); % cell lineage and centroids

% compute centroids for each cell
C = DG_calc_cell_centroids(C, V);

% reformat data structure. Remove unused fields in C and make entries in V arrays
C = cellfun(@(s) rmfield(s, {'edges', 'cells'}), C, 'UniformOutput', false);
V = vertices_to_array(V);

% create graphs at each time point using cell-cell adjacencies
[G, inview_cells] = build_graph(C, E);

% center V and cell_lineage_centroids based on first time point
mean_xyz0 = mean(V{1}, 1, 'omitnan');
V = cellfun(@(A) A - mean_xyz0, V, 'UniformOutput', false);
cell_lineage_centroids = cell_lineage_centroids - reshape(mean_xyz0, [3, 1, 1]);

% final frame is the length of the cells cell array
final_frame = length(C);

% create list of cell ids based on their order in the C cell array
cell_ids = (1:size(inview_cells, 1)).';

% number of subgraphs to sample per frame
samples_per_frame = 5e4; 

% number of subgraphs to sample in each batch. Smaller values will be slower but compute standard error estimates more often
% larger values can be faster but may over sample some cells
batch_size = 1e3; 

% what fraction of subgraphs to randomly discard from the previous frame pair. Subgraphs which are not
% connected will be discarded on top of the dropout fraction
dropout = 0.2; % (dropout of 1 means all subgraphs are re-sampled)

% size of each subgraph
subgraph_size = 20;

tic;
% here we sample subgraphs, compute their incremental deformation, and average subgraph
% incremental deformation to cells
[incF_cells, cell_frequency, max_stde] = ...
    sample_incremental_deformation(G, cell_lineage_centroids, inview_cells, lineage, ...
                                   final_frame, subgraph_size, samples_per_frame, batch_size, dropout);
toc;

% cell incremental deformation is integrated over time 
F_cells = batch_integrate_inc_def_gradient(incF_cells);

% using the deformation gradient we can compute the strain tensors, stretch, and stretch eigenvectors
[Strain_cells, RightCauchy_cells, LeftCauchy_cells, Stretch_cells, Eigenvector_cells] = compute_strain_stretch(F_cells, true);

% save output
% you can pick and choose what you want to save as output. Strain (and other derived
% quantities) can be computed from F_cells relatively quickly so it doesn't have to be saved.
% It is relatively cheap to re-compute later
save('output.mat', 'F_cells', 'cell_frequency', 'max_stde', 'inview_cells', 'subgraph_size', ...
    'samples_per_frame', 'batch_size', 'dropout');

% export a sequence of vtk files to be viewed in paraview later
export_vtk('output_vtk', 'BINARY', C, V, inview_cells, cell_frequency, F_cells, ...
    Strain_cells, Stretch_cells, Eigenvector_cells, max_stde, cell_ids, 1:final_frame);


