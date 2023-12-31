# Subgraph Strain

This code in this repository will reproduce the analysis in "Computing Whole Embryo Strain Maps During Gastrulation".

## Setup

### Software Requirements:
- MATLAB R2012a or newer
- MATLAB image processing toolbox addon

### System Requirements:
- At least 8 GB of RAM

Download/Clone this repository into a location of your choosing. MATLAB version R2022a (tested on Windows 10 and Ubuntu 22.04) was used to develop this code, but versions newer than R2012a should work as well. Additionally, the image processing toolbox addon for MATLAB was used for the function 'regionprops'.

While not necessary to run the code, for visualization we recommend the tool Paraview, which can be downloaded here: [https://www.paraview.org/download/](https://www.paraview.org/download/). 

## Prerequisite data

This project utilizes the published segmentation data found in "[Deconstructing gastrulation at single-cell resolution](https://www.sciencedirect.com/science/article/pii/S0960982222003293?via%3Dihub)". Their published data and code can be downloaded here: [https://doi.org/10.6084/m9.figshare.18551420.v2](https://doi.org/10.6084/m9.figshare.18551420.v2).

The files 'Deconstructing Gastrulation - Data.zip' and 'Deconstructing Gastrulation - Code.zip' should be extracted into the ./data directory in this repository.

The folder structure before running the code should be: 
```
├── data
│   ├── Img_1620 (intercalations)
│   │   ├── ...
│   ├── Img_1830 (divisions)
│   │   ├── ...
│   ├── DG_calc_cell_areas.m
│   ├── ...
│   ├── full_lineage_preprocess.mat
├── output_vtk
│   ├── ...
├── ...
└── README.md
```

## Running

In MATLAB, the script 'process_movie.m' will load the segmented data and compute the deformation gradient and strain for each cell in the embryo. The results will be saved to a file called 'output.mat' which contains:

```F_cells``` - the deformation gradient tensor for each cell at each time point

```cell_frequency``` - the number of samples containing a particular cell at each time point

```max_stde``` - the maximum standard error of each cell at each frame

```inview_cells``` - a logical array indicating if the cell is "in view". Cells which are not "in view" have invaginated, divided, or are future daughter cells.

```subgraph_size``` - the size of each subgraph

```samples_per_frame``` - the number of subgraph samples per frame

```batch_size``` - the number of subgraphs to sample in batch, which speeds up the sampling process

```dropout``` - the fraction of subgraphs to forget after each frame

The Green strain can be computed using the function 'compute_strain_stretch.m'. This function also can return the right and left Cauchy-Green deformation tensors, the principal stretches, and their corresponding eigenvectors. 'compute_strain_stretch.m' is called at the end of 'process_movie.m'.

By default the parameters ```samples_per_frame``` and ```batch_size``` are set to 50,000 and 1,000 respectively. Once you have tested to make sure the code is running ```samples_per_frame``` should be increased to 100,000 or even 1,000,000 and ```batch_size``` should be increased to 5,000 or 10,000.

Typical run time for ```samples_per_frame = 100,000``` should be around one hour, depending on your device.

## Visualization

'process_movie.m' will export a sequence of .vtk files which can be loaded by Paraview. Each vtk file contains:

```Cell_IDs``` - the id number for each cell

```Stretch_[1-3]``` and ```Eigenvector_[1-3]``` - the principal stretches and their corresponding eigenvectors. Note: Stretch_3 is the largest magnitude stretch and Stretch_1 is the smallest.

```F_[1-3]_XYZ``` and ```Strain_[1-3]_XYZ``` - the rows of the deformation gradient and strain tensors are exported as vectors for use in Paraview. For example the $E_{XX}$ component would be the X component of Strain_1_XYZ and $E_{YZ}$ would be the Z component of Strain_2_XYZ.

```Samples_per_cell``` - the number of subgraph samples per cell in each frame

```F_norm``` - the L2-norm of the deformation gradient for each cell

```Max_STDe``` - the maximum standard error for each cell in each frame



