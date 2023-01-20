# Exploring properties of structural and functional connectomes

In recent years, understanding how the brain functions has been a focus of intense study. The connectome, or the network of connections between cells, regions, and pathways in the brain, is one lens through which the brain is studied. Mapping the connectome in great detail is a challenging task due to the complexity and sensitivity of the brain. The connectome can be broken down into structural and functional connectivity. Structural connectivity refers to the physical connections between different regions of the brain, such as axon fiber bundles. Functional connectivity refers to the temporal interdependence of neural activity in cortical regions or cell populations. Both structural and functional connectivity can be studied using different sub-modalities of magnetic resonance imaging (MRI) scans. In this project, the team aims to use these different MRI scans to model the relationship between structural and functional connectivity maps.

# Core Task 1: Graph properties of structural and functional connectomes

To build an understanding of the full structural connectivity map there are a few steps required to obtain a detailed graphical representation of the brain’s structure. I use high resolution structural MRI images, parcelate the cortex into anatomically sub-regions, and then use tractography to assemble edges between each pair of regions into an abstract graph. In this project I will be using a dataset with a set of multi-modal imaging data from a single healthy adult and structural connectomes derived using different anisotropy thresholds. The pre-processing of the MRI imaging was done by the TractoR software package. The most important files used to reconstruct a structural connectome are listed.

```
- structural/refT1.nii.gz: T1-weighted structural image
- structural/parcellation.nii.gz: A cortical parcellation, where each voxel is labelled with a cortical region (int)
- structural/parcellation.lut: A table of information about the labelled regions in the parcellation
- diffusion/data.nii.gz: The dMRI dataset, a 4D volume with gradient direction varying along the fourth dimension
- diffusion/dti FA.nii.gz: A map of fractional anisotropy
- diffusion/parcellation.nii.gz: The same parcellation as above but translated for the diffusion data
- functional/data.nii.gz: The rs-fMRI dataset, a 4D volume with time varying along the fourth dimension to align the individual volumes and apply temporal filtering 
- functional/parcellation.nii.gz: The same parcellation as above, converted for the functional data
```
# Core Task 2: Modelling the relationships between connectomes

Next, the relationship between structural and functional connectomes is explored by using various linear models to predict functional connectivity weights, fij between regions i and j from the structural connectivity weights, sij. The sparse nature of the structural connectivity matrix is taken into account by considering indirect structural connectivity weights, tij, which represent connections between two regions that pass through one middle region first. These connections are computed by using the equation tij = maxk(min(sik, skj)), where sij, skj ̸= 0. The linear models are applied to the data with and without these indirect connections.

$$f_{ij} = \alpha_{ij} + \beta_{ij}s_{ij}$$

$$f_{ij} = \alpha_{ij} + \beta_{ij}s_{ij} + \gamma_{ij} s_{ij}^{2}$$

$$f_{ij} = \alpha_{ij} + \beta_{ij}+ t_{ij}$$

$$f_{ij} = \alpha_{ij} + \beta_{ij}+ t_{ij} + \gamma_{ij} t_{ij}^{2}$$

$$f_{ij} = \alpha_{ij} + \beta_{ij}+ t_{ij} + \gamma_{ij} t_{ij}$$

# Core Task 3: An Advance Task 

To further understand the complexity of the interconnections between different regions of the brain, an analysis of the structural and functional connectivity graphs is conducted, with a focus on the density of each vertex.

To model the relationship between structural and functional connectivity maps, LASSO is used as the model due to the sparse nature of the structural connectivity graph. LASSO is a linear regression that uses a shrinking parameter, λ, to penalize the loss function and encourage simpler models with small weights. The aim is to relate each functional connection to a subset of structural connections. The method is implemented using a 3-fold cross validation and iterating over every functional connection with a lasso regression. Where, the lowest λ in the array will act closer to an OLS approach. The average of the SSE across all voxels is computed to be 0.07848 over all 3 folds, which is low and suggests that the LASSO approach is a good performing model.
