# Exploring properties of structural and functional connectomes

In recent years, understanding how the brain functions has been a focus of intense study. The connectome, or the network of connections between cells, regions, and pathways in the brain, is one lens through which the brain is studied. Mapping the connectome in great detail is a challenging task due to the complexity and sensitivity of the brain. The connectome can be broken down into structural and functional connectivity. Structural connectivity refers to the physical connections between different regions of the brain, such as axon fiber bundles. Functional connectivity refers to the temporal interdependence of neural activity in cortical regions or cell populations. Both structural and functional connectivity can be studied using different sub-modalities of magnetic resonance imaging (MRI) scans. In this project, the team aims to use these different MRI scans to model the relationship between structural and functional connectivity maps.

# Core Task 1: Graph properties of structural and functional connectomes

To build an understanding of the full structural connectivity map there are a few steps required to obtain a detailed graphical representation of the brain’s structure. I use high resolution structural MRI images, parcelate the cortex into anatomically sub-regions, and then use tractography to assemble edges between each pair of regions into an abstract graph. In this project I will be using a dataset with a set of multi-modal imaging data from a single healthy adult and structural connectomes derived using different anisotropy thresholds. The pre-processing of the MRI imaging was done by the TractoR software package. The most important files used to reconstruct a structural connectome are listed.






Exploration of the relationship between structural and functional connectomes by using different linear models to predict functional connectivity weights, fij between the regions i and j from the structural connectivity weights, sij
