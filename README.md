# Clouds
# Clouds.jl

A package to **analyze 2D and 3D distributions (clouds)** and treat them as discrete images. A cloud is segmented into cells (**pixel or voxels**), for each **cell** the local gradient, laplacian and curvatures are computed. Cells are grouped into Nodes accordingly with the gradient. The center of the **node** is the local cell with the highest content. Adjacent nodes are grouped into a **graph**.

See example of usage in notebook/clouds_example.jl Pluto Notebook or its [static version](./notebooks/clouds_example.jl.html).
