![build status](https://github.com/phydev/angio/actions/workflows/docker-image.yml/badge.svg)

----

## Angiogenesis Project

This project is hosted by the [CFisUC](http://cfisuc.fis.uc.pt/) at the [University of Coimbra](www.uc.pt) 
and consists in a phase-field model for tumor angiogenesis. The model is based on the paper published by
[Travasso et al. (2011)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0019989) with some physical and computational improvements.

###  What is new on v6.0.s

- Coupled the blood flow and the hypoxic cell deactivation;
- The hypoxic cells have volume;
- The system is modeled in three dimensions;
- tools: Routines to measure the number of branches, anastomoses and the vessels diameter;

### Running with Docker

1. Build the image:
```docker build --platform linux/x86_64 -t angio docker/``` (the ```--platform``` flag is necessary for the image to run on Apple Silicon)
2. Run the container:
```docker run -it -v $(pwd):/code angio /bin/bash run.sh <run_id>```. Please make sure that an input file named ```inp<run_id>``` is present in the root directory of the project.

### Parameters

We provide typical parameters for simulations in [`input_file`](https://github.com/phydev/angio/blob/master/input_file). 
See the description of each paramter in the table bellow and the mapping with experimental data in [Supplementary material](https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-018-27034-8/MediaObjects/41598_2018_27034_MOESM1_ESM.pdf).


| Reference value | Variable               | Description                            |
|-----------------|-----------------------|----------------------------------------|
| 4.00 | cell_radius | Radius of individual cells |
| 100.0 | diffusion_const | Diffusion constant for VEGF |
| 1.00 | interface_width | Width of the phase-field interface |
| 0.30 | vegf_p | VEGF concentration for maximum proliferation |
| 0.09 | vegf_c | VEGF concentration for branching |
| 20.0 | diff_oxy_length | Diffusion length for oxygen |
| 6.25 | vegf_rate | Rate of VEGF uptake by cells |
| 1.00 | vegf_source_conc | Concentration of VEGF at source |
| 1.00 | prolif_rate | Rate of proliferation of endothelial cells |
| 5.00 | vessel_radius | Initial radius of blood vessels |
| 150000 | total_time_step | Total time steps for simulation |
| 0.0010 | dt | Time step size |
| 800.00 | chi_chemiotactic_resp | Chemotactic response of endothelial cells |
| 100, 100, 50,  1, 1, 1 | Lx_Ly_Lz_dx_dy_dz | Simulation domain size and grid spacing |
| -754333222 | random_seed | Seed for random number generation - must be negative |
| 20 | number_of_boundary_points | Number of boundary points to keep track |
| 10000 | source_max | Maximum number of VEGF source points |
| 0.01 | vegf_grad_min | Threshold for VEGF gradient |
| 0.03 | vegf_max | Maximum VEGF concentration |
| 2.00 | depletion_weight | Energy cost to avoid overlap between vessels and hypoxic cells |
| 2000 | output_period | Time period between outputting results |
| 40000 | extra_steps | Additional time steps for simulation |
| 4000 | max_number_of_tip_cells | Maximum number of tip cells allowed |
| T | thinning_FT | Flag for using thinning algorithm |
| F | periodic_FT | Flag for using periodic boundary conditions |
| F | flow_FT | Flag for computing the blood flow |



Beware that the input file `inp001` is solely used to run the docker image with Github actions. The parameters were changed to produce a short simulation with small grid, few iterations, and no sources of VEGF. Please do __not__ base future studies on this file. 

### Postprocessing tools

In the folder `tools` we have scripts to postprocess the data generated with the simulations. The documentation for each script is found [here](https://github.com/phydev/angio/blob/master/tools/about.md).

### Publications

M. Moreira-Soares, R. Coimbra, L. Rebelo, J. Carvalho & R. D. M. Travasso. [Angiogenic Factors produced by Hypoxic Cells are a leading driver of Anastomoses in Sprouting Angiogenesis–a computational study](https://www.nature.com/articles/s41598-018-27034-8). Scientific Reports 8, 8726 (2018)


![alt tag](https://moreirasm.files.wordpress.com/2018/06/300ms.gif?w=364&h=335)

### Acknowledgements

Special thanks to [João Simões](https://github.com/joaoplay) for supporting the long term reproducibility of this code by dockerizing the repository. 

This project was funded by the National Council of Technological and Scientific Development (CNPq - Brazil) under the grant 235101/2014-1. 
