# Tools

## remove_space.sh
This script is supposed to run inside the data folder to remove all spaces in the files names.

**Caution**: you must use it only **ONCE** for each folder, if you run the second time you will lose your data.

## run_get.F90 and get_data.F90
These two files provide routines to compute the diameter and number of branches as a function of time. Compilation is performed with the shell script `compile_run_get.sh`. 


## Running the run_get

You should give to the program an input file (e.g. `inp_get`):

`./run_get < inp_get`

The program will create some output and auxiliar files (.aux). The important files are:

| file | description |
| --- | --- |
| diameterXXX.pdf   | graphic of the diameter in function of time |
| branchesXXX.pdf   | graphic of the number of branches in function of time |
| dataXXX.dat       | file containing the number of branches and diameter |



