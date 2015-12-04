## remove_space.sh
This script is supposed to be runned inside the data folder in order to remove all spaces in the files names.
Caution: you must use it only **ONE** time in each folder, if you do that for the second time you will loose your data.

## run_get.F90 and get_data.F90
In these two files you can find the routines to obtain the diameter and number of branches in function of time. You can compile executing the shell script compile_run_get.sh. 


## Running the run_get
You should give to the program an input file (e.g. inp_get):
###### ./run_get < inp_get

The program will create some output and auxiliar files (.aux). The importan files are:
###### diameterXXX.pdf   | graphic of the diameter in function of time
###### branchesXXX.pdf   | graphic of the number of branches in function of time
###### dataXXX.dat       | file containing the number of branches and diameter



