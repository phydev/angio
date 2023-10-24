## Angiogenesis Project
This project is hosted by the [CFisUC](http://cfisuc.fis.uc.pt/) at the [University of Coimbra](www.uc.pt) 
and consists in a phase-field model for tumor angiogenesis. The model is based on the paper published by
[Travasso et al. (2011)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0019989) with some physical and computational improvements.

### Version 6.0.s:
- Coupled the blood flow and the hypoxic cell deactivation;
- The hypoxic cells have volume;
- The system is modeled in three dimensions;
- tools: Routines to measure the number of branches, anastomoses and the vessels diameter;

### Running with Docker:

1. Build the image:
```docker build --platform linux/x86_64 -t angio docker/``` (the ```--platform``` flag is necessary for the image to run on Apple Silicon)
2. Run the container:
```docker run -it -v $(pwd):/code angio /bin/bash run.sh <run_id>```. Please make sure that an input file named ```inp<run_id>``` is present in the root directory of the project.

### Publications:
M. Moreira-Soares, R. Coimbra, L. Rebelo, J. Carvalho & R. D. M. Travasso. [Angiogenic Factors produced by Hypoxic Cells are a leading driver of Anastomoses in Sprouting Angiogenesis–a computational study](https://www.nature.com/articles/s41598-018-27034-8). Scientific Reports 8, 8726 (2018)


![alt tag](https://moreirasm.files.wordpress.com/2018/06/300ms.gif?w=364&h=335)
