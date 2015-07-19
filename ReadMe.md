#ALLELIC DIVERSITY AND ADAPTIVE POTENTIAL

This repository includes code scripts and parameter files for the simulation analyses showed in the paper: Vilas et. al. **Allelic diversity for neutral markers retains a higher adaptive potential for quantitative traits than expected heterozygosity** currently under revision on *Molecular Ecology* (MEC-15-0617). Please refer to the text in the paper in order to understand the following descriptions.

## Requirements

This scripts and programs were developed for a linux cluster using bash, python, R and C so these tools are required (Bash, R, python > 2.6 and a C compiler). Some of the scripts also require to be run in a machine with a directory /scratch so it would be useful to make it before or modify the appropriate lines in the scripts.

The program selection.c needs to be compiled in advance using the libraries genlib.c ranlib.c com.c and linpack.c.

The third-part programs metapop (Pérez-Figueroa et al. 2008, download: https://dl.dropboxusercontent.com/u/3144341/metapopweb/metapop.14.08.zip Unzip the file and run make. This will create an executable file called metapop) and quantiNemo (v1.0.4, Neuenschwander et al., 2008, download: http://www2.unil.ch/popgen/softwares/quantinemo/) need to be compiled and the binaries added into the $PATH

## Files
These are the descriptions of the different files in this repository, divided into scripts and templates directories.

### Scripts

- premetapop.sh: This script run the first part of the simulations. Using quantiNemo initialize the base population, subdivide it in eigth subpopulations and run these populations for 50 generations. The script takes two arguments, the number of qtls and the number of neutral markers to be simulated and generates all the directories and configuration files according to those parameters.

- metapoping.sh: This script run the second part of the simulations. Using metapop with the data from premetapop.sh obtains the number of individuals of each subpopulation needed to build the synthetic population in function of the desired optimization. It takes four arguments: the number of qtls and the number of neutral markers (these should match a previous run of premetapop.sh), the type of loci subject of optimization (qtl ot ntrl) and the optimization method (following the metapop nomenclature: NA for O_A, He for O_H).

- selection.sh: This script run the third and last part of the simulations. Using synthetic.py creates the synthetic populations given the numbers obtained by metapop in the previous step. Then, using selection.c run those populations for artificial selection. After all these steps, stats for different parameters are gathered and summarized.  It takes four arguments: the number of qtls and the number of neutral markers (these should match a previous run of premetapop.sh), the type of loci subject of optimization (qtl ot ntrl) and the optimization method (following the metapop nomenclature: NA for O_A, He for O_H).

- genAllelicFile.sh and genAllelicMArkerFile.sh: These scripts automatically create the appropriate allelic files used by quantiNemo given the number of qtls and markers specified in premetapop.sh.

- fs2mp.py: This is a converter between Fstat file (used by quantiNemo) and metapop files. It is internally used in premetapop.sh

- synthetic.py: This python script builds the synthetic populations given the genotypes from quantiNemo and the poolfile (number of individuals of each subpopulation) from metapop. It is used in selection.sh

- selection.c: This C program (to be compiled as selection) simulates the evolution of population across several generation of artificial selection to increase an hypothetical quantitative trait determined by the qtls. For each generation yields the phenotipic mean of the population for the trait, as well as heterozygosity and allelic diversity for qtls and neutral markers.

### Templates

There are some template files with some default parameters for quantiNemo, metapop and selection that are automatically tuned in the scripts in order to adapt them to the specific simulations. These templates include some SET_* words which would be substituted automatically by the scripts to the right values.

## Procedure

In order to run any of the cases described in the paper you need a working directory with both the directories in this repository (scripts and templates). Then the scripts should be called from the working directory with the appropriate parameters. New directories will be created automatically by the scripts. To create the base populations (for 10 qtl and 100 neutral markers, by instance):
```
scripts/premetapop.sh 10 100
```

Then, to obtain the contribution of the subdivided populations to obtain a synthetic population using both methods of optimization over the two kind of loci:
```
scripts/metapoping.sh 10 100 ntrl NA
scripts/metapoping.sh 10 100 ntrl He
scripts/metapoping.sh 10 100 qtl NA
scripts/metapoping.sh 10 100 qtl He
```
After the synthetic populations are defined tou can run the last script to generate them and run the replicates of selection:
```
scripts/selection.sh 10 100 ntrl NA
scripts/selection.sh 10 100 ntrl He
scripts/selection.sh 10 100 qtl NA
scripts/selection.sh 10 100 qtl He
```
Results from selection for this case would be in the directory 10/100/stats with three files (data, averages, errors) for each type of loci and optimization.

To run other case just substitute the number of qtl and neutral markers. This is true for the described as scenario A. For other scenarios scripts should be modified in order to fit the SIZE_SUBPOP variable values accordingly to the values in the scenario.


## Contact

All questions and inquiries about the code included in this repository should be addressed to:
Andrés Pérez-Figueroa (anpefi@uvigo.es)

