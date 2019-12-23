# Code for ["On the numerical solution of the exact factorization equations"](https://aip.scitation.org/doi/10.1063/1.5090802)

This code includes all the config files and plotting scripts to replicate the figures and results in the EF paper linked above. 

### Compiling
This code requires C++11 and OpenMP to be compiled. A makefile is provided so simply navigate to the root directory and type

```
make
```

and the main solver executable will be compiled. (The crude set of tests is compiled separately)

### Running the code
To replicate these results simply run 

```
./paper_runs.sh
```

which will loop through the different config files, with various things switched on or off. 


### Making plots
Assuming all these run without error the  following script in ShinMetiui/plots/ can be run

```
./do_plots.sh
```

which will first run the Gnuplot code to generate a bunch of TeX and then convert it to PDF. 

The folder ShinMetiu/ext_exact_tdnd contains the exact nuclear density, solving the full 2D TDSE, for this model from a separate code and this is what is is used in the plotting comparisons.

### Further info/instruction/discussion will be forthcoming