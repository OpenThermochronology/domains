---
### ABOUT domains

Inversion code to determine diffusion-domain structure of a multi diffusion domain sample (usually K-feldspar).

This program is basically Oscar Lovera's autoarr program with an updated console user interface plus integration of
plotting using gmt (Generic Mapping Tools). I made a few minor changes; blame any issues on me, not Oscar.

---
### AUTHOR

Peter Zeitler, Lehigh University, Bethlehem, PA USA

---
### COMPILATION

Compile `domains` like this (you MUST do it this way, using the `fno-automatic` and the `fallow-argument-mismatch` flags):

`gfortran domains-130.f90 -o domainsM2 -fno-automatic -O2  -fallow-argument-mismatch -w; rm -f *.mod`

(substitute your local source file name for 'domains-130.f90', and your preferred executable name for 'domainsM2').

Using O2 optimization seems to work reliably.

The two '-f...' flags are required to make this complex legacy code compile: newer compilers are more strict about F77
code that used to generate warnings. These now throw errors that block compilation, so the two flags are required to work around that.

The second command (`rm...`) just cleans up a compile-time file that is not needed following compilation.

For MacOS, a good source for gcc and gfortran installer packages can be found at:

[hpc.sourceforge.net](https://hpc.sourceforge.net)

---
### REQUIREMENTS FOR PLOTTING

Although domains can be used without the plotting option, it's not advisable and it will be tedious.

For the plotting option to work, you need to have an installation of either gmt 5 or gmt 6. Installation packages and instructions can be found at:

[www.generic-mapping-tools.org/download/](https://www.generic-mapping-tools.org/download/)

---
### USAGE

`./domains SAMPLENAME SIZEFILENAME PLOTFLAG`

- SAMPLENAME is a string of 10 or fewer characters (extra will be ignored) to name your output

- SIZEFILENAME is the name of a file in your current working directory that is in .size format (see below)

- PLOTFLAG is an integer: 1 results in use of gmt plotting commands to display output (gmt must be installed)

*NOTE: You can place place a copy of your domains executable in any directory in your PATH, and then
 do work with this code in any other directory containing data (i.e., you won't need to drag around copies of the executable).
 In this case, in any directory containing your datafile, USAGE becomes simply:*

`domains SAMPLENAME SIZEFILENAME PLOTFLAG`


---
### HELP AND EXAMPLES

Consult the quick-start parts of the guide for an example of how to run `domains`. There is also an EXAMPLES directory including some sample input files and examples of the output.

---
