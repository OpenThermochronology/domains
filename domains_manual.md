---
### CONTENTS

This document, about the domains inverse code, serve three purposes.

* A quick start guide that will help you run one of the sample input files (sn17.size)
* A more detailed description of some code options
* Some more detailed background about this code and how it relates to Oscar Lovera's original autoarr code.

Note: A good understanding of the MDD diffusion model is required for use of this code and associated notes.

#### Sections
1. CONTENTS  
2. QUICK START GUIDES AND NOTES
3. CODE INPUT AND OUTPUTS
4. PLOT ELEMENTS
5. OTHER USER OPTIONS
6. MORE INFORMATION ABOUT THE CODE AND HOW IT OPERATES

---
### QUICK START (simple auto-run)

1. Navigate to a directory containing an in put file in .size format

2. Run domains: `domainsM2 test sn17.size 1`

3. `Choose diffusion geometry (0 = slab, 1 = sphere):` **1**

4. `Enter upper temperature cut-off for modeling (˚C):` **1050**  
*For feldspars, enter temperature above which you think the sample began to partially melt*

5. `To override defaults, enter "y". Otherwise type "n" to proceed:` **n**

6. `To select E & Ro reference, enter "y". Or enter "n" to proceed:` **n**

7. ` Model finished with no major worries.` (a kinetic plot should appear!)

### QUICK START (choosing regression for activation energy)

1. Navigate to a directory containing an in put file in .size format

2. Run domains: `domainsM2 test sn17.size 1`

3. `Choose diffusion geometry (0 = slab, 1 = sphere):` **1**

4. `Enter upper temperature cut-off for modeling (˚C):` **1050**

5. `To override defaults, enter "y". Otherwise type "n" to proceed:` **n**

6. `To select E & Ro reference, enter "y". Or enter "n" to proceed:` **y**

7. ` Select an option for activation energy determination:` **y**

8. `To select E & Ro reference, enter "y". Or enter "n" to proceed:`
  * `0 -- use the values chosen by this program` (same as auto-run)
  * `1 -- specify a range of heating steps to regress` (only use if you have seen plot!)
  * `2 -- select from a table of values`
  
  `Enter option:` **2**

9. `You've elected to choose E and the Ro reference from a table of models:`

  `Enter minimum steps for regression:` **7**
  
10. Two tables appear, sorted by activation energy and by fit. Select from them.
  
  `Enter preferred model:` **70**

11. ` Model finished with no major worries.` (a kinetic plot should appear!)

### QUICK START (choosing regression based on set of contiguous data)

1. Navigate to a directory containing an in put file in .size format

2. Run domains: `domainsM2 test sn17.size 1`

3. `Choose diffusion geometry (0 = slab, 1 = sphere):` **1**

4. `Enter upper temperature cut-off for modeling (˚C):` **1050**

5. `To override defaults, enter "y". Otherwise type "n" to proceed:` **n**

6. `To select E & Ro reference, enter "y". Or enter "n" to proceed:` **y**

7. ` Select an option for activation energy determination:` **y**

8. `To select E & Ro reference, enter "y". Or enter "n" to proceed:`
  * `0 -- use the values chosen by this program` (same as auto-run)
  * `1 -- specify a range of heating steps to regress` (only use if you have seen plot!)
  * `2 -- select from a table of values`
  
  `Enter option:` **1**

9. `You've chosen to regress a single set of contiguous points.`

  `Enter step number, start of regression interval:` **3**  
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`end of regression interval:` **14**  
  
10. ` Model finished with no major worries.` (a kinetic plot should appear!)

### QUICK START -- NOTES

1. Most models are quick, but a few can be slow, up to minutes. You should see a report of model progress.

2. If several minutes have gone by or you do not see progress notes, the inversion might be stuck. Or, you might see a crash or an error report.

3. If this happens on the first attempted run, woe is you. Examine your input data for problems or unusual features.

4. If the code hangs after a run following an earlier one from which you've managed to get a plot, look for odd features in the data like convex-upward arrhenius trends, an outlier, or an odd first point, then try regressing using different combinations numbers of steps.

5. See below for information about other options available, some of which are hidden behind other options. The nature of the code made it necessary to have this flow; simplifying this would require a major rewrite.

---
### CODE INPUT AND OUTPUTS


#### INPUT FILE

- Remember that your input file must be in UNIX file format. Data format:

- number-of-modeled-steps (note: one less than total since we can't deal with f=1.000)

temp(degC)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fractional loss&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;time(minutes)
&nbsp;&nbsp;&nbsp;*(as many lines as there are modeled heatings steps)*

*It's up to you but a useful convention is to use `.size` as the file suffix for you input files.*

---
#### OUTPUT FILES

Various text output files record both the observed data and model results, and run information. Additionally, the file `domains-SAMPLENAME.in` is in the proper format to use with the arvert inversion (you just need to rename the file as `domains.in` to use with arvert).

All results are placed in a subdirectory named `results.SAMPLENAME`. By intention, this subdirectory is **overwritten** if you rerun a model with the same name. My reasoning is that when assessing a sample, you will usually make several model runs to hone in on a preferred model, so it seems better not to pollute the world with more directories. If you want to retain intermediate model results, you will either need to rename the output subdirectory or run models with different sample names.

*XXX will be replaced by samplename*
| OUTPUT FILES  | UNIT | Description |
| ------------- |:-------------:| ------------- |  
| domains-XXX.in      | 14  | input file for Arvert inversion |
| report-XXX.out      | 28  | summarizes model results |
| ener-XXX.out      | 30  | Ea values determined during search for best value |
| arr-model-XXX.dat      | 16  | model arrhenius data, for X-Y plotting |
| logr-model-XXX.dat      | 16  | model logRRo data, for line plotting as spectrum |
|  arr-observed-XXX.dat      | 18  | observed arrhenius data, for X-Y plotting |
|  logr-observed-XXX.dat      | 18  | observed logRRo data, for line plotting as spectrum |
| logr-XXX.dat      | 20  | observed logRRo data, as simple list for autocorrelation |
| kinetics-XXX.pdf      | N/A  | plotted kinetic and logRRo data, plus run info |


---
### PLOT ELEMENTS

The plot file `kinetics-SAMPLENAME.pdf` summarizes a great deal of information about the sample and the model run. Here's a description of the elements on the plot.

*The top plot* is an arrhenius plot with the observed data shown as blue circles and the model predictions as red stars. The faint gray lines show the arrhenius trends of each individual diffusion domain; the darker thicker line shows the r<sub>o</sub> reference line.

*The bottom plot* is a log(r/r<sub>o</sub>) plot showing the relative retentivity of the sample as a function of gas release (assumed here to be <sup>39</sup>Ar). Blue lines show the observed data and red lines the model predictions. The right-hand axis gives values of log<sub>10</sub>(r/r<sub>o</sub>); the left-hand axis expresses these values as sizes (radii) relative to the r<sub>o</sub> reference. The thick gray bars illustrate the domain structure, showing volume fraction and relative size, but it is important to realize that in reality all domains contribute some gas release during heating - using the <sup>39</sup>Ar-loss axis is just a plotting contrivance. The pinkish horizontal bar shows the region of fractional loss across which the data were **not** modeled. Finally, at the lower right of this plot, the fit of the data across the modeled interval are reported (I believe this is the unreduced chi square value).

*The tables at upper right* report summary run information, and a summary of the domain structure including each domain's retentivity expressed as a closure temperature for a 10˚C/m.y. cooling rate.

---
### OTHER USER OPTIONS

This section is meant to provide a guide to the other options available to the user. It's organized by the program flow. If an option in the sequence is not mentioned, it's part of the standard flow covered above in the quick-start guides.

`To override defaults, enter "y". Otherwise type "n" to proceed:` **y** &nbsp;&nbsp;*Note that once you choose to override defaults you must enter values for each option*  
- `Enter weighting option for Arrhenius-plot regression:` - 0 &nbsp;&nbsp;*applies no weighting (default; uncertainties will be good); 1 - chooses the legacy weighting by step size (uncertainties not reliable)*
- `Maximum number of domains (<= 15) (Default is 10):` &nbsp;&nbsp;*use this if you have reason to allow more domains than the default, or have reason to limit the total to fewer than the default*
- `Minimum number of domains (>= 3) (Default is 3):` &nbsp;&nbsp;*self explanatory*
- `To keep Do fixed type 1, otherwise enter 0:` &nbsp;&nbsp;*0 is the default; if you choose 1, later in the input flow you'll be asked for a fixed maximum value of Do (input is separated because of the nature of the code, sorry!). After you've chosen a regression for E, you will be asked:*
	- `Enter maximum value for log(r/ro) (default=free (0)):` &nbsp;&nbsp;*Enter a value for the maximum value of log<sub>10</sub>(D<sub>o</sub>/r<sub>o</sub><sup>2</sup>) you want to see for the largest domain (usually something like 2 or 3) - this avoids situations where the last bits of data cause the code to shoot up to a silly large value at the end.*

### MORE INFORMATION ABOUT THE CODE AND HOW IT OPERATES

#### Algorithm
This program inverts for a domain-size distribution by minimizing the misfit between predicted and observed logR/Ro data. After selection of an Ro reference line with associated activation energy, the code guesses a starting distribution which is then modified using the Levenberg-Marquardt method, a form of nonlinear least-squares fitting.
The core routines used in autoarr (mrqmin() and mtqcof()) were taken from the Numerical Recipes volume by Press et al.

---
#### Major and minor modifications by Zeitler

Here are most of the changes I made to the original Lovera code. The following includes some technical description of issues you might want to be aware of, especially if you are thinking about how to propagate uncertainty in the domain structure into *tT* inversions.

 *  added command-line inputs that identify model run, input file, and plotting option  
	- samplename is used to label output files and create a results-SAMPLENAME output directory 	
	- added an option to have the code produce plots as a final step (this uses system() calls to gmt, so the user must have gmt installed).

 * switch to single input file having the format: TempC floss delta-time(m)

 * code reformatted to compile under f90 and prettied up (alignment, white space, comments, no loop labels)

 * expanded the number of steps over which the auto-code looks for a maximum E early in the gas release (originally was 20 steps; now it is 30 to better accommodate more detailed heating schedules that might include many isothermal replicates)

 * expanded options for determining E:
      - use legacy code's auto-suggestion routine, OR
     - let user specify the regression range, in terms of contiguous steps, OR
     - let user inspect a table showing highest-E and best-fit higher-E regressions before selecting a model

---
**Extended discussion about selection of E and its uncertainty.** There are three significant issues in the way the original autoarr worked to choose E. You **must** understand this
in order to get consistent results that will match those from other software and from prior and future MDD work ythat ou do.

> *First,* the original routine in Oscar's autoarr code had several criteria for finding the rollover on the Arrhenius plot. Once it does, the code _averages_ the last and next-to-last values it determined for E and the y-intercept. This is ok as there are many approaches one could advocate in finding the "best" E, but in practice this recipe means that you will never get the same E and intercept as you would using other approaches that simply regress sequences of steps. This is not a huge deal, but it will complicate comparison between methods. Usually the difference is not large, but just be aware of this. 

> *Second,* the original auto-routine always starts with the first step -- it does not allow results from tiny and spurious early steps to be omitted, other than editing the input file. I prefer to keep the input file as a pristine record of the observations, and not be forced to maintain and track two separate files.

> *Third,* the original autoarr approach weights the input data for both the line-fitting and also the non-linear least-squares routine for domain-finding: the regression is NOT a simple linear regression of Y on X. What I am calling the "legacy" weights are just the step size. I am not sure what went into this decision. Most generally, one could argue that while  small early steps might have more error, these key early steps determine the Ro reference and address a substantial fraction of the sample's temperature-constraining power, so I do not think that de-emphasizing them is a good idea.

> Simple regressions of Y on X tend to be biased to their component data's larger values. Note also that the regression we are doing to get E and Do are in logD- 1/K space, so any weighting using f-loss will be non-linear with respect to the actual regressed data. Traditionally, one would weight using uncertainties (squared inverse of the variance), but given how conventional step heating is done, that would require additional calculations be made by mass-spec data-reduction software to propagate beam and other errors into the cumulative f determination. One could be cheap and just do a generic "typical" determination of these but that's not rigorous at all. In either case the f-loss errors would have to be propagated through the logD calculation.

> What to do?

> To recap, the legacy code worked as follows: the weighting by step size was used for both the linear regression in the Arrhenius plot and the non-linear least-squares routine used to find the domain structure. As a consequence, the uncertainty data returned by the Fit() function would not be meaningful.

> **_I changed the default behavior of the code__** and added an option in the "advanced" input section that allows the user to make choices about line-fitting. The user now has the option to use

> * *an unweighted regression to get E (**DEFAULT**), OR
> * the Lovera legacy weights, OR
> * *a weight proportional to the 1/logD value

> The default weighting produces a simple linear regression of y on x, minimizing the y misfit. True propagation of uncertainties would require that a two-error regression be used that accounts for both temperature and f-loss uncertainties. The legacy weighting is just by step size.

> There are two important things to note. First, no matter what weighting is used for determining E, the stepsize weights are passed to the non-linear least-squares routine, for compatibility with the legacy code.  Second, if you use the legacy stepsize weighting for E, it is not possible to get uncertainty estimates because we are not supplying uncertainties and the weights make no statistical sense. Using the unweighted regression, it is possible get error estimates using a bootstrap approach that uses the scatter around the best-fit line to guess at the experimental uncertainty. When you do this, you can get estimates for the slope and intercept error, but any attempts at quantifying goodness of fit for the line would be meaningless (impossible).

---
#### Some information about variables used in the code

If you start to examine the code, this table might help you get started. This is not at all a complete list!

| Variable  | Description |
| ------------- |:-------------:|
| samplename      | text string (10 char or less) used to name files and directories (obtained from command line)     |
| gmt      | flag for using gmt for  plotting (1) or not (any other integer) (obtained from command line)     |
| e      | activation energy (kcal/mol)     |
| ord      | log<sub>10</sub>(D<sub>o</sub>/r<sub>o</sub><sup>2</sup>)      |
| c(j)      | volume fraction of jth domain     |
| rp(j)      | plateau size (log(r/r<sub>o</sub>) plot)     |
| ni      | number of lab heating steps modeled     |
| telab()      | heating step temperature (K) (read from file as ˚C)     |
| tilab()      | heating step duration (minutes)     |
| tinv      | inverse laboratory temperature (10 000 / K)     |
| dzx      | -log(D<sub>o</sub>/r<sup>2</sup>) (1/sec)  |
| f(j)*100      | cumulative percent <sup>39</sup>Ar released     |
| avt      | (f(j)+f(j-1))*50     |
| xlogr      | log(r/r<sub>o</sub>)     |
| rad(j)      | size of the jth diffusion domain     |
| c(j)      | volume fraction of jth diffusion domain     |
