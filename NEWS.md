
# FMM 0.4.1

Enhancements:

- Improved convergence speed of the backfitting algorithm, requiring fewer iterations to reach a critical solution.

Changes:

- Fixed bugs related to the ‘omegaGrid’ and ‘omegaMin’ arguments.

- Internal functions for better structure and maintainability:
    - `mobius()` function needed for the backfitting algorithm.
    - Introduced `calculateDesignMatrix()` and `calculateSingleCosPhi()` for cleaner and more efficient code.

# FMM 0.4.0

Enhancements:

- Improved computational efficiency: The grid search has been optimized to avoid redundant calculations, and a profile likelihood approach has been incorporated into the post-optimization procedure. These refinements lead to a more efficient convergence, reducing the need for explicit parallelization (parallelize argument) and iterative grid refinement (numReps argument).

- Refinements in `plotFMM()`: Bug fixes and improvements in prediction computation, now performed on a denser, evenly spaced grid in addition to the original time points.

- Code simplifications: The arguments `numReps` and `parallelize` have been deprecated, and the default value of `showProgress` has been changed from `TRUE` to `FALSE`.

New features:

- Additional arguments in `fitFMM()`: New arguments `omegaMin`, `omegaMax`, and `omegaGrid` have been added to provide more flexibility.

- New warnings messages: Additional warnings notify users about excessively large amplitude estimates and significant gaps in time points.

Changes:

-  Repository update: The GitHub repository link in DESCRIPTION has been updated to FMMGroupVa/FMM.

-  Vignette updates: The vignette has been revised to align with the latest modifications.

-  Internal functions restructuring: 
    - Functions `bestStep1()` and `getApply()` have been relocated from FMM_internal to FMM_internal_restr.
    - The previous version of `step1FMM()` has been renamed step1FMMRestr() and moved to FMM_internal_restr.
    - Two new internal functions, `checkArguments()` and `checkSolution()`, have been introduced to improve input validation and ensure the correctness of computed results.
 

# FMM 0.3.1

Minor improvements and fixes:

-   A vignette has been added.
-   The default value of the argument showTime of fitFMM() function has been changed to FALSE.
-   `R2()` and `alwaysFalse()` functions are listed in `NAMESPACE` file.
-   A link to the GitHub repository has been added to `DESCRIPTION` file.
-   An example for plotting the original time scale on the x-axis of the ggplots has been added to the documentation of the `plotFMM()` function.
-   Add automated unit tests for `fitFMM()` function.
-   `fitFMM()` prints informative error messages when the number of observation is insufficient or data are constant.

# FMM 0.3.0

Enhancements:

-   Optimization of the estimation procedure for the monocomponent, multicomponent and restricted FMM models. An embedded parallelized procedure is available for all the models for faster estimations.
-   Code has been refined: it is now more object-oriented, code duplicity has been reduced, unnecessary orders have been suppressed, API usability has been improved, and the documentation is now automatically created with `roxygen2`.

New features:

-   In addition to the estimation by blocks previously available for the restricted model, an exact implementation has been added for more accurate fits in expense of more computational time.

# FMM 0.2.0

New features:

-   Optimized monocomponent and multicomponent model computation.

# FMM 0.1.2

Enhancements:

-   Added a slot in the FMM object with information on the number of iteraction of the fitting process.
-   Renamed the `FMMPeaks` function to `getFMMPeaks`.
-   Renamed the `generate_FMM` function to `generateFMM`. Also parameters `a`, `b` and `w` to `alpha`, `beta` and `omega`, respectively.

New features:

-   New `show()` method to print basic information on an FMM object.

# FMM 0.1.1

Fixed the issues pointed out in the review of the CRAN submission:

-   Modified the beginning of the description.
-   Added details about internal functions in the documentation.
-   Replaced the use of `cat()` by `warning()`.
-   `options(warn=-1)` is no longer used.
-   No example in the documentation of the functions uses more than two processor cores.

# FMM 0.1.0

-   Submitted version 0.1.0 to CRAN.
