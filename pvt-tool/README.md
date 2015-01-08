#Notes for pvt-tool code design and implementation

##Measurement Groups
* Current implementation is overly complex:
    * meas groups are defined as integer IDs in the PVT data file and then stringIDs are defined for each one in the input file
    * When read in, the integers are replaced with their stringIDs and stored that way
    * The error model then fits the errorbar adjustment factors by considering the *unique* stringIDs in the measGrpID list
* Problem:
    * There is a lot of overhead storing the measGrps as strings rather than just the input integer values (cellstr, string operations, etc)
    * Since the input list of measurement grp stringIDs is not stored, measurement groups are considered based on the *unique* subset
    * This causes them to be *alphabetically sorted*, rather than using the order defined in the input file!
* Solution (to be implemented):
    * Store meas groups as simple integers (same as PVT data file)
    * Also store meas group ID names
    * Upon reading the input file, verify that the measurement group IDs are non-negative integers
        * Zero values could be used as a flag to ignore (mask-out) certain data points
    * error Model parameters are considered in order according to the integer meas grp IDs (intuitive)
    * Final fitted values are easily examined and interpreted, since the measGrpID strings are stored in the PVTeval structure!

##Role of Temperature Errors
* Random temperature errors are often thought to play an important role in total propagated error when fitting P-V-T datasets
* If we fit our data in the more sensible Vmark-V-T space, error propagation changes considerably:
    * errors are entirely independent
    * errors in Vmark and Vsamp contribute significantly
    * random errors in temperature only affect the pressure misfit in proportion to the thermal pressure differences between marker and sample
    * BUT thermal pressures for materials at high pressure are VERY SIMILAR, even for neon and perovskite (very different materials!)
    * So the random errors in T contribute negligibly to the total propagated pressure error (typically less than 0.1 GPa)
* However, temperature uncertainties ARE important!:
* Temperature gradients are often serious in experiments, and these increase uncertainties dramatically:
    * T-gradients lead to P-gradients, which induces some flow at high temperature to relieve stress
    * Relaxation causes V-gradients, which can broaden/blur XRD lines
    * T value for sample and marker may be different, leading to potential systematic errors
* Temperature errors prescription:
    * Propagate random temperature errors as usual for Vmark-V-T analyses (even though contribution is small)
    * If non-neglible T-gradients are suspected (i.e. with laser-heating), separate out heated and unheated measurements into different measGroups
    * This allows separate adjustment of errorbars for heated and unheated samples, enabling the extra empirical scatter of heated data to speak for itself
