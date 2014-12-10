#pvt-tool:  re-organized code for eos fitting PVT data


* Use best-fit value of Debye Temp for Tange's data (1000K)
* Robust/Normal fit on Tange
* Robust/Normal fit on Wolf
* Send the fitted results PVTdatasets from robust fits

# NEED to fix bug in  fitErrModPVTeval.m:
* dimension mismatch error does not allow single measgrp data

# TODO
* Increase by 1 decimal place the reported values
* P vs V
* Fitting changing prior on gamma0
* Fit in Std mode the Wolf2014ErrCorrDataR.md file
    * least squares
    * robust

* PVT figures
*  reduced isotherm figures
*  EOS fit files for both
* ideal mixing model/updated LLSVP figure
* Refit En87 maker mode, all one data group

