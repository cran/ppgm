# ppgm 1.1.0

* Parallel processing added to ppgm in order to reduce processing time when analyzing large tree datasets
* Added fix to getBioclimVars to assign fossils correct paleoclimate
* Added vignette for new users
* Added labelling to plotAnimatedPPGM and plotAnimatedPPGMMultiPhylo for each time slice
* Changed gif saving in AnimatedPPGM functions from animation::saveGIF (which requires imageMagick) to gifski::save_gif (works with base R)
* Corrected details of getLineageClimate function
* Added evolutionary models to ppgm and ppgmConsensus (mtrend - mean_trend from fitContinuous; rtrend - rate_trend from fitContinuous)
* Added fix in addFossil if fossil age is exactly between two paleoclimate ages - will select earlier age
* Made sure fossils plot in MESS maps
* Updated output of richnesscount in ppgm function to give all time periods

# ppgm 1.0.3

* Fixed .Rd files for CRAN resubmission

# ppgm 1.0.2

* Fixed plot.TraitGram not plotting
* Fixed MESS maps not including fossils
* Modified so that input fossil data can include min and max age
* Added layerAge argument to functions to add user to input paleoclimate where layer does not equal age
* Adds fossils to tree based on min and max age within ppgm and ppgmConsensus functions
* Fixed examples
* Added testing suite

# ppgm 1.0.1

* Fixed issues for CRAN resubmission

# ppgm 1.0.0

* Initial CRAN submission.

