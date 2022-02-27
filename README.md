# AdaptiveMRMC
Adaptive designs in Multi-Reader Multi-Case (MRMC) clinical trials of imaging devices

## License Statement
Please read our [license statement](https://github.com/DIDSR/AdaptiveMRMC/blob/main/LICENSE-STATEMENT.md)

## Cite our work
The methodology behind the code shared here is in the following paper:

Huang Z, Samuelson F, Tcheuko L, Chen W. Adaptive designs in multi-reader multi-case clinical trials of imaging devices. Statistical Methods in Medical Research. 2020;29(6):1592-1611. doi:10.1177/0962280219869370

## Basic User Manual

1. Install iMRMC R package from CRAN or [github](https://github.com/DIDSR/iMRMC). Note that this package requires Java 1.7 or above.

2. Install [AdaptiveMRMC](https://github.com/DIDSR/AdaptiveMRMC) R pakcage by downloading the source package [released on github](https://github.com/DIDSR/AdaptiveMRMC/releases/download/1st_release/AdaptiveMRMC_1.0.0.tar.gz)
3. Run the workflow example under *your_package_path/AdaptiveMRMC/Workflow_example/demoAdaptiveMRMC.R* This workflow example shows the four steps in designing and analyzing an adaptive MRMC study using adaptive method 1 (adaptively re-sizing readers) as described in our paper. A similar workflow can be followed using adaptive method 2 (adaptively re-sizing both readers and cases).

  `library(AdaptiveMRMC)`  
  `source(paste0(.libPaths()[1], '/AdaptiveMRMC/Workflow_example/demoAdaptiveMRMC.R'))` 
  
4. The code for simulation studies in our paper is included in the package under *your_package_path/AdaptiveMRMC/Simulations*. Please note these are running codes customized to our computer cluster and may not run on your computer. The purpose of including these raw research code is for reference, e.g., by reading the code one may find how the adaptive methods are applied to the simulated data.
