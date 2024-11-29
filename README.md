# Export_BW_HPC
Export_BW function adapted from user Baboon61 (Bastien Hervé) who in turn adapted the function from Arch_R.

I implemented paralellization so that the function can be easily manipulated based on the number of cores called by an Rstudio session. 
The original program vesion utilizes future to split the fragment files and run them in parallel. I have decided to use doParallel instead as my experience when running the code in an R session has shown that future did not parallelize the process. Also made fixes to minor errors encountered at runtime.

## Citation
SplitFragments RCode adopted and modified from Bastien Herve (Baboon61) 
Apr 4, 2023. [https://github.com/stuart-lab/signac/issues/28]
