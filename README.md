# CICPT_1831_Rwanda

1. Prior to any analysis, test is done to check if samples are put on desired plates and chips. For this test, samples sheet prepared to run chips will be compared with the sample sheet of the output data (after running chips). 
Run `samples_on_plates_check.R` file. All the samples need to be at desired locations as per input sample sheet.
2. Discordant sex check need to be done to check if the gender on input sample sheet matched with the data generated. Run `Disccordant_Sex_check.R`. There should be not discordant sex and if exists then problematic samples need to be removed from the analysis. 
3. Cross-hybridizing probes file which is used for quality control to remove the probes can be prepared by running `MkCrossHybridData.R`. The cross hybrid probe files are download from internet for epic chips. 
4. Run `PreprocessNoob.R` for preprocessing using `Noob` method.
5. Quality control test can be performed by running `Quality_control_RawandaEWAS.R`. This will try to remove the missing values, perform normalization and plot density plots.
6. For PCA to find the variation in the data run `Combat_PCA.R`. This is done to remove batch effects.
7. To run analysis on pre-processed and quality controlled data, run `mcSEA_on_Top_variable_Probes_clean.R` file https://doi.org/10.2217/epi-2021-0310. 
8. Use `Ridge_regression.R` to perform Ridge regression analysis. 
