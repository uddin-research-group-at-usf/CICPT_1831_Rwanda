# CICPT_1831_Rawanda
## Steps:
1. Prior to any analysis, test is done to check if samples are put on desired plates and chips. For this test, samples sheet prepared to run chips will be compared with the sample sheet of the output data (after running chips). 
Run `samples_on_plates_check.R` file. 
2. Discordant sex check need to be done to check if the gender on input sample sheet matched with the data generated. Run `Disccordant_Sex_check.R`
3. Cross-hybridizing probes file which is used for quality control to remove the probes can be prepared by running `MkCrossHybridData.R`
4. Run `PreprocessNoob.R` for preprocessing using `Noob` method.
5. Once platting testing is done, quality control is performed by running `Quality_control_RawandaEWAS.R`.
