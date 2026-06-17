**TFBScount**

R code to count the number of transcription factor binding sites (TFBS) in JASPAR output

**Introduction**

JASPAR is an online website that takes in a particular regulatory element (e.g. promoter or enhancer) and uses binding site matrices for different transcription factors to calculate all the possible TFBS and their relative profile score. Users provide a regulatory sequence and select the transcription factors to include and set a relative profile score threshold. We assume that comparable lengths of regulatory sequences are analyzed for a set of species with the results saved in separate files.
This program reads in those files and then counts the number of TFBS for each sequence above a series of relative profile score thresholds. This version is set up to compare TFBS counts for fourteen transcription factors (TFs) from JASPAR analysis of short wavelength sensitive 2a (SWS2a) opsin promoters. It compares counts for five cichlid species and to five other fishes. It uses Welch statistics to determine if the cichlid values are less than or greater than those of the other fish.
Prior to running this program, the regulatory sequences need to be analyzed by JASPAR. For this example, we added 36 binding matrices for a set of 14 TFs to the cart : Crem, Foxb1, Nfix, Rara::Rxrg, Rarg, Rax2, Rorc, Rxrg, Six3, Sox5, Sox6, Tbx2, Tbx4, and Thrb. These were then analyzed and the tab delimited results saved to a file, with one file for each of ten species. No header is needed.
In this program, the JASPAR results for 10 species (10 files) are each read into a table. The columns of the JASPAR output are:
"MatrixID", "Name", "Score", "Relative_score", "SequenceID", "Start", "End", "Strand", and "PredictedSequence". The function File_TFnames is used to convert the Name column, which contains the MatrixID.TFname to just TFname and makes it all capitals.
The JASPAR results are used to count the number of TFBS for each TF for increasing relative profile score thresholds. For example, JASPAR usually considers relative scores that are 0.8 or above. We count the number of TFs with TFBS above 0.8, 0.82, 0.84 … up to 1.0. The counts for different species can be plotted versus relative score threshold, Welch statistics can be calculated, and the location of TFBS across species can be plotted.
