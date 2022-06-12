The output of step 4 was used as the BLSSpeller input. BLSSpeller was set up on the UGent-HPC. Both the motif table output from the first step of the BLSSpeller algorithm and the BED file output of the second step of BLSSpeller are retreived and used in subsequent analyses.

processBlsOutput.sh combines the results, selects for OG count cutoff and performs QC checks.
getUniqueTrueAndShuffledMotifs.sh divides the motifs of the shuffled OGs (biological negative control) and the true data, into the overlapping and unique motifs for both motif sets.
