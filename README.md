# microarray_transformation_evaluation
# the input for the variance analyses (wether it is avg based or fraction based) is a csv that was extracted from the Tableau of a model version, for a specific comparison.
# You can look for examples in the "/source_csvs" folder.
# The mandatory columns are:
#1 Cell Typeold - sometimes it appreas in the Tableau as "Cell Type", so you need to manualy change it to "Typeold"
#2 Platform - you need to set it for one of 2 optional values per row - "microarray", "rnaseq".
#3 Experiment ID
#4 Estimate Up Down
#5 significant - You need to set this column manually and mark as TRUE/FALSE each row, based on the FDR threshold that you want. If you want to include all the cttest results without any limitations - set all the rows to "TRUE". 