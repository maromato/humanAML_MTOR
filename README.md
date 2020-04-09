# humanAML_MTOR

This is the survival analysis with AML patients derived microarray data and mTORC1 related gene signatures.

# 1. Data sets

AML transcriptome data was originally derived from Blood 2008 Nov 15;112(10):4193-201 (PMID: 18716133) and the data was deposited at GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12417)(GSE12417). Expression data (GSE12417_mod_4.txt) and meta data (GSE12417-GPL96_meta.txt) was created based on the site as above. 

# 2. Gene signature sets

mTORC1 realed gene signatures were used 

1.RAPAMYCIN_SENSITIVE_GENES (MSigDB)

2.HALLMARK_MTORC1_SIGNALING (MSigDB) 

3.CUNNINGHAM_RAPAMYCIN_DN  (MSigDB)

4.mouse DEG (OUR DATA SETS) -- Not currently included. 

# 3. Overview of the Process

1. AML patients sample will be stratified with hierachial clustering based on the expression of mTORC1 relaed gene signatures. 
2. mTORC1 high and low culsters will beidentified and survival curve is drawn for these two clusters.
3. p-value will be caluculated at day 300, 600, 900 and overall.
4. Differentially expressed gene (q-value < 0.05) will be identified.

# 4. Programs files (gene signature/R code)

1.RAPAMYCIN_SENSITIVE_GENES (MSigDB) :AML_MTORP_RAP.R

2.HALLMARK_MTORC1_SIGNALING (MSigDB) :AML_MTORP_Hallmarks.R

3.CUNNINGHAM_RAPAMYCIN_DN  (MSigDB) :AML_MTORP_CUNNING.R

4.mouse DEG (OUR DATA SETS) :AML_MTORP_DEG.R

# 5. Software and miscellaneous 
The programs were written by R 3.6.1.

This is a beta-test version for educational purposes and please let me know if there are any questions (toshihiko.oki@gmail.com).


