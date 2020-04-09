# humanAML_MTOR

This is the survival analysis with AML patients derived microarray data and mTORC1 related gene signatures.

(Data sets)

AML transcriptome data was originally derived from Blood 2008 Nov 15;112(10):4193-201 (PMID: 18716133) and the data was deposited at GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12417)(GSE12417). Expression data and meta data was created based on the site as above. 

mTORC1 realed gene signatures were used 

1.RAPAMYCIN_SENSITIVE_GENES
2.HALLMARK_MTORC1_SIGNALING
3.CUNNINGHAM_RAPAMYCIN_DN
4.mouse DEG 

(Process)

1. AML patients sample were stratified with hierachial clustering based on the expression of mTORC1 relaed gene signatures. 
2. mTORC1 high and low culsters are identified and survival curve is drawn for these two clusters.
3. p-value will be caluculated at day 300, 600, 900 and overall
4. Differentially expressed gene (q-value < 0.01) will be identified. 

(Software)
The program was written by R 3.6.1

