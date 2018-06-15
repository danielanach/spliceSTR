python /home/dana/workspace/gtex-estrs/Scripts/Linear_Regression/RegressionAssociationTest.py \
	--expr CFTR_PSI.csv  \
	--exprannot /storage/dana/spliceSTR/expression/gencode.v19.annotations_exons.csv \
	--distfromgene 50 \
	--strgt Allele_Gentype_CFTR.table \
	--out CFTR_regression_no_anova.csv \
	--chrom chr7 \
	--norm  --quadratic --linear --alleles --tmpdir tmp --ingene

