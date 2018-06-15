
python /home/dana/workspace/gtex-estrs/Scripts/Linear_Regression/RegressionAssociationTest.py \
	--expr ATXN2_PSI.csv  \
	--exprannot /storage/dana/spliceSTR/expression/gencode.v19.annotations_exons.csv \
	--strgt Allele_Gentype_ATXN2.table \
	--distfromgene 50 \
	--out ATXN1_regression.csv \
	--chrom chr16 \
	--linear --norm --alleles --tmpdir tmp --ingene
