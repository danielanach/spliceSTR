snp_genos='/storage/szfeupe/Runs/650GTEx_estr/SNP_Analysis/chr3.tab'

python /home/dana/workspace/gtex-estrs/Scripts/Linear_Regression/oldLinRegAssociationTest.py \
	--expr FXR1_PSI.csv  \
	--exprannot /storage/dana/spliceSTR/expression/gencode.v19.annotations_exons.csv \
	--distfromgene 50 \
	--strgt $snp_genos \
	--out FXR1_regression.csv \
	--norm
	--chrom chr3 \
	--tmpdir tmp
