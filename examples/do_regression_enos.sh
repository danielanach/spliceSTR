python /home/dana/workspace/gtex-estrs/Scripts/Linear_Regression/RegressionAssociationTest.py \
        --expr enos_PSI.csv  \
        --exprannot /storage/dana/spliceSTR/expression/gencode.v19.annotations_exons.csv \
        --distfromgene 50 \
        --strgt Allele_Gentype_enos.table \
        --out enos_regression.csv \
        --chrom chr7 \
        --norm  --quadratic --linear --alleles --tmpdir tmp --ingene
