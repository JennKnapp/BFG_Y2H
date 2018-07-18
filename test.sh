# test for create fasta
#/home/rothlab/rli/py/bin/python2.7 ./src/main.py --ad ./summary/example_AD.csv --db ./summary/example_DB.csv --pfasta ./path_to_fasta/

# test for build
#/home/rothlab/rli/py/bin/python2.7 ./src/main.py --pfasta ./path_to_fasta/

# test for alignment
#/home/rothlab/rli/py/bin/python2.7 ./src/main.py --pfastq /home/rothlab/rli/01_ngsdata/180709_merged_yAD1DB4/ --output ./test_output/

# test for analysis
/home/rothlab/rli/py/bin/python2.7 ./src/main.py --ad ./summary/example_AD.csv --db ./summary/example_DB.csv --r1 yAD1DB4_GFP_high_R1_AD_BC.sam --r2 yAD1DB4_GFP_high_R2_DB_BC.sam 


