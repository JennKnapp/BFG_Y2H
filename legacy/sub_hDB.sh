
#!/bin/bash
#$ -cwd
#$ -V

/home/rothlab/rli/py/bin/python2.7 ./src/sam_rc.py --r1sam $1 --r2sam $2 --prefix $3

