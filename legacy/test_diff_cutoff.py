import argparse
import os

def assign_cutoff(cut_off, r1_csv, r2_csv, redirect_dir, sh_dir, name):
    """
    cut_off: list of int specify the cut-offs
    r1_csv:
    r2_csv:
    redirect_dir: directory to save the read count files
    using different curoffs for alignemnt files 
    """
    for c in cut_off:
        c = str(c)
        sh_sub = "#!/bin/bash \n" + \
              "#$ -N " + name + "-" +  str(c) +"\n" + \
              "#$ -cwd \n"
        # make bash file for submitting
        python_cmd = "/home/rothlab/rli/py/bin/python2.7 /home/rothlab/rli/02_dev/08_bfg_y2h/bfg_analysis/main.py " +\
                "--r1 "+ r1_csv + \
                " --r2 " + r2_csv + \
                " -r " + redirect_dir + \
                " -cut " + c + \
                " -prefix " + name + \
                " --mode yeast "
        sh_fname = name+"_sub_"+str(c)+".sh"
        sh_f = os.path.join(sh_dir, sh_fname)
        with open(sh_f, "w") as sh:
            sh.write(sh_sub)
            sh.write(python_cmd)
        cmd = "chmod 755 "+sh_f
        os.system(cmd)
        #cmd = "qsub " + sh_f
        #print(cmd)
       # os.system(cmd)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="")
    # input of the script:
    parser.add_argument("-dir", help="input dir contains r1 csv and r2 csv", required=True)
    parser.add_argument("-out", help="output dir to save all the read counts", required=True)

    parser.add_argument("-prefix", help="prefix to add to the names", required=True)

    args = parser.parse_args()
    sh_dir = "/home/rothlab/rli/02_dev/08_bfg_y2h/tmp_sh/"

    # find r1_csv and r2_csv in dir
    for f in os.listdir(args.dir):
        if "R1" in f and "_noh.csv" in f:
            if "high" in f:
                r1_high = os.path.join(args.dir,f)
            elif "med" in f:
                r1_med = os.path.join(args.dir,f)
            else:
                r1_pre = os.path.join(args.dir,f)

        if "R2" in f and "_noh.csv" in f:
            if "high" in f:
                r2_high = os.path.join(args.dir,f)
            elif "med" in f:
                r2_med = os.path.join(args.dir,f)
            else:
                r2_pre = os.path.join(args.dir,f)
            
    # test
    cut_off=range(3,21)
    print(cut_off)
    name=args.prefix

    # GFP high
    assign_cutoff(cut_off, r1_high, r2_high, args.out, sh_dir, name+"-high")
    # GFP pre
    assign_cutoff(cut_off, r1_pre, r2_pre, args.out, sh_dir, name+"-pre")
    # GFP med
    assign_cutoff(cut_off, r1_med, r2_med, args.out, sh_dir, name+"-med")

