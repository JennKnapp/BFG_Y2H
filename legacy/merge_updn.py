
# merge all the files (upup + dndn) together in the same input dir

import os 
import glob

def merge(d):
    
    list_files = os.listdir(d)
    for f in list_files:
        if "upup" in f and "R1" in f:
            print(f)
            f_r2 = f.replace("R1", "R2")
            f_id = f.split("_")[:3]
            f_id = "_".join(f_id).replace("upup", "dndn")
            # cat 2 files together
            full_dn_R1 = glob.glob(os.path.join(d, f_id+"*R1*"))[0]
            print(glob.glob(os.path.join(d, f_id+"*R1*")))
            full_dn_R2 = glob.glob(os.path.join(d, f_id+"*R2*"))[0]
            print(glob.glob(os.path.join(d, f_id+"*R2*")))
            output_r1 = os.path.join(d, f.replace("-upup", ""))
            output_r2 = os.path.join(d, f_r2.replace("-upup", ""))
            cmd = "cat "+os.path.join(d, f) +" " + full_dn_R1 + " > " + output_r1

            os.system(cmd)
            cmd = "cat "+os.path.join(d, f_r2) +" " + full_dn_R2 + " > " + output_r2
            os.system(cmd)
        else:
            continue

if __name__ == "__main__":

    d = "/home/rothlab/rli/01_ngsdata/200729_updn/"
    merge(d)
        
