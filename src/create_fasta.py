import pandas as pd
import os
import param


def create_fasta(AD_summary, DB_summary, output_path, group_spec=False, AD="G0", DB="G0"):

    """
    Generate fasta file from summary
    See example_summary as templete
    """
    # set summary in param.py
    AD_summary = pd.read_table(AD_summary, sep="\t")
    DB_summary = pd.read_table(DB_summary, sep="\t")
    
    if group_spec: # make a group specific fasta file
        # select the group from AD and DB
        AD_summary = AD_summary[AD_summary.Group==AD]
        DB_summary = DB_summary[DB_summary.Group==DB]
        # fasta filename
        f_ad = "h_AD_"+AD+".fasta" 
        f_db = "h_DB_"+DB+".fasta"
    else:
        f_ad = "h_AD_G0.fasta"
        f_db = "h_DB_G0.fasta"

    with open(os.path.join(output_path,f_ad), "w") as ad:
        # grep sequence and sequence name from summary file
        # >G1;YDL169C_BC-1;7;up
        # CCCTTAGAACCGAGAGTGTGGGTTAAATGGGTGAATTCAGGGATTCACTCCGTTCGTCACTCAATAA
        for index, row in AD_summary.iterrows():
            up_seq_name = ">"+row.Group+";"+row.Locus+";"+str(row.Index)+";"+"up"
            ad.write(up_seq_name+"\n")
            ad.write(param.AD_Up1+row.UpTag_Sequence+param.AD_Up2+"\n") # add padding sequences

            dn_seq_name = ">"+row.Group+";"+row.Locus+";"+str(row.Index)+";"+"dn"
            ad.write(dn_seq_name+"\n")
            ad.write(param.AD_Dn1+row.DnTag_Sequence+param.AD_Dn2+"\n")

    with open(os.path.join(output_path, f_db), "w") as db:
            
        for index, row in DB_summary.iterrows():
            up_seq_name = ">"+row.Group+";"+row.Locus+";"+str(row.Index)+";"+"up"
            db.write(up_seq_name+"\n")
            db.write(param.DB_Up1+row.UpTag_Sequence+param.DB_Up2+"\n")

            dn_seq_name = ">"+row.Group+";"+row.Locus+";"+str(row.Index)+";"+"dn"
            db.write(dn_seq_name+"\n")
            db.write(param.DB_Dn1+row.DnTag_Sequence+param.DB_Dn2+"\n")


# Modified version to make fasta files for Miha's ORFs. 
def create_fasta_miha(AD_summary, DB_summary, output_path):
    AD_summary = pd.read_table(AD_summary, sep="\t")
    DB_summary = pd.read_table(DB_summary, sep="\t")
    mihas_AD = AD_summary[AD_summary.Plate.str.contains("Miha")]
    mihas_DB = DB_summary[DB_summary.Plate.str.contains("Miha")]
    
    f_ad = "y_AD_M.fasta"
    f_db = "y_DB_M.fasta"
    
    with open(os.path.join(output_path,f_ad), "w") as ad:
        for index, row in mihas_AD.iterrows():
            up_seq_name = ">"+row.Group+";"+row.Locus+";"+str(row.Index)+";"+"up"
            ad.write(up_seq_name+"\n")
            ad.write(param.AD_Up1+row.UpTag_Sequence+param.AD_Up2+"\n") # add padding sequences
            dn_seq_name = ">"+row.Group+";"+row.Locus+";"+str(row.Index)+";"+"dn"
            ad.write(dn_seq_name+"\n")
            ad.write(param.AD_Dn1+row.DnTag_Sequence+param.AD_Dn2+"\n")

    with open(os.path.join(output_path, f_db), "w") as db:
        for index, row in mihas_DB.iterrows():
            up_seq_name = ">"+row.Group+";"+row.Locus+";"+str(row.Index)+";"+"up"
            db.write(up_seq_name+"\n")
            db.write(param.DB_Up1+row.UpTag_Sequence+param.DB_Up2+"\n")
            dn_seq_name = ">"+row.Group+";"+row.Locus+";"+str(row.Index)+";"+"dn"
            db.write(dn_seq_name+"\n")
            db.write(param.DB_Dn1+row.DnTag_Sequence+param.DB_Dn2+"\n")


def build_index(fasta_file, output_dir):
    """
    Bowtie build for fasta_file
    """
    basename = os.path.basename(fasta_file).split(".")[0]
    #print basename
    cmd = param.BOWTIE2_BUILD+fasta_file+" "+os.path.join(output_dir,basename)
    #print cmd
    os.system(cmd)

if __name__ == "__main__":

    AD_summary = param.AD_summary
    DB_summary = param.DB_summary
    # set for our experiment 
    # AD G1-4
    # DB G1-4
    for i in range(1,5):
        create_fasta(AD_summary, DB_summary, param.REF_PATH, group_spec=True, AD="G"+str(i), DB="G"+str(i))

    # example of create fasta for all
    # create_fasta(AD_summary, DB_summary, fasta_output)
    
    # for Miha's ORFs
    #create_fasta_miha(AD_summary, DB_summary, param.REF_PATH)

    # bowtie build
    list_fasta = os.listdir(param.REF_PATH)
    for fasta in list_fasta:
        create_fasta.build_index(os.path.join(param.REF_PATH, fasta), param.REF_PATH)

