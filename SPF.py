# SPF.py
# single hashes for info. double hashes for debugging/alternate code
# code to run superfamily.pl script automatically for files in a given a list 

# IMPORTANT: this code works IF it is in the same parent directory as the superfamily executable dir 

import os, shutil, sys, PDB_modules as PDBm
root= os.getcwd() 

counter = 0

def run_spf(pdb_dir,file_list=None):

# if no file_list given, run for all
    if accession_list==None:
        accession_list=[]
        for file in os.listdir(pdb_dir):
            file_path=os.path.join(pdb_dir,file)
            if os.path.isfile(file_path) and file.endswith('pdb'):
                if file not in file_list:
                    accession_list.append(file.strip())
                
    total = len(file_list)    

    # overall output directory
    out_dir= os.path.join(root,'/protein_results') 
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # running through all files
    for file in os.listdir(pdb_dir):
        file_path= os.path.join(pdb_dir,file)

        #file validation
        if not (os.path.isfile(file_path) and file.endswith('.PDB')):
            continue
        elif file not in file_list:
            if file.split('-')[1] not in file_list:
                continue


        # generate specific protein output directory
        file_out= os.path.join(out_dir,file[:-4])
        if not os.path.exists(file_out):
            os.mkdir(file_out)    
        
        work_dir=os.path.join(root,'/superfamily')
        PDBm.get_solo_fasta(file_path,work_dir) # make .fa from .pdb
        fa_name=file[:-4] 

        os.chdir(work_dir)
        os.system("./superfamily.pl " + file) # run superfamily on the .fa

        # move result files
        shutil.move(f"{work_dir}/{fa_name}.ass", file_out)
        shutil.move(f"{work_dir}/{fa_name}.html", file_out)
        shutil.move(f"{work_dir}/scratch/{fa_name}.res", file_out)
        shutil.move(f"{work_dir}/{fa_name}", file_out)

        # remove the _torun fasta file
        os.remove(f"/scratch/{fa_name}_torun.fa")

    print(f"completed {counter}/{total} files")
    print("\n\n")

if __name__ == "__main__":
    if len(sys.argv)==4 and sys.argv[1] =='--list':
        in_dir= sys.argv[2] 
        file_list= sys.argv[3]

        if not os.path.isdir(in_dir):
            print("Invalid PDB dir path")
            exit(1)

        accession_list=[]
        if os.path.isfile(file_list) and file_list[-4:] in ['.txt', '.tsv', '.csv']:
            with open(file_list, 'r') as ref:  # this works if file_list is a file with a list of filenames
                for line in ref:
                    if line.split()[0] not in accession_list:
                        accession_list.append(line.split('\n')[0]) # array of file names
        else:
            accession_list.append(file_list.strip())
        
        run_spf(in_dir,accession_list)

    elif len(sys.argv)==3 and sys.argv[1] =='--all':
        in_dir= sys.argv[2] 

        if not os.path.isdir(in_dir):
            print("Invalid PDB dir path")
            exit(1)

        run_spf(in_dir)
    
    else:
        print("SPF.py → USAGE:")
        print("For a list of files:")
        print("\tpython SPF.py --list <path to input dir> <protein filename>\t OR")
        print("\tpython SPF.py --list <path to input directory> <path to text file with all filenames>")
        print("For all files in the folder:")
        print("\tpython SPF.py --all <path to input dir>")
