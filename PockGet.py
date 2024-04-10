# PockGet.py
# single hashes for info. double hashes for debugging/alternate code
# code to run fpocket, sitehound and pocketdepth for files in a given a list 

# IMPORTANT: this code works IF sitehound, fpocket & pocketdepth folders
# as well as the directory with input proteins and this code are in the same root folder

import os, shutil, sys
root= os.getcwd() 

# path for executables
fpocket_path=os.path.join(root,'/fpocket')
sitehound_path=os.path.join(root,'/EasyMIFs_SiteHound_Linux')
pocketdepth_path=os.path.join(root,'/PocketDepth')

#---------------------------------------------------------------------------------------------------
# POCKET GENERATION PIPELINES AND FILE HANDLING

def fpocket(file_path, op_folder):

    file_name=os.path.split(file_path)[-1] # just file name

    if os.path.exists(op_folder) and os.listdir(op_folder):
        print(f"FPOCKET RESULT ALREADY EXISTS FOR {file_name}")
        return
    else:
        if os.path.exists(op_folder):
            os.rmdir(op_folder)
    
    print("RUNNING FPOCKET...")
    shutil.copy(file_path, fpocket_path) # moves to fpocket directory

    os.chdir(fpocket_path)
    os.system("fpocket -f "+ file_path)  # runs fpocket
    
    os.rename(fpocket_path+'/'+file_name[:-4]+'_out',op_folder)
    os.remove(fpocket_path+'/'+file_name)

    print()
    print("FPOCKET COMPLETE FOR",file_name)
    print("-------------------------------------------------------------")
    print()

    return

def pocketdepth(file_path, op_folder):

    file_name=os.path.split(file_path)[-1] # just file name

    if os.path.exists(op_folder) and os.listdir(op_folder):
        print(f"POCKETDEPTH RESULT ALREADY EXISTS FOR {file_name}")
        return
    else:
        if os.path.exists(op_folder):
            os.rmdir(op_folder)

    shutil.copy(file_path,pocketdepth_path) # move file to pocketdepth directory
    os.chdir(pocketdepth_path)              # change env to pocketdepth directory
    os.system("python PD.py "+ file_name)  # runs pocketdepth
 
    os.rename(pocketdepth_path+'/'+file_name[:-4],op_folder)
    os.remove(pocketdepth_path+'/'+file_name)

    print()
    print("POCKETDEPTH COMPLETE FOR",file_name)
    print("-------------------------------------------------------------")
    print()

    return
   
    
def sitehound(file_path, op_folder):
    
    file_name=os.path.split(file_path)[-1] # just file name

    if os.path.exists(op_folder) and os.listdir(op_folder):
        print(f"SITEHOUND RESULT ALREADY EXISTS FOR {file_name}")
        return
    
    else:
        if os.path.exists(op_folder):
            os.rmdir(op_folder)
            print(f"UPDATING SITEHOUND...")


    shutil.copy(file_path,sitehound_path)           # move file to sitehound directory
    os.chdir(sitehound_path)                        # change env to sitehound directory
    os.system("python2.7 sitehound.py "+ file_name) # runs sitehound

    os.rename(sitehound_path+'/'+file_name[:-4],op_folder)
    os.remove(sitehound_path+'/'+file_name)

    print()
    print("SITEHOUND COMPLETE FOR", file_name)
    print("------------------------------------------------------")  
    print() 

    return

#---------------------------------------------------------------------------------------------------
# AUTOMATION 

def main(pdb_paths, accession_list=None):

    # if no file_list given, run for all
    if accession_list==None:
        accession_list=[]
        for file in os.listdir(pdb_paths):
            file_path=os.path.join(pdb_paths,file)
            if os.path.isfile(file_path) and file.endswith('pdb'):
                if file not in file_list:
                    accession_list.append(file.strip())
        
    output_path= os.path.join(root,'/protein_results') 
    if not os.path.exists(output_path):
        os.mkdir(output_path)

        for accession in accession_list:
            ##print(accession)
            for file in os.listdir(pdb_paths):
            
                if file.endswith('.pdb') and accession in file:
                    file_address=os.path.join(pdb_paths,file)
                    ##print(file_address)
                    ##print()
                    protein_folder=os.path.join(output_path,file[:-4])

                    if not os.path.exists(protein_folder):
                        os.mkdir(protein_folder)

                    # fpocket
                    new_folder= protein_folder + '/' + file[:-4]+'_fpk'
                    fpocket(file_address, new_folder) 

                    # sitehound
                    new_folder= protein_folder + '/' + file[:-4]+'_sth'
                    sitehound(file_address, new_folder) 

                    # pocketdepth
                    new_folder= protein_folder + '/' + file[:-4]+'_pkd'
                    pocketdepth(file_address, new_folder)  
            

#---------------------------------------------------------------------------------------------------
# MAIN PROGRAM STARTS HERE

if __name__ == "__main__":

    if len(sys.argv)==4 and sys.argv[1] =='--list':
        in_dir= sys.argv[2] # path to folder containing input pdb files 
        file_list= sys.argv[3] # path to text file with list of desired accession numbers

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
        
        main(in_dir,accession_list)

        
    elif len(sys.argv)==3 and sys.argv[1]=='--all':
        in_dir= sys.argv[2] # path to folder containing input pdb files 

        if not os.path.isdir(in_dir):
            print("Invalid PDB dir path")
            exit(1)

        main(in_dir)  

    else:
        print("PockGet.py → USAGE:")
        print("For a list of files:")
        print("\tpython PockGet.py --list <path to input dir> <protein filename>")
        print("\tpython PockGet.py --list <path to input directory> <path to text file with all filenames>")
        print("For all files in the folder:")
        print("\tpython PockGet.py --all <path to input dir>")
