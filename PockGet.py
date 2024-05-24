# PockGet.py
# code to run fpocket, sitehound and pocketdepth for files in a given a list 
# single hashes for info. double hashes for debugging/alternate code

# IMPORTANT: this code works IF sitehound, fpocket & pocketdepth folders
# as well as the directory with input proteins and this code are in the same root folder

import os, shutil, sys
root= os.getcwd() 

# path for executables
fpocket_path=os.path.join(root,'fpocket')
sitehound_path=os.path.join(root,'EasyMIFs_SiteHound_Linux')
pocketdepth_path=os.path.join(root,'PocketDepth')

#---------------------------------------------------------------------------------------------------
# POCKET GENERATION PIPELINES AND FILE HANDLING

def fpocket(file_path, op_folder):

    #print(op_folder)

    file_name=os.path.split(file_path)[-1] # just file name

    if os.path.isdir(op_folder):
        if os.listdir(op_folder):
            print(f"FPOCKET RESULT ALREADY EXISTS FOR {file_name}")
            return
        else:
            os.rmdir(op_folder)
            print(f"UPDATING FPOCKET...")
            
    
    print(F"RUNNING FPOCKET FOR {file_name}...")
    shutil.copy(file_path, fpocket_path) # moves to fpocket directory
    os.chdir(fpocket_path)
    os.system("fpocket -f "+ os.path.join(fpocket_path,os.path.split(file_path)[-1]))  # runs fpocket

    os.rename(os.path.join(fpocket_path,file_name[:-4]+'_out'),op_folder)
    os.remove(os.path.join(fpocket_path,file_name))

    print()
    print("FPOCKET COMPLETE")
    print("-------------------------------------------------------------")
    print()

    return

def pocketdepth(file_path, op_folder):

    file_name=os.path.split(file_path)[-1] # just file name

    if os.path.exists(op_folder):
        if os.listdir(op_folder):
            print(f"POCKETDEPTH RESULT ALREADY EXISTS FOR {file_name}")
            return
        else:
            os.rmdir(op_folder)
            print(f"UPDATING POCKETDEPTH...")

    print(F"RUNNING POCKETDEPTH FOR {file_name}...")
    shutil.copy(file_path,pocketdepth_path) # move file to pocketdepth directory
    os.chdir(pocketdepth_path)              # change env to pocketdepth directory
    os.system("python PD.py "+ file_name)  # runs pocketdepth
 
    os.rename(os.path.join(pocketdepth_path,file_name[:-4]),op_folder)
    os.remove(os.path.join(pocketdepth_path,file_name))

    print()
    print("POCKETDEPTH COMPLETE FOR",file_name)
    print("-------------------------------------------------------------")
    print()

    return
   
    
def sitehound(file_path, op_folder):
    
    file_name=os.path.split(file_path)[-1] # just file name

    if os.path.exists(op_folder):
        if os.listdir(op_folder):
            print(f"SITEHOUND RESULT ALREADY EXISTS FOR {file_name}")
            return
    
        else:
            os.rmdir(op_folder)
            print(f"UPDATING SITEHOUND...")

    print(F"RUNNING SITEHOUND FOR {file_name}...")
    shutil.copy(file_path,sitehound_path)           # move file to sitehound directory
    os.chdir(sitehound_path)                        # change env to sitehound directory
    os.system("python2.7 sitehound.py "+ file_name) # runs sitehound

    os.rename(os.path.join(sitehound_path,file_name[:-4]),op_folder)
    os.remove(os.path.join(sitehound_path,file_name))

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
            if os.path.isfile(file_path) and file.endswith('.pdb'):
                if file not in accession_list:
                    accession_list.append(file.strip())
        
    output_path= os.path.join(root,'PocketResults') 
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
                new_folder= os.path.join(protein_folder,file[:-4]+'_fpk')
                fpocket(file_address, new_folder)                 

                # sitehound
                new_folder= os.path.join(protein_folder,file[:-4]+'_sth')
                sitehound(file_address, new_folder) 

                # pocketdepth
                new_folder= os.path.join(protein_folder,file[:-4]+'_pkd')
                pocketdepth(file_address, new_folder)  
            

def main2(pdb_paths, mode, accession_list=None):
    
    if accession_list==None: # basically accounts for the all protein files in the path if a list is not given
        accession_list = list(set([file for file in os.listdir(pdb_paths) if os.path.isfile(os.path.join(pdb_paths,file)) and file.endswith('.pdb')]))
    
    output_path= os.path.join(root,f"{mode[2:]}_results")
    count=0
    while os.path.exists(output_path):
        count+=1
        output_path= os.path.join(root,f"{mode[2:]}_results({count})")
    os.mkdir(output_path)

    def themode(pockmode,filepath,op_path):
        if pockmode == '--fpocket':
            return fpocket(filepath,op_path)
        elif pockmode == '--sitehound':
            return sitehound(filepath,op_path)
        elif pockmode == '--pocketdepth':
            return sitehound(filepath,op_path)


    for accession in accession_list:
        for file in os.listdir(pdb_paths): 
            if file.endswith('.pdb') and accession in file:
                file_address=os.path.join(pdb_paths,file)
                ##print(file_address)
                ##print()

                protein_folder=os.path.join(output_path,f'{file[:-4]}_fpk') 
                themode(mode,file_address,protein_folder)

#---------------------------------------------------------------------------------------------------
# TERMINAL MODULE STARTS HERE

if __name__ == "__main__":

    if (len(sys.argv)==4 and sys.argv[2]=='--dir') or (len(sys.argv)==5 and sys.argv[2]=='--list'):

        # accounting for folder containing PDB files
        input_dir=sys.argv[3]
        if not os.path.isdir(input_dir):
            print("Invalid PDB dir path. Exiting Code...")
            exit(1)

        # accounting for list of protein accession numbers
        if len(sys.argv)==4:
            accession_list=None
        else:
            file_list=sys.argv[4]

            if not (os.path.isfile(file_list) and file_list.endswith('.txt')):
                print("Invalid Reference File Path. Exiting Code...")
                exit(1)

            with open(file_list, 'r') as ref:
                lines=ref.readlines()
            accession_list=sorted(list(set([i.split()[0] for i in lines])))

        # accounting for mode
        runmode= sys.argv[1]

        if runmode== '--all':
            main(input_dir,accession_list)
        elif runmode in ('--fpocket','--sitehound','--pocketdepth'):
            main2(input_dir,runmode,accession_list)       
        

    else:
        print("PockGet.py → USAGE:")
        print("To run all three pocket prediction tools ")
        print("\tFor a list of files:")
        print("\t\tpython PockGet.py --all --list <path to input directory> <path to reference file>")
        print("\tFor all files in the folder:")
        print("\t\tpython PockGet.py --all --dir <path to input dir>")
        print()
        print("To run a specific pocket prediction tool ")
        print("(replace <mode> with '--fpocket','--sitehound' or '--pocketdepth')")
        print("\tFor a list of files:")
        print("\t\tpython PockGet.py <mode> --list <path to input directory> <path to reference file>")
        print("\tFor all files in the folder:")
        print("\t\tpython PockGet.py <mode> --dir <path to input dir>")
        print()
