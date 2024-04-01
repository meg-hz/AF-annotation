# PGet.py
# single hashes for info. double hashes for debugging/alternate code
# code to run fpocket, sitehound and pocketdepth for files in a given a list of acession IDs

# IMPORTANT: this code works IF sitehound, fpocket & pocketdepth folders
# as well as the directory with input proteins and this code are in the same root folder

import os, shutil, sys
root= os.getcwd() 

# path for executables
fpocket_path=os.path.join(root,'/fpocket')
sitehound_path=os.path.join(root,'/EasyMIFs_SiteHound_Linux')
pocketdepth_path=os.path.join(root,'/PocketDepth')

def fpocket(file_path, op_folder):

    file_name=file_path.split('/')[-1] # just file name

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
    
    os.rename(fpocket_path+'/'+file[:-4]+'_out',op_folder)
    os.remove(fpocket_path+'/'+file_name)

    print()
    print("FPOCKET COMPLETE FOR",file_name)
    print("-------------------------------------------------------------")
    print()

    return

def pocketdepth(file_path, op_folder):

    file_name=file_path.split('/')[-1] # just file name

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
    
    file_name=file_path.split('/')[-1] # just file name

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
# MAIN PROGRAM STARTS HERE

if __name__ == "__main__":
    if len(sys.argv) == 3:

        pdb_paths= sys.argv[1] # path to folder containing input pdb files 
        filename_list= sys.argv[2] # path to text file with list of desired accession numbers

        # path to where the output folders go
        output_path= os.path.join(root,'/result_pockets') 
        count=0
        while os.path.exists(output_path):
            count+=1
            output_path= os.path.join(root,'/result_pockets'+str(count))
        os.mkdir(output_path)

        accession_list=[]

        if filename_list.endswith('.txt'):
            with open(filename_list,'r') as ref:    # this works if filename_list is a .txt with a list of filenames
                
                for line in ref:
                    if line.split('\n')[0] not in accession_list:
                        accession_list.append(line.split('\n')[0]) # array of accession numbers
        else:
            accession_list.append(filename_list.strip())
        

        for accession in accession_list:
            ##print(accession)
            for file in os.listdir(pdb_paths):
            
                if file.endswith('.pdb') and accession in file:
                    file_address=os.path.join(pdb_paths,file)
                    ##print(file_address)
                    ##print()
                    protein_folder=output_path+'/'+file[:-4]
            
                    # fpocket
                    new_folder= protein_folder + '/' + file[:-4]+'_fpk'
                    fpocket(file_address, new_folder) 

                    # sitehound
                    new_folder= protein_folder + '/' + file[:-4]+'_sth'
                    sitehound(file_address, new_folder) 

                    # pocketdepth
                    new_folder= protein_folder + '/' + file[:-4]+'_pkd'
                    pocketdepth(file_address, new_folder)  
    else:
        print("Usage: python PGet.py <path to input dir> <protein accession/filename>")
        print("or: python PGet.py <path to input directory> <path to filename_list.txt>")
