# FAuto.py
# move files to flapp folder + generate pairs.txt + run FLAPP
# single hashes for info. double hashes for debugging/alternate code

# IMPORTANT: this code only works if it is in the same directory as the FLAPP directory

import os, shutil, sys

root= os.getcwd()
FLAPP_path=os.path.join(root,'FLAPP') 

# obtaining list of addresses of files 
def address_list(folderpath):
    arr_of_add=[]
    for file in os.listdir(folderpath):
        file_path=os.path.join(folderpath,file)
        ##print(file_path)
        if os.path.isfile(file_path) and file_path.endswith('.pdb'):
            arr_of_add.append(file_path)

    return arr_of_add

# writing pairwise matrix of files
def write_pairs(basename, list1, list2):

    filepath = os.path.join(FLAPP_path,f'{basename}.txt')
    count=0
    while os.path.exists(filepath):
        count+=1
        filepath = os.path.join(FLAPP_path,f'{basename}{str(count)}.txt')

    # write pair file
    with open(filepath,'w') as pairwise:
        for file1 in list1:
            file1_name=os.path.split(file1)[-1]

            for file2 in list2:
                file2_name=os.path.split(file2)[-1]

                pairwise.write(file1_name + '\t' + file2_name +'\n')   

    print(f"{os.path.split(filepath)[-1]} generated")
    return filepath

def main(mode, pocket_path, template_path, cores):

    # foldernames based on comparison mode
    if mode== '--f2self':
        basename='SELF_'
    elif mode=='--f2temp':
        basename='TEMP_'
    
    # for BindingSite
    count=0
    BS_dir=os.path.join(FLAPP_path,f'{basename}BindingSites')
    while os.path.exists(BS_dir):
        count+=1
        BS_dir=os.path.join(FLAPP_path,f'{basename}BindingSites{str(count)}')
    os.mkdir(BS_dir)

    BS_dir_asvar=os.path.split(BS_dir)[-1]
    print(f"{BS_dir_asvar} Folder Generated")

    # for SiteVector
    count=0
    SV_dir=os.path.join(FLAPP_path,f'{basename}SiteVector')
    while os.path.exists(SV_dir):
        count+=1
        SV_dir=os.path.join(FLAPP_path,f'{basename}SiteVector{str(count)}')
    SV_dir_asvar=os.path.split(SV_dir)[-1]

    # get address of all necessary files
    af_temp=address_list(pocket_path)
    template_temp=address_list(template_path)

    # writing the pairfile
    pair_file=os.path.split(write_pairs(basename, af_temp,template_temp))[-1] 
    
    # moving necessary files
    for i in af_temp:
        shutil.copy(i,BS_dir)
    print('Pocket files moved')

    for i in template_temp:
        shutil.copy(i,BS_dir)
    print('Template files moved')
     
    count=0
    outfile_path = os.path.join(FLAPP_path,f'{basename}Outfile.txt')
    while os.path.exists(outfile_path):
        count+=1
        outfile_path=os.path.join(FLAPP_path,f'{basename}Outfile'+str(count)+'.txt')
    outfile_path=os.path.split(outfile_path)[-1]
    
    os.chdir(FLAPP_path)
    print("Running FLAPP...")
    print('Generating Binary Files...')
    os.system(f'conda run -n FLAPP')
    os.system(f'python CreateVectors.py {BS_dir_asvar} {SV_dir_asvar}')
    print('Generating Alignment File...')
    os.system(f'python FLAPP.py {SV_dir_asvar} {pair_file} {outfile_path} {cores}')

    return

if __name__ == "__main__":
    if len(sys.argv) ==5 and sys.argv[1] =='--f2temp':
        mode=sys.argv[1]
        af_dir = sys.argv[2]    # path to pocket PDB files
        temp_dir = sys.argv[3]  # path to template directory
        core_no= sys.argv[4]    # no. of cores
        
        if not (os.path.isdir(af_dir) and os.path.isdir(temp_dir)):
            print("Invalid path name(s). Exiting Code.")
            exit(1)
    
        print("Running FLAPP...")
        main(mode,af_dir,temp_dir,core_no)
    
    elif len(sys.argv)==4 and sys.argv[1] == '--f2self':
        mode=sys.argv[1]
        af_dir = sys.argv[2]
        core_no= sys.argv[3]
    
        if not os.path.isdir(af_dir):
            print("Invalid path name(s). Exiting Code.")
            exit(1)

        print("Running FLAPP...")
        main(mode,af_dir,af_dir,core_no)

    else:
        print("FAuto.py → USAGE:")
        print("To run FLAPP:")
        print("\tpython FAuto.py --f2temp <path to pocket-file-dir> <path to template directory> <no of cores>\t OR")
        print("\tpython FAuto.py --f2self <path to pocket-file-dir> <no of cores>\n")
