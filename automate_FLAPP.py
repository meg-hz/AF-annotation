# FLAPP.py
# single hashes for info. double hashes for debugging/alternate code
# move files to flapp folder + generate pairs.txt + run FLAPP

# IMPORTANT: this code only works if it is in the same directory as the FLAPP directory
# to change this, modify root variable to path of dir containing FLAPP folder

import os, shutil, sys

root= os.getcwd() 
FLAPP_path=os.path.join(root,'FLAPP') 

# obtaining list of addresses of pocket files 
def AF_add(ip_path):
    arr_of_add=[]
    for accession in os.listdir(ip_path):
        ##print(accession)
        site1_path= os.path.join(ip_path, accession, 'all_match')
        for file in os.listdir(site1_path):
            file_path=os.path.join(site1_path,file)
            ##print(file_path)
            if os.path.isfile(file_path) and file_path.endswith('.pdb'):
                arr_of_add.append(file_path)

    return arr_of_add

# obtaining list of addresses of CRD files 
def CRD_add(site2_path):
    arr_of_add=[]
    for file in os.listdir(site2_path):
        file_path=os.path.join(site2_path,file)
        ##print(file_path)
        if os.path.isfile(file_path) and file_path.endswith('.pdb'):
            arr_of_add.append(file_path)

    return arr_of_add

# writing pairwise matrix of files
def write_pairs(list1, list2):

    filepath = os.path.join(FLAPP_path,'Pairs.txt')

    if os.path.exists(filepath):
        response=input("Pairs.txt already exists. Do you want to rewite it (Y/N)? ").strip()
        if  response in ['N','n']:
            print("Rewriting Cancelled. Running FLAPP with existing Pairs.txt file...")
            return
        
        elif response not in ['Y', 'y']:
            print("Invalid response. Exiting Code")
            sys.exit(1)
    
    with open(filepath,'w') as pairwise:
        for file1 in list1:
            file1_name=file1.split('/')[-1].strip()
            for file2 in os.listdir(list2):
                file2_name=file2.split('/')[-1].strip()
                pairwise.write(file1_name + '\t' + file2_name +'\n')         
    print("Pairs.txt generated")
    return

# Making a New Folder with all binding sites
def FLAPP_BS():
    count=0
    BS_dir=os.path.join(FLAPP_path,'BindingSites')
    while os.path.exists(BS_dir):
        count+=1
        BS_dir=os.path.join(FLAPP_path,'BindingSites'+str(count))
    os.mkdir(BS_dir)
    return BS_dir 

def main(pocket_path, crd_path):

    # get the file addresses
    af_temp=AF_add(pocket_path)

    nrsite_path=os.path.join(crd_path,'CustomDatas','nrsiteActual')
    crd_temp=CRD_add(nrsite_path)

    write_pairs(af_temp,crd_temp)

    bindingsite_dir=FLAPP_BS()
    print(f"{bindingsite_dir.split('/')[-1]} Folder Generated")

    for i in af_temp:
        shutil.copy(i,bindingsite_dir)
    print('Pocket files moved')

    for i in crd_temp:
        shutil.copy(i,bindingsite_dir)
    print('Template files moved')

    count=0
    SV_dir = os.path.join(FLAPP_path,'SiteVector')

    while os.path.exists(SV_dir):
        count+=1
        SV_dir=os.path.join(FLAPP_path,'SiteVector'+str(count))
    os.mkdir(SV_dir)
    
    BS_dir_asvar=bindingsite_dir.split('/')[-1]
    SV_dir_asvar=SV_dir.split('/')[-1]
    
    os.chdir(FLAPP_path)

    print('Running FLAPP...')
    
    os.system('conda run -n FLAPP python3.9 CreateVectors.py ' + BS_dir_asvar + SV_dir_asvar)
    os.system('conda run -n FLAPP python3.9 FLAPP.py ' + SV_dir_asvar + ' Pairs.txt Outfile.txt 4')

    return

if __name__ == "__main__":
    if len(sys.argv) ==3:
        af_dir = sys.argv[1] # path to pocket PDB files
        crd_dir = sys.argv[2] # path to CRD directory
        
        main(af_dir,crd_dir)
    else:
        print("Usage: python automate_FLAPP.py <path to pocket PDB files> <path to SiteDesign directory>")




    

        