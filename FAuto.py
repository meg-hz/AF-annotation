# FAuto.py
# single hashes for info. double hashes for debugging/alternate code
# move files to flapp folder + generate pairs.txt + run FLAPP

# IMPORTANT: this code only works if it is in the same directory as the FLAPP directory

import os, shutil, sys

root= os.getcwd()
FLAPP_path=os.path.join(root,'FLAPP') 

#---------------------------------------------------------------------------------------------------
# MODULES FOR RUNNING FLAPP

# obtaining list of addresses of pocket files 
def input_add(ip_path):
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
def template_add(site2_path):
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
    count=0
    while os.path.exists(filepath):
        count+=1
        filepath = os.path.join(FLAPP_path,'Pairs'+str(count)+'.txt')
    
    with open(filepath,'w') as pairwise:
        count=0
        for file1 in list1:
            count+=1
            count1=0
            file1_name=os.path.split(file1)[-1]
            for file2 in list2:
                file2_name=os.path.split(file2)[-1]
                pairwise.write(file1_name + '\t' + file2_name +'\n')   
                count1+=1
            print(f" List2: {count1}")
        print(f" List1: {count}")

    print(f"{os.path.split(filepath)[-1]} generated")
    return filepath

# Making a New Folder with all binding sites
def FLAPP_BS():
    count=0
    BS_dir=os.path.join(FLAPP_path,'BindingSites')
    while os.path.exists(BS_dir):
        count+=1
        BS_dir=os.path.join(FLAPP_path,'BindingSites'+str(count))
    os.mkdir(BS_dir)
    return BS_dir 

def main(mode,pocket_path, template_path, cores):

    # get the file addresses
    if mode == '--f2temp':
        af_temp=input_add(pocket_path)
        template_temp=template_add(template_path)
    elif mode == '--f2self':
        af_temp=input_add(pocket_path)
        template_temp=input_add(template_path)
    
    elif mode == '--f2misc':
        af_temp=template_add(template_path)
        template_temp=template_add(template_path)

    pair_file=os.path.split(write_pairs(af_temp,template_temp))[-1] #pairs.txt

    bindingsite_dir=FLAPP_BS()
    print(f"{bindingsite_dir.split('/')[-1]} Folder Generated")

    for i in af_temp:
        shutil.copy(i,bindingsite_dir)
    print('Pocket files moved')

    for i in template_temp:
        shutil.copy(i,bindingsite_dir)
    print('Template files moved')

    count=0
    SV_dir = os.path.join(FLAPP_path,'SiteVector')
    while os.path.exists(SV_dir):
        count+=1
        SV_dir=os.path.join(FLAPP_path,'SiteVector'+str(count))
    os.mkdir(SV_dir)

    BS_dir_asvar=os.path.split(bindingsite_dir)[-1]
    SV_dir_asvar=os.path.split(SV_dir)[-1]

    count=0
    outfile_path = os.path.join(FLAPP_path,'Outfile.txt')
    while os.path.exists(outfile_path):
        count+=1
        outfile_path=os.path.join(FLAPP_path,'Outfile'+str(count)+'.txt')
    outfile_path=os.path.split(outfile_path)[-1]
    
    os.chdir(FLAPP_path)
    print("Running FLAPP...")
    print('Generating Binary Files...')
    os.system(f'conda run -n FLAPP')
    os.system(f'python CreateVectors.py {BS_dir_asvar} {SV_dir_asvar}')
    print('Generating Alignment File...')
    os.system(f'python FLAPP.py {SV_dir_asvar} {pair_file} {outfile_path} {cores}')

    return

#---------------------------------------------------------------------------------------------------
# MODULES FOR PARSING THE OUTFILES

def Fmin_cutoff(outfile,threshold):
    outfile_name=os.path.split(outfile)[-1]
    file_name= outfile_name[:-4]+'_cutoff'+str(threshold)+'.txt'

    with open(file_name,'w') as write_file:
        # Read lines from input file
        with open(outfile, 'r') as read_file:
            lines = read_file.readlines()

            write_file.write(lines[0]) #writing the header line

            for line in lines[1:]:
                #ignoring lines where there's no alignment
                if line.split()[-2]=='NoResidueMatch':
                    continue 

                # Extract the Fmin (6th component) in terms of %
                value = int(float(line.split()[5])*100)
                if value >= threshold: # Check if the value exceeds the threshold
                    write_file.write(line)

    print(f"{os.path.split(file_name)[-1]} generated.")
    return file_name

def Fmin_all(outfile):
    outdir=outfile[:-4]+str("_all")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    count=0
    
    for threshold in (40,50,60,70,80,90):
        print(f"Generating threshold files for {threshold}...")
        shutil.move(Fmin_cutoff(outfile,threshold),outdir)
        count+=1
    
    print(f"All files moved to {os.path.split(outdir)}")
    return
             
def cutoff_analysis_set(outfolder):
    
    for file in os.listdir(outfolder):
        print(f"Reading {file}...")
        filepath=os.path.join(outfolder,file)
        cutoff_analysis(filepath)

    print('\nAll analysis files generated')

def cutoff_analysis(filepath):
    print("Reading File...")
    with open(filepath,'r') as read_file:
        filelines= read_file.readlines()
            
    pocket_file = filelines[1].split()[0] #name of the pocket file
    protein=pocket_file.split('-')[1] #accession number of protein
    value=0 #number of occurances

    tsvname=filepath[:-4]+'_analysis.tsv'
    with open(tsvname,'w')as write_tsv:
        write_tsv.write("Accession Number\tPocket File\tNumber of Occurances\n")
        for line in filelines[1:]:
            if line.split()[0]== pocket_file:
                value+=1
            else:
                write_tsv.write(f'{protein}\t{pocket_file}\t{value}\n')
                pocket_file=line.split()[0]
                protein=pocket_file.split('-')[1]
                value=1
        print(f'{os.path.split(tsvname)[-1]} generated')


#---------------------------------------------------------------------------------------------------
# TERMINAL PARSING

if __name__ == "__main__":
    if len(sys.argv) ==5 and sys.argv[1] in ['--f2temp','--f2misc']:
        mode=sys.argv[1]
        af_dir = sys.argv[2] # path to pocket PDB files
        temp_dir = sys.argv[3] # path to template directory
        core_no= sys.argv[4] #no. of cores
        
        if not (os.path.isdir(af_dir) and os.path.isdir(temp_dir)):
            print("Invalid path name(s). Exiting Code.")
            sys.exit(1)
    
        print("Running FLAPP...")
        main(mode,af_dir,temp_dir,core_no)
    
    elif len(sys.argv)==4 and sys.argv[1] in ['--f2self','--f2misc']:
        mode=sys.argv[1]
        af_dir = sys.argv[2]
        core_no= sys.argv[3]
    
        if not os.path.isdir(af_dir):
            print("Invalid path name(s). Exiting Code.")
            sys.exit(1)

        print("Running FLAPP...")
        main(mode,af_dir,af_dir,core_no)

    elif len(sys.argv)==4 and sys.argv[1]=='--fmin':
        outfile_path = sys.argv[2]
        cutoff= sys.argv[3]

        if not (os.path.isfile(outfile_path) and outfile_path.endswith('.txt')):
            print("Invalid file path. Exiting Code.")
            sys.exit(1)

        elif not (cutoff.isdigit() and 0<int(cutoff)<100):
            if cutoff== "-all":
              Fmin_all(outfile_path)
              sys.exit(0)

            print("Invalid cutoff Value. Exiting Code")
            sys.exit(1)
        
        Fmin_cutoff(outfile_path,int(cutoff))  
    
    elif sys.argv[1]=='--fset':
        if len(sys.argv)==3: #in case of single file
            single_file=sys.argv[2]

            if not (os.path.isfile(single_file) and single_file.endswith('.txt')):
                print("Invalid file path. Exiting Code.")
                sys.exit(1)
            
            cutoff_analysis(single_file)
        
        elif len(sys.argv)==4 and sys.argv[3]=='-all':
            folderr=sys.argv[2]

            if not os.path.isdir(folderr):
                print("Invalid file path. Exiting Code.")
                sys.exit(1) 
            
            cutoff_analysis_set(folderr)         

    else:
        print("FAuto.py → USAGE:")
        print("To run FLAPP:")
        print("\tpython FAuto.py --f2temp <path to pocket-file-dir> <path to template directory> <no of cores>\t OR")
        print("\tpython FAuto.py --f2self <path to pocket-file-dir> <no of cores>\n")
        print("To get entries within a particular fmin cutoff:")
        print("\tpython FAuto.py --fmin <path to alignment file> <cutoff>\t OR")
        print("\tpython FAuto.py --fmin <path to alignment file> -all\n")
        print("To get analysis reports of threshold outfiles:")
        print("\tpython FAuto.py --fset <path to alignment file>\t OR")
        print("\tpython FAuto.py --fset <path to alignment folder> -all\n")
