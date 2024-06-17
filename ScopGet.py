# ScopGet.py
# To generate input file SUPFAM web tool and to analyse output
# single hashes for info. double hashes for debugging/alternate code

import os, sys, PDB_modules as PDBm

# generate input combined fasta file for SUPERFAM
def fasta_for_spf(mode,in_dir,out_dir=os.getcwd(),ref_file=None):
  
    # getting list of all PDB files in the folder
    print("Reading Directory...")

    pdb_list=[]
    for PDBfile in os.listdir(in_dir):
        PDB_filepath=os.path.join(in_dir,PDBfile)
        if os.path.isfile(PDB_filepath) and PDBfile.endswith('.pdb'):
            pdb_list.append(PDB_filepath)
                
    ##print(pdb_list)        

    # getting list of PDB files for which fasta needs to be generated
    listoffiles=[]  
    if mode=='--all':
        listoffiles=pdb_list

    elif mode=='--list':
        print("Reading Reference File")
        with open(ref_file,'r') as fastalist:
            filelist=[i.strip() for i in fastalist.readlines()]

        for ref in pdb_list:
            code1= os.path.split(ref)[-1].split('-')[1].strip()
            for code2 in filelist:
                if code1==code2:
                    ##print(f"{code1}??{code2}")
                    listoffiles.append(ref)
            
    ##print(listoffiles)
    print(f"{len(listoffiles)} files considered")

    # iterating through PDB file paths to write fasta files
    count=0
    out_file_name='SPF_input_1p1.fa'
    if os.path.exists(os.path.join(out_dir,out_file_name)):
        count+=1
        out_file_name='SPF_input'+str(count)+'p1.fa'
    out_file_path=os.path.join(out_dir,out_file_name)

    print("Writing Fasta file...")
    with open(out_file_path, "w") as f:
        line_count=0
        file_count=1
        for file_path in listoffiles:
            f.write(PDBm.get_fasta(file_path))
            line_count += 1
            if line_count >= 1000:
                # Close current file and open a new one
                f.close()
                file_count+=1
                line_count=0

                print("1000 Fastas exceeded. Generating new pdb file")
                out_file_name='SPF_input'+str(count)+'p'+str(file_count)+'.fa'
                out_file_path=os.path.join(out_dir,out_file_name)

                f = open(out_file_path, "w")
    
    print("Fasta Files written")


# parsing SUPERFAM output    
def ScopParse(inputSCOP,spfREF,outpath=None):    

    print("Reading reference files")

    with open(inputSCOP,'r') as scopIN:     # SUPERFAM output file
        ref1= scopIN.readlines()

    with open(spfREF, 'r') as spfIN:        # file containing family names
        ref2= spfIN.readlines()

    dict={} # to store all protein names with assignment and corresponding residues
    SPFdict={} # to store all family names for corresponding SUPFAM ID numbers

    # running through SUPERFAM o/p
    for line in ref1:
        protein = line.split()[0]   # accession number
        if protein not in dict:
            dict[protein]={}

        scop_ref = int(line.split()[-1])    # SUPERFAM id number

        # amino acid ranges for assignment
        if ',' in line.split()[2]:  
            res_range = line.split()[2].split(',')
        else:    
            res_range = [line.split()[2]]

        dict[protein][scop_ref]=res_range

        # running through SUPERFAM assignment file
        for line in ref2[5:]:    
            if scop_ref == int(line.split('\t')[8]): # SUPERFAM id number
                if scop_ref not in SPFdict:
                    SPFdict[scop_ref]= [str(line.split('\t')[6]),str(line.split('\t')[-2])]
                    # ^^ for the superfamily ID number, we get the family and superfamily name in a dictionary

    print("Dictionaries completed!")
        
    # dict= {acession: {SUPERFAM ID1: [aa range(s)], SUPERFAM ID2: [aa range(s)] } }
    # SPFdict= {SUPERFAM ID1: [Superfamily name, Family name], SUPERFAM ID2: [Superfamily name, Family name]}

    if outpath==None:
        outpath=os.getcwd()

    outfile=os.path.join(outpath,'Supfam_Assignment.txt')
    count=0
    while os.path.isfile(outfile):
        count+=1
        outfile=f'{outfile[:-4]}({count}).txt'

    print(f"Writing {outfile}...")

    ##print(SPFdict)
    with open(outfile,'w') as out_f:
        out_f.write("Protein Accession\tStart Res\tEnd Res\t")
        out_f.write("Protein Family\tProtein Superfamily\tSuperfam ID\n")

        for protein in dict: 
            for one_scop in dict[protein]:
                for res_range in dict[protein][one_scop]:
                    if one_scop in SPFdict:
                        out_f.write(f"{protein}\t{res_range.split('-')[0]}\t{res_range.split('-')[1]}\t")
                        out_f.write(f"{SPFdict[one_scop][0]}\t{SPFdict[one_scop][-1]}\t{one_scop}\n")
                    else:
                        out_f.write(f"{protein}\t{res_range.split('-')[0]}\t{res_range.split('-')[1]}\tNA\tNA\t{one_scop}\n")
         
    print("Completed!\n")

if __name__=="__main__":

    # for input for SUPFAM web tool
    if len(sys.argv) in (4,5,6) and sys.argv[1]=='--input':
        mode=sys.argv[2]
        input_dir=sys.argv[3]

        # accounting for output dir
        if len(sys.argv) in (4,5) and mode=='--all':
            file_list=None
            if len(sys.argv)==4:
                output_dir=None

            elif not os.path.isdir(sys.argv[4]):
                print("Invalid Output Path. Exiting Code...")
                exit(1)

            else:
                output_dir=sys.argv[4]

        elif len(sys.argv) in (5,6) and mode=='--list':
            file_list=sys.argv[4]
            if not (os.path.isfile(file_list) and file_list.endswith('.txt')):
                print("Invalid text file path. Exiting Code...")
                exit(1)

            elif len(sys.argv)==5:
                output_dir=None

            elif not os.path.isdir(sys.argv[5]):
                print("Invalid Output Path. Exiting Code...")
                exit(1)
            else:
                output_dir=sys.argv[5]   
    
        else: 
            print("Syntax Error. Exiting Code...")
            exit(1)

        # accounting for validity of inputpath
        if not os.path.isdir(input_dir):
            print("Invalid input folder. Exiting Code...")
            exit(1)

        fasta_for_spf(mode,input_dir,output_dir,file_list)

    # for output for SUPFAM web tool
    elif len(sys.argv) in (4,5) and sys.argv[1]=='--output':
        spf_file=sys.argv[2]
        ref_file=sys.argv[3]

        for i in (spf_file,ref_file):
            if not (os.path.isfile(i) and i.endswith('.txt')):
                print('Invalid File Path')
                exit(1)
        
        if len(sys.argv)==4:
            output_dir=None
        
        elif not os.path.isdir(sys.argv[4]):
            print("Invalid Output Path. Exiting Code...")
            exit(1)
        
        else:
           output_dir=sys.argv[4]
        
        ScopParse(spf_file,ref_file,output_dir)

    else:
        print("ScopGet.py → USAGE:")
        print("To Get Fasta file input for Superfamily Web Tool ")
        print("\tFor all files in a folder:")
        print("\t\tpython ScopGet.py --input --all <input dir>\t OR")
        print("\t\tpython ScopGet.py --input --all <input dir> <output dir>")
        print("\tFor a list of files:")
        print("\t\tpython ScopGet.py --input --list <input dir> <text file containing list>\t OR")
        print("\t\tpython ScopGet.py --input --list <input dir> <text file containing list> <output dir>")
        print()
        print("To Parse Superfamily Web Tool Output ")
        print("\t\tpython ScopGet.py --output <SUPERFAM output file> <reference file>\t OR")
        print("\t\tpython ScopGet.py --output <SUPERFAM output file> <reference file> <output dir>")
        print() 
