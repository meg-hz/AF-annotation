import os, sys, PDB_modules as PDBm

def fasta_for_spf(mode,in_dir,out_dir,ref_file=None):

    listoffiles=[]

    print("Reading Reference File")
    with open(ref_file,'r') as fastalist:
        filelist=[i.strip() for i in fastalist.readlines()]
    
    ##print(filelist)
    
    print("Reading Directory...")
    # getting list of all valid PDB files
    pdb_list=[]
    for PDBfile in os.listdir(in_dir):
        PDB_filepath=os.path.join(in_dir,PDBfile)
        if os.path.isfile(PDB_filepath) and PDBfile.startswith('AF-')and PDBfile.endswith('.pdb'):
            pdb_list.append(PDB_filepath)
                
    ##print(pdb_list)        

    # getting list of PDB files for which fasta needs to be generated
    if mode=='-all':
        listoffiles=pdb_list
    elif mode=='--list':
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
        for file_path in listoffiles:
            f.write(PDBm.get_fasta(file_path))
            line_count += 1
            if line_count >= 1000:
                # Close current file and open a new one
                f.close()
                file_count += 1
                line_count = 0

                print("1000 Fastas exceeded. Generating new pdb file")
                out_file_name='SPF_input'+str(count)+'p'+str(file_count)+'.fa'
                out_file_path=os.path.join(out_dir,out_file_name)

                f = open(out_file_path, "w")
    
    print("Fasta Files written")


def spf_out(spf,outdir=None):
    if outdir==None:
        outdir = os.path.dirname(spf)

    with open(spf,"r") as spf_file:
        SPFoutput = [line for line in spf_file.readlines()]
    
    SPF_dict={}
    for line in SPFoutput:
        SPFarr=[int(i) for i in line.replace(',', ' ').split()[1:] if i != '_gap_' ]
        SPFarr= sorted(list(set(SPFarr)))
        SPF_dict[line.split('\t')[0]]=SPFarr

    SPF_dict=dict(sorted(SPF_dict.items()))

    count=0
    SPF_out_name='SPF_output.txt'
    if os.path.exists(os.path.join(outdir,SPF_out_name)):
        count+=1
        out_file_name='SPF_input'+str(count)+'.txt'
    out_file_path=os.path.join(outdir,SPF_out_name)

    with open(out_file_path,'w') as file:
        file.write("Accession Number\tSuperFam ID\n")
        for protein in SPF_dict:
            file.write(f"{protein}\t{SPF_dict[protein]}\n")



if __name__=="__main__":
    if len(sys.argv) in (3,4) and sys.argv[1]=='--all':
        mode = sys.argv[1]
        input_dir=sys.argv[2]
        if not os.path.isdir(input_dir):
            print("Invalid Input folder path. Exiting Code.")
            exit(1)
        elif len(sys.argv) ==3:
            fasta_for_spf(mode,input_dir,input_dir)
        elif len(sys.argv) ==4:
            output_dir=sys.argv[3]
            if not os.path.isdir(output_dir):
                print("Invalid Output folder path. Exiting Code.")
                exit(1)
            fasta_for_spf(mode,input_dir,output_dir)

    elif len(sys.argv) in (4,5) and sys.argv[1]=='--list':
        mode = sys.argv[1]
        input_dir=sys.argv[2]
        file_list=sys.argv[3]

        if not os.path.isdir(input_dir):
            print("Invalid Input folder path. Exiting Code.")
            exit(1)
        elif not os.path.isfile(file_list) and file_list.endswith('.txt'):
            print("Invalid text file path. Exiting Code.")
            exit(1)

        elif len(sys.argv) ==4:
            output_dir=input_dir
        elif len(sys.argv) ==5:
            output_dir=sys.argv[4]

        if not os.path.isdir(output_dir):
            print("Invalid Output folder path. Exiting Code.")
            exit(1)
        
        fasta_for_spf(mode,input_dir,output_dir,file_list)
    
    elif len(sys.argv)in (3,4) and sys.argv[1]=='--spf':
        input_file=sys.argv[2]
        if not os.path.isfile(input_file) and input_file.endswith('.txt'):
            print("Invalid text file path. Exiting Code.")
            exit(1)
        elif len(sys.argv) ==3:
            spf_out(input_file)
        elif len(sys.argv) ==4:
            output_dir=sys.argv[3]
            if not os.path.isdir(output_dir):
                print("Invalid Output folder path. Exiting Code.")
                exit(1)
            spf_out(input_file,output_dir)

    else:
        print("ScopGet.py → USAGE:")
        print("1. To Get Fasta file input for Superfamily Web Tool:")
        print("For all files in a folder")
        print("\tpython ScopGet.py --all <input dir>\t OR")
        print("\tpython ScopGet.py --all <input dir> <output dir>")
        print("For a list of files:")
        print("\tpython ScopGet.py --list <input dir> <text file containing list>\t OR")
        print("\tpython ScopGet.py --list <input dir> <text file containing list> <output dir>")
        print("\n2. To analyse output from Superfamily Web Tool")
        print("\tpython ScopGet.py --spf <input file>\t OR")
        print("\tpython ScopGet.py --spf <input file> <output dir>\t")
