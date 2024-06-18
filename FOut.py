# FOut.py
# parsing of FLAPP outfiles and providing cytoscape input
# single hashes for info. double hashes for debugging/alternate code

# IMPORTANT: this code only works if it is in the same directory as the FLAPP outfiles

import shutil, os, sys

### THRESHOLD FILES ###
# get a set of output files for inputted level of Fmin cutoff
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
                ##print(line.split())
                fmin = int(float(line.split()[5])*100)
                fmax = int(float(line.split()[6])*100)
                if fmin >= threshold and fmax<100: # Check if the value exceeds the threshold
                    write_file.write(line)

    print(f"{os.path.split(file_name)[-1]} generated.")
    return file_name

# get a set of output files for different levels of Fmin cutoff
def Fmin_all(outfile):
    outdir=outfile[:-4]+str("_all")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    count=0
    
    for threshold in (40,50,60,70,80,90):
        print(f"Generating threshold files for {threshold}...")
        shutil.move(Fmin_cutoff(outfile,threshold),outdir)
        count+=1
    
    print(f"All files moved to {os.path.split(outdir)[1]}")
    return

### FREQUENCY ANALYSIS FILES ###
# to get a frequency analysis file from a threshold outfile
def cutoff_analysis(filepath):
    print("Reading File...")
    with open(filepath,'r') as read_file:
        filelines= read_file.readlines()

    pocketlist=[i.split()[0] for i in filelines[1:]]
    pocketdict={}

    for i in pocketlist:
        if i not in pocketdict:
            pocketdict[i]=1
        else:
            pocketdict[i]+=1 

    tsvname=filepath[:-4]+'_analysis.tsv'
    with open(tsvname,'w')as write_tsv:
        write_tsv.write("Accession Number\tPocket File\tNumber of Occurences\n")
        for pocket in pocketdict:
            protein=pocket.split('-')[1]
            write_tsv.write(f'{protein}\t{pocket}\t{pocketdict[pocket]}\n')
        print(f'{os.path.split(tsvname)[-1]} generated')

# to get a frequency analysis file from a set of threshold outfiles
def cutoff_analysis_set(outfolder):
    
    for file in os.listdir(outfolder):
        filepath=os.path.join(outfolder,file)
        if filepath.endswith('.txt') and '_cutoff' in filepath: 
            print(f"Reading {file}...")
            cutoff_analysis(filepath)

    print('\nAll analysis files generated')

### CYTOSCAPE INPUT FILES ###
# to get a cytoscape input file from a threshold outfile
def cytoscape_input(in_filepath):
    all_pairs=[]
    with open (in_filepath,"r") as f:
        for line in f:
            pock_arr=list(set(sorted([i for i in line.split()[:2]])))
            if len(pock_arr) ==2 and pock_arr not in all_pairs:
                all_pairs.append(pock_arr)
        
    out_filepath=in_filepath[:-4]+"_out.csv"
    count=0
    while os.path.isfile(out_filepath):
        count+=1
        out_filepath=out_filepath[:-4]+str(count)+".csv"

    with open(out_filepath,"w") as f2:
        f2.write(f"Source, Target\n")
        for couple in all_pairs[1:]: 
            ##print("hi")
            f2.write(f"{couple[0]},{couple[1]}\n")
    print(f'{os.path.split(out_filepath)[-1]} generated')

# to get a cytoscape input file from a set of threshold outfiles
def cytoscape_input_set(outfolder):
    for file in os.listdir(outfolder):
        filepath=os.path.join(outfolder,file)
        if filepath.endswith('.txt') and '_cutoff' in filepath:
            print(f"Reading {file}...")
            cutoff_analysis(filepath)

    print('\nAll analysis files generated')

### SYSTEM PARSING ###    
if __name__ == "__main__":

    if len(sys.argv)==4 and sys.argv[1]=='--fmin':
        cutoff= sys.argv[2]
        outfile_path = sys.argv[3]

        if not (os.path.isfile(outfile_path) and outfile_path.endswith('.txt')):
            print("Invalid file path. Exiting Code.")
            sys.exit(1)

        elif not (cutoff.isdigit() and 0<int(cutoff)<100):
            if cutoff== "-all":
                Fmin_all(outfile_path)
                exit()

            print("Invalid cutoff Value. Exiting Code")
            exit(1)
        
        Fmin_cutoff(outfile_path,int(cutoff))  
    
    elif len(sys.argv) in (3,4) and sys.argv[1]=='--fset':

        threshold_f=sys.argv[-1]

        if len(sys.argv)==3:    # in case of single file
            if not (os.path.isfile(threshold_f) and threshold_f.endswith('.txt')):
                print("Invalid file path. Exiting Code.")
                sys.exit(1)
            cutoff_analysis(threshold_f)
        
        elif len(sys.argv)==4 and sys.argv[2]=='-all':
            if not os.path.isdir(threshold_f):
                print("Invalid file path. Exiting Code.")
                sys.exit(1) 
            
            cutoff_analysis_set(threshold_f) 

    elif len(sys.argv) in (3,4) and sys.argv[1]== '--fcytos':
        threshold_f=sys.argv[-1]

        if len(sys.argv)==3:
            if not os.path.isfile(threshold_f) or not threshold_f.endswith('.txt'):
                print("Invalid Outfile")
                sys.exit(1)
            cytoscape_input(threshold_f) 
            
        elif len(sys.argv)==4 and sys.argv[2]=='-all':
            if not os.path.isdir(threshold_f):
                print("Invalid Output Path")
                sys.exit(1)

            cytoscape_input_set(threshold_f)  

    else:
        print("FOut.py → USAGE:")
        print("To get entries within a particular fmin cutoff:")
        print("\tpython FOut.py --fmin <cutoff> <path to alignment file>\t OR")
        print("\tpython FOut.py --fmin -all <path to alignment file>\n")
        print("To generate Cytoscape input file from threshold outfiles:")
        print("\tpython FOut.py --fcytos <path to threshold file>\t OR")
        print("\tpython FOut.py --fcytos -all  <path to folder with threshold files>")
        print("To get analysis reports of threshold outfiles:")
        print("\tpython FOut.py --fset <path to threshold file>\t OR")
        print("\tpython FOut.py --fset -all <path to folder with threshold files>\n")

