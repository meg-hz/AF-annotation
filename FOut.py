# FOut.py
# single hashes for info. double hashes for debugging/alternate code
# parsing of FLAPP outfiles

# IMPORTANT: this code only works if it is in the same directory as the FLAPP directory

import shutil, os, sys

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
        write_tsv.write("Accession Number\tPocket File\tNumber of Occurences\n")
        for line in filelines[1:]:
            if line.split()[0]== pocket_file:
                value+=1
            else:
                write_tsv.write(f'{protein}\t{pocket_file}\t{value}\n')
                pocket_file=line.split()[0]
                protein=pocket_file.split('-')[1]
                value=1
        print(f'{os.path.split(tsvname)[-1]} generated')

    
if __name__ == "__main__":

    if len(sys.argv)==4 and sys.argv[1]=='--fmin':
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
        print("FOut.py → USAGE:")
        print("To get entries within a particular fmin cutoff:")
        print("\tpython FOut.py --fmin <path to alignment file> <cutoff>\t OR")
        print("\tpython FOut.py --fmin <path to alignment file> -all\n")
        print("To get analysis reports of threshold outfiles:")
        print("\tpython FOut.py --fset <path to threshold file>\t OR")
        print("\tpython FOut.py --fset <path to folder with threshold files> -all\n")
