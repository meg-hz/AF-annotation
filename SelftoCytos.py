# SelftoCytos.py
# single hashes for info. double hashes for debugging/alternate code
# to generate a cytoscape input file from a given alignment file

import os, sys

def cytoscape_input(in_filepath, out_dirpath=None):
    all_pairs=[]
    with open (in_filepath,"r") as f:
        for line in f:
            pock_arr=list(set(sorted([i for i in line.split()[:2]])))
            if len(pock_arr) ==2 and pock_arr not in all_pairs:
                all_pairs.append(pock_arr)
    
    if out_dirpath==None:
        out_dirpath=os.path.dirname(in_filepath)
        
    out_name=os.path.split(in_filepath)[-1][:-4]+"_out.csv"
    out_filepath=os.path.join(out_dirpath,out_name)
    
    count=0
    while os.path.isfile(out_filepath):
        count+=1
        out_name=out_name[:-4]+str(count)+".csv"
        out_filepath=os.path.join(out_dirpath,out_name)


    with open(out_filepath,"w") as f2:
        f2.write(f"Source, Target\n")
        for couple in all_pairs[1:]: 
            ##print("hi")
            f2.write(f"{couple[0]},{couple[1]}\n")

if __name__=="__main__":
    if len(sys.argv) in (2,3):
        cutoff_file=sys.argv[1]
        if not os.path.isfile(cutoff_file) or not cutoff_file.endswith('.txt'):
            print("Invalid Outfile")
            sys.exit(1)
        elif len(sys.argv)==2:
            op_file_dir=None
            
        elif len(sys.argv)==3:
            op_file_dir=sys.argv[2]
            if not os.path.isdir(op_file_dir):
                print("Invalid Output Path")
                sys.exit(1)

        cytoscape_input(cutoff_file,op_file_dir)  
    else:
        print("SelfToCytos.py → USAGE:")
        print("\tpython SelfToCytos.py <path to alignment file>\t OR")
        print("\tpython SelfToCytos.py <path to alignment file> <path to output dir>\t OR")
        
