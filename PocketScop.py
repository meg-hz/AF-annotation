# PocketScop.py
# get the family assignment for pockets using the SUPERFAM assignment file generated for the original protein files
# single hashes for info. double hashes for debugging/alternate code

import os, sys

def pocket_scop(folderpath,reffile):
    print("Reading reference files")
    pocket_dict={}
    for file in os.listdir(folderpath):
        filepath= os.path.join(folderpath,file)

        with open(filepath,'r') as pocketfilepath:
            lines=pocketfilepath.readlines()
            lines=sorted(list(set([str(x[17:26].strip()) for x in lines if x.startswith('ATOM')])))
            # ^^ get residue, chain name and position from pocket PDB file

        accession=file.split('-')[1].strip()    # gettting accession number from file name

        pocket_dict[accession]=lines

    with open(reffile,'r') as referencepath:
        lines= referencepath.readlines()
    reference= [x for x in lines]

    # generating a new output file
    basename='Superfam_forPockets.txt'
    outfile=os.path.join(os.getcwd(),basename)
    count=0
    while os.path.isfile(outfile):
        count+=1
        outfile=os.path.join(os.getcwd(),f"{basename[:-4]}{count}.txt")


    print(f"Writing {os.path.split(outfile)[-1]}...")

    with open(outfile,'w') as f:
        # writing header
        f.write(f"Protein Accession\tSuperfamily Name\tSuperfam ID\t")
        f.write(f"Fraction of Residues\tAssignment Notation\tAll Residues\n")

        
        for protein in pocket_dict:
            scops=[x for x in reference[1:] if protein== x.split('\t')[0]]
            for assignment in scops:
                arr=assignment.split('\t')
                start= int(arr[1])
                end= int(arr[2])
                res_set=[]
                res_not=''    

                for residue in pocket_dict[protein]:
                    if start < int(residue[-4:]) < end:
                        res_set.append(residue)
                        res_not+='x'
                    else:
                        res_not+='-'

                if res_set:
                    res_ratio=f"{len(res_set)}/{len(pocket_dict[protein])}"
                    f.write(f"{protein}\t{arr[4]}\t{arr[5].strip()}\t{res_ratio}\t{res_not}\t{res_set}\n")

    print("Writing Complete!")

if __name__ == "__main__":
    if len(sys.argv)==3: 
        dir_to_pocket = sys.argv[1]
        protein_scopFile = sys.argv[2]

        if not (os.path.isdir(dir_to_pocket) and os.path.isfile(protein_scopFile)):
            print("Invalid paths. Exiting...")
            exit(1)

        elif not os.path.isdir(dir_to_pocket):
            print("Invalid Directory Path. Exiting...")
            exit(1)
        
        elif not (os.path.isfile(protein_scopFile) and protein_scopFile.endswith('.txt')):
            print("Invalid SUPERFAM File Path. Exiting...")
            exit(1)

        pocket_scop(dir_to_pocket,protein_scopFile)
    
    else:
        print("PocketScop.py → USAGE:")
        print("python PocketScop.py <path to directory with pocket files> <path to SCOPGet.py output file>")
