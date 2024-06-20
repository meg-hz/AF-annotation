# FLigand.py
# obtain all unique ligands for a protein from an alignment file
# single hashes for info. double hashes for debugging/alternate code

import os, sys, requests
from bs4 import BeautifulSoup

def get_ligand(molcode):

    molcode=str(molcode).strip() #removing whitespaces 

    response = requests.get(f'https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{molcode[:3]}')
    data = response.text
    soup= BeautifulSoup(data,'html.parser')

    ref=soup.get_text().split("\n")
    ref=[i.strip() for i in ref if i.strip() !='']
    mol= [i for i in ref if ref[ref.index(i)-1]=='Molecule name']

    if mol:
        return mol[0].title()
    else:
        return 'N/A'

def findLigands(infile,outpath=None):
    with open(infile,'r') as f:
        lines=f.readlines()

    pocketdict={} # dictionary for the ligands of a pocket
    pocketres={} # dictionary for the number of residues in a pocket

    for x in lines[1:]:
        our_pocket=x.split()[0]
        no_res= int(x.split()[3])

        pocketres[our_pocket]=no_res 

        if our_pocket not in pocketdict:
            pocketdict[our_pocket]={}
    
        temp_pocket=x.split()[1].split('_')[1].strip() 
        pock_ratio=int(x.split()[4])

        count=1
        while temp_pocket in pocketdict[our_pocket]:
            count+=1 
            temp_pocket=x.split()[1].split('_')[1].strip() + "("+str(count)+")"
        
        pocketdict[our_pocket][temp_pocket]=pock_ratio

    filename=os.path.split(infile)[-1][:-4]+'_Ligands.txt'
    if outpath==None:
        outfile=os.path.join(os.getcwd(),filename)
    else:
        outfile=os.path.join(outpath,filename)
    
    count=0
    while os.path.isfile(outfile):
        count+=1
        outfile=outfile[:-4]+str(count)+'.txt'
            
    with open(outfile,'w') as f:
        f.write('Protein\tPocketFile\tNo. of Residues\tNo of Pockets\tLigand Code\tLigand Name\t No. of BS Residues\tRatio of residues of Prot Pocket/BS Pocket\n')
        for i in pocketdict:
            ligand_dict=pocketdict[i]
            protein=i.split('-')[1]
            for ligand in ligand_dict:
                f.write(f'{protein}\t{i}\t{pocketres[i]}\t{len(ligand_dict)}\t')
                f.write(f'{ligand}\t{get_ligand(ligand)}\t{ligand_dict[ligand]}\t{pocketres[i]/ligand_dict[ligand]}\n')
        
    print(f'{os.path.split(outfile)[-1]} generated.')

if __name__ == "__main__":
    if len(sys.argv) in (2,3):
        inpath=sys.argv[1]
        if not os.path.isfile(inpath):
            print("Invalid Input file path. Exiting...")
            exit(1)

        elif len(sys.argv) ==2:
            outpath=None

        elif len(sys.argv) ==3:
            outpath=sys.argv[2]
            if not os.path.isdir(outpath):
                print("Invalid Output Dir path. Exiting...")
                exit(1)
            
        findLigands(inpath,outpath)
    
    else:
        print("FLigand.py → USAGE:")
        print("\tFLigand.py <path to alignment/threshold file>")
        print("\tFLigand.py <path to alignment/threshold file> <path to output dir>")
