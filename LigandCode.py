# LigandCode.py
# get the corresponding Ligand names for the given list of 3 letter PDB codes
# single hashes for info. double hashes for debugging/alternate code

import requests, sys, os
from bs4 import BeautifulSoup

def get_ligand(molcode):

    molcode=str(molcode).strip() #removing whitespaces 

    response = requests.get(f'https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{molcode[:3]}')
    data = response.text
    soup= BeautifulSoup(data, 'html.parser')

    ref=soup.get_text().split("\n")
    ref=[i.strip() for i in ref if i.strip() !='']
    mol= [i for i in ref if ref[ref.index(i)-1]=='Molecule name']
    mol_type= [i for i in ref if ref[ref.index(i)-1]=='Type description']

    if mol and mol_type:
        return f"{molcode}\t{mol[0].title()}\t{mol_type[0].title()}"
    else:
        return f"{molcode}\tN/A\tN/A"

def main(filename,outpath=None):
    with open(filename,'r') as readfile:
        lines=readfile.readlines()

    count=0
    if outpath == None:
        outpath=os.getcwd()

    outname=os.path.split(filename)[-1][:-4]+ "_ligands"

    outfile=os.path.join(outpath,f'{outname}.txt')
    while os.path.isfile(outfile):
        count+=1
        outfile=outname+str(count)+'.txt'

    count=0
    print("Writing File...")
    with open(outfile,'w') as writefile:
        writefile.write(f"Ligand Code\tLigand Name\tAdditional Info\n")
        for i in lines:
            count+=1
            writefile.write(f"{get_ligand(i)}\n")
        
        if count%50==0:
            print(f"{count} ligands written")
 
    print("Complete!")

if __name__=="__main__":
    if len(sys.argv) in (2,3):
        listpath=sys.argv[1]
        outpath=None
        if not os.path.isfile(listpath):
            print("Invalid List File path. Exiting...")
            exit(1)
        elif len(sys.argv)==3:
            outpath=sys.argv[2]
            if not os.path.isdir(outpath):
                print("Invalid Output Directory Path. Exiting...")
                exit(1)
        main(listpath,outpath)

    else:
        print("LigandCode.py → USAGE:")
        print("\tLigandCode.py <path to list file>")
        print("\tLigandCode.py <path to list file> <path to output dir>")
