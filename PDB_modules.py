# PDB_modules.py
# main function produces output csv with information about all files in alphafold folder
# single hashes for info. double hashes for debugging/alternate code

# dependencies
import re, csv, os, sys, numpy as np, matplotlib.pyplot as plt 
from collections import defaultdict 
from math import floor, log10
from textwrap import wrap

ref = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

#---------------------------------------------------------------------------------------------------
# FUNCTIONS FOR BASIC PDB FILE PARSING // FOR ALL PDB FILES

# reads file. returns an array of file lines based on the first word/index
def get_lines(path, index):    
    extracted=[]

    with open(path,"r") as f:
        pdb_read= f.readlines()

        for line in pdb_read:

            if index in line.split()[0]:
                extracted.append(line)

    if extracted:
        return extracted
    else:
        return None

# get file title
def get_title(path):
    if get_lines(path, "TITLE"):
        title=' '.join([line[10:79].strip() for line in get_lines(path, "TITLE")])
        return title.strip()
    else:
        return '-'

# get name of molecule and chain
def get_compound(path):
    if get_lines(path,'COMPND'):
        compound={}
        str=''.join([line[11:79].strip() for line in get_lines(path,'COMPND')])
        arr=str.split(";")

        for i in arr:
            if "MOLECULE:" in i: 
                key= i.split(":")[-1].strip()
            elif "CHAIN:" in i:
                chain= i.split(":")[-1].strip().split(',')
                compound[key] = [s.strip().replace(';','') for s in chain]

        if len(compound)==1:
            af_cpd=list(compound.keys())[0]
            af_ch=compound[af_cpd][0]
            return [af_cpd,af_ch]

        else:
            return compound
    else:
        return '-'

# get accession IDs of molecules in file
def get_accession(path):
    if get_lines(path,'DBREF'):
        PDB_accession=get_lines(path,'DBREF')[0][7:11]
        UP_accession= [line[33:].split()[0] for line in get_lines(path,'DBREF')]

        if (PDB_accession)=='XXXX':
            return UP_accession[0]
        else:
            return PDB_accession
    else: 
        return os.path.split(path)[-1].split('.')[0]


# get fragment number from filename
# note that this only works because Alphafold filename includes the fragment number
def get_fragment(path):
    filename = os.path.basename(path)
    fragment = re.search(r'-F(\d+)', str(filename))
    if fragment:
        return int(fragment.group(1))
    else:
        return 0
  
#---------------------------------------------------------------------------------------------------
# FUNCTIONS FOR OBTAINING SEQUENCE DATA AS FASTA FILE // FOR ALL PDB FILES
    
# returns sequence data
def get_fasta(path):

    # gets required atom data from specific type of lines.
    atom_data=[]

    for line in get_lines(path, 'ATOM'):
        res = line[17:20].strip()
        chain = line[21].strip()
        pos = int(line[22:26].strip())

        # from ATOM lines retrieve array= [residue, chain, position]
        rcp= [res, chain, pos] 

        if rcp not in atom_data:
          atom_data.append(rcp)  

    seq={} ; starting = None
    r_val=''

    for element in atom_data:
        residue=ref[element[0]]
        chain = element[1]
        position = element[2]
    
        if chain not in seq:
            seq[chain]=''
            starting = element[2]

        if position == starting:
            seq[chain]+=residue

        else:
            dashes=position-starting
            seq[chain]=seq[chain]+'-'*dashes+residue
            starting += dashes

        starting += 1

    for i in seq:
        r_val+=(f'>{get_accession(path)} | {get_title(path).split("(")[0].strip()}; Chain_{i}\n{seq[i]}\n')

    return r_val

# writes fasta output for a single file
def get_solo_fasta(pdb_path,out_path=None):

    if out_path==None:
        out_path = os.getcwd()
    filename=os.path.split(pdb_path)[-1][:-4]+".fa"
    output = os.path.join(out_path, filename)

    with open(output, "w") as fasta_file:
        fasta_file.write(get_fasta(pdb_path))


# writes out fasta for all files in the dir
def get_all_fasta(directory_path,out_path=None):
    # note that for get_solo_fasta input function before this- i/p is a file path but here it's the dir path 
    if out_path==None:
        out_path = os.getcwd()
    
    for file in os.listdir(directory_path):
        file_path=os.path.join(directory_path,file)

        if not (os.path.isfile(file_path) and file.endswith(".pdb")):
            continue

        # output file info
        filename=os.path.split(file)[-1][:-4]+".fa"
        output = os.path.join(out_path, filename)

        with open(output, "w") as fasta_file:
            fasta_file.write(get_fasta(file_path))
           

# writes fasta output for proteins with fragment PDB files
def get_combined_fasta(directory_path,out_path=None):

    if out_path==None:
        out_path =os.getcwd()

    # STEP 1: actually get the duplicates 
    def get_duplicates(directory_path):
        # returns a dictionary {accession id: no' of fragment files} for all proteins with more than one PDB file

        files_by_code = defaultdict(list)

        #STEP 1.1: get accession IDs of proteins that contain fragment files   
        for file_name in os.listdir(directory_path):
            file_path = os.path.join(directory_path, file_name)
            if os.path.isfile(file_path) and file.endswith(".pdb"):
                code=get_accession(file_path)
                files_by_code[code].append(file_name)

        ##non_duplicates={code: file_name for code,file_name in files_by_code.items() if len(file_name)==1}
        duplicates={code: file_name for code,file_name in files_by_code.items() if len(file_name)>1}

        # len(os.listdir(directory_path) is the total number of files
        # len(non_duplicates) is the number of proteins with single file entry
        # len(duplicates) is the number of proteins with multiple file entry
        # sum of len(file_names) in duplicates would be the number of fragment files

        # STEP 1.2: sort the filenames on the basis of order of fragment
        
        for file_list in duplicates.items():
            sorted(file_list, key=get_fragment(file_list))

        # STEP 1.3: return a dictionary with all accession IDs with fragments 
        return duplicates


    duplicate_dict=get_duplicates(directory_path) 

    # STEP 2: if more than one protein with fragments, make a folder where o/p fasta files will go
    if len(duplicate_dict)>1:
        out_folder=os.path.join(out_path,'fasta_fragments')
        os.makedirs(out_folder)
    else:
        out_folder=out_path

    # tracker
    folder_count=0
    file_count=defaultdict(int)

    # STEP 3: obtaining fasta file per protein/ per accession code
    for code in duplicate_dict.items():
        file_list= duplicate_dict[code]
        
        # STEP 3.1: obtain an array with sequences from all fragment files // note that this list is already ordered
        seq=[]
        for file in file_list:

            ip_file=os.path.join(directory_path,file)

            line=get_fasta(ip_file).split("\n")[1]
            seq.append(line)

            file_count[str(code)]+=1
            # for a single code, this value should be equal to the length of the fragment files in the code

        print(f"{file_count[str(code)]} files read for accession {code}")
    
        # STEP 3.2: check for overlapping regions and append to string accordingly
        dup=seq[0]
        for i in range(1,len(seq)):

            if i != (len(seq)):
                for j in range(len(seq[i])):

                    lim=len(seq[i])-j
                    if seq[i][:lim] in dup:
                        dup+=str(seq[i][lim:])
                        break
        
        # STEP 3.3: output the resultant sequence onto a file
        # Filename to fit the AF pdb naming format
        ret_file=file_list[0].split('-')
        op_filename=ret_file[0]+'-'+ret_file[1]+"-Full"+ret_file[3][:-4]+'.fa'
        output=os.path.join(out_folder,op_filename)

        with open(output,'a') as fasta_file:
            file_input=">"+str(code)+" | "+ get_title(ip_file).split("(")[0]+ "\n"+dup
            fasta_file.write(file_input)

        print(f"{op_filename} generated with {len(dup)} residues. moving on...")
        print()
        folder_count+=1

    print(f"{folder_count} files generated in total")

#---------------------------------------------------------------------------------------------------
# FUNCTIONS FOR OBTAINING CONFIDENCE SCORE RELATED DATA // FOR APLHAFOLD PDB FILES

# get an array ofper sequence confidence scores for a protein 
def get_confscore(path):
    seq_scores=[]

    for line in get_lines(path, "ATOM"):
        pos=line[22:26].strip()
        if len(seq_scores)==int(pos) - 1:
            seq_scores.append(float(line[61:66].strip()))

    return seq_scores

# get an array with total number of residues and number of residues within the corresponding confidence threshold scores
def get_coverage(path, threshold):

    ##all_conf= len(get_confscore(path))
    conf=0
            
    for res in get_confscore(path):
        conf += res > threshold

    return conf

# get a CSV with all info on the file
def get_metadata(directory_path,op_path=None):

    if op_path==None:
        op_path=directory_path
    
    csv_file_path=os.path.join(op_path,'Metadata.csv')

    files = os.listdir(directory_path)  # List all files in the folder

    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)

        csv_writer.writerow(['File Name','AF Accession No.', 'Fragment', 'File Title', 'Molecule Name', 'Chain', 
                            'Total No. Of Residues', 'Residues with 90+ Confidence', 'Residues with 70+ Confidence', 'File Size'])

        count=0

        # Iterate through each file
        for file in files:

            ##print(f"{file}")

            count+=1
            file_path = os.path.join(directory_path, file)

            # Get file data
            accession_no= get_accession(file_path)
            fragment=get_fragment(file_path)
            title= get_title(file_path).rsplit(' ', 1)[0]
            mol= get_compound(file_path)[0].strip(';')
            chain= get_compound(file_path)[1]

            res= get_coverage(file_path, 0)
            res90= get_coverage(file_path, 90)
            res70= get_coverage(file_path, 70)

            size = os.path.getsize(file_path)
        
            # Write information to CSV
            csv_writer.writerow([file, accession_no, fragment, title, mol, chain, res, res90, res70, size])

    print(count," files detected")
    print(f'Information saved to {csv_file_path}')


# get graph showing residue position vs score 
def get_conf_graph(pdb_path,out_path=None):

    if out_path==None:
        out_path = os.getcwd()
    
    arr=get_confscore(pdb_path)
    i= len(arr) 
    scale = 10 * (floor(log10(abs(i)))-1)

    # make data
    x = np.linspace(0, i, i)
    y = np.array(arr)

    # plot
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.plot(x, y, linewidth=2.0)

    title_text = get_title(pdb_path)
    wrapped_title = '\n'.join(wrap(title_text, 50))  # Wrap at 20 characters
    ax.set_title(wrapped_title)

    ax.set_xlabel('Residue Position')
    ax.set_ylabel('Confidence Score Per Residue')
    ax.set(xlim=(0, i), xticks=np.arange(0,i+scale,scale),  ylim=(0, 100), yticks=np.arange(0,110,10))

    ax.grid(True)

    
    filename=os.path.split(pdb_path)[-1][:-4]+"_conf"
    filepath=os.path.join(out_path,filename)

    plt.savefig(filepath,dpi=600)
    print(f"{filename} generated")

    return 

#---------------------------------------------------------------------------------------------------
# FUNCTIONS FOR OBTAINING BINDING SITE // FOR PDB FILES CONTAINING HETATM 

# gets all protein residues within 4.5A of given ligand
def get_dist(ligand, path, out_path=None):
    if out_path==None:
        out_path = os.getcwd()

    ligand_data=[]
    cls_res=[]

    for array in get_lines(path, "HETATM"):
        mol = array[17:20].strip() 
        if ligand == mol:
            ligand_data.append(array)

    for i in ligand_data: 
        xyz_i= [float(i[30:37].strip()), float(i[38:45].strip()), float(i[46:53].strip())] 
        atom_lines=get_lines(path,"ATOM")
        for j in atom_lines:  
            xyz_j= [float(j[30:37].strip()), float(j[38:45].strip()), float(j[46:53].strip())] 
            pt1= np.array(xyz_i)
            pt2= np.array(xyz_j)
            temp = pt1 - pt2
            sum_sq = np.dot(temp.T, temp) 

            if np.sqrt(sum_sq) < 4.5:
                cls_res.extend([k for k in atom_lines if k[21].strip()==j[21].strip() and k[23:26].strip() == j[23:26].strip()])                             

    cls_res=list(set(cls_res)) 
    cls_res.sort() 
    
    filename = os.path.split(path)[-1][:-4] + "_" + ligand +'.pdb'
    output = os.path.join(out_path, filename)

    with open(output, 'w') as f:

        for element in cls_res:
            f.write(element)
    
    return

#---------------------------------------------------------------------------------------------------
# TERMINAL PARSING

if __name__=="__main__":

    # for fasta  
    if len(sys.argv) in (4,5) and sys.argv[1] == '--fasta':

        in_path=sys.argv[3]
        if len(sys.argv)==5:
            out_path = sys.argv[4]
            if not os.path.isdir(out_path):
                print("Invalid Output Path. Exiting Code")
        else: 
            out_path=None
        
        if sys.argv[2]== '-file':
            if not (os.path.isfile(in_path) and in_path.endswith('.pdb')):
                print("Invalid File Path. Exiting Code")
                sys.exit(1)
            get_solo_fasta(in_path,out_path)
        
        elif sys.argv[2]== '-all':
            if not os.path.isdir(in_path):
                print("Invalid Directory Path. Exiting Code")
                sys.exit(1)
            get_all_fasta(in_path,out_path)
        
        elif sys.argv[2]== '-frag':
            if not os.path.isdir(in_path):
                print("Invalid Directory Path. Exiting Code")
                sys.exit(1)
            get_combined_fasta(in_path,out_path)
    
    # for metadata and confidence score  
    elif len(sys.argv) in (3,4) and sys.argv[1] in ('--metadata','--conf'):
        
        in_path=sys.argv[2]

        if len(sys.argv)==3:
            out_path = sys.argv[4]
            if not os.path.isdir(out_path):
                print("Invalid Output Path. Exiting Code")
                sys.exit(1)
        else: 
            out_path=None


        if sys.argv[1] == '--metadata':
            if not os.path.isdir(in_path):
                print("Invalid Directory Path. Exiting Code")
                sys.exit(1)
            get_metadata(in_path,out_path)
        
        elif sys.argv[1] == '--conf':
            if not os.path.isfile(in_path):
                print("Invalid File Path. Exiting Code")
                sys.exit(1)
            get_conf_graph(in_path,out_path)
        
    # for getting binding site
    elif  len(sys.argv) in (4,5) and sys.argv[1]=='--bsite':
        ligand=sys.argv[2]
        in_path=sys.argv[3]

        if not os.path.isfile(in_path):
            print("Invalid File Path. Exiting Code")
            sys.exit(1)
        elif len(sys.argv)==4:
            out_path = sys.argv[4]
            if not os.path.isdir(out_path):
                print("Invalid Output Path. Exiting Code")
                sys.exit(1)
        else: 
            out_path=None

        get_dist(ligand, in_path, out_path=None)
    
    else:
        print("PDB_modules.py → USAGE:")
        print("To generate Fasta from PDB files:")
        print("\tpython PDB_modules.py --fasta -file <path to pdb file> <path to output dir/ optional>")
        print("\tpython PDB_modules.py --fasta -all <path to dir with pdb files> <path to output dir/ optional>")
        print("\tpython PDB_modules.py --fasta -frag <path to dir with fragment files> <path to output dir/ optional>")
        print("To obtain information of all files in the folder:")
        print("\tpython PDB_modules.py --metadata <path to dir with pdb files> <path to output dir/ optional>")
        print("To generate Residue v/s Confidence graph for a .pdb file:")
        print("\tpython PDB_modules.py --conf <path to pdb file> <path to output dir/ optional>")
        print("To generate an output PDB with only binding site region for a ligand:")
        print("\tpython PDB_modules.py --bsite <HETATM code> <path to pdb file> <path to output dir/ optional>")
