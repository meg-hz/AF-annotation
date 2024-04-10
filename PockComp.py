# PockComp.py
# single hashes for info. double hashes for debugging/alternate code
# code to obtain text fike about matching pockets bw sitehound, fpocket and pocketdepth 
# as well as pocket files with residues matching among all three

import os,sys, shutil

# returns a dictionary of unique residues from a .pdb file
def pocket_dict(filepath):

    pockets={}

    for PDB_file in os.listdir(filepath):
        residues=[]
        if '.pdb' in PDB_file:
            with open(filepath+'/'+PDB_file, 'r') as file:
                PDB_lines=file.readlines()

                for line in PDB_lines:
                    if line.startswith('ATOM') and line[17:26].strip() not in residues:
                        residues.append(line[17:26].strip())
                    
            pockets[PDB_file]=residues

    ##for i in pockets:     
        ##print(f'{i}:{(pockets[i])}')
    return pockets
    
# returns a dictionary with matching files + array of matching residues if main path to pocket folders is given
def pocket_matches(pock_dir):
    all_matches={}

    print('Extracting protein information...')

    for protein_folder in os.listdir(pock_dir):

            ##filename=str(protein_folder.split('-')[1]+'-'+protein_folder.split('-')[2])

            # path to folders containing pocket files
            accession_filepath=pock_dir+'/'+protein_folder
            SH_filepath=accession_filepath+'/'+ protein_folder+'_sth/pockets'
            FP_filepath=accession_filepath+'/'+ protein_folder+'_fpk/pockets'
            PD_filepath=accession_filepath+'/'+ protein_folder+'_pkd'

            # dictionary of {filename1: [array of residues], filename2: [array of residues]}
            SH_pocket=pocket_dict(SH_filepath)
            FP_pocket=pocket_dict(FP_filepath)
            PD_pocket=pocket_dict(PD_filepath)

            match_pairs=[]

            for SH_key in SH_pocket:
                ##print (f'{SH_key}:{SH_pocket[SH_key]}')
                ref= SH_pocket[SH_key] # the array of residues for one SH pocket file
                for FP_key in FP_pocket:
                    count=[]
                    for residue in ref:
                        if residue in FP_pocket[FP_key]:
                            count.append(residue)
                            # all common residues between sitehound and pocketdepth

                    for PD_key in PD_pocket:
                        count2=[]
                        for residue in count:
                            if residue in PD_pocket[PD_key]:
                                count2.append(residue)
                                # all common residues between all   
                
                        if len(count2)>=4:
                            match_pairs.append((str(SH_filepath+'/'+SH_key),str(FP_filepath+'/'+FP_key),str(PD_filepath+'/'+PD_key),count2))
                            ##print((str(SH_filepath+'/'+SH_key),str(FP_filepath+'/'+FP_key)))

            all_matches[protein_folder]=match_pairs 
            # ^ dictonary will look something like this-
            # { accssion1: [ (SH_pocket1, FP_pocket1, PD_pocket1, [matched residues]),
            #                (SH_pocket2, FP_pocket2, PD_pocket2, [matched residues]) ]
            #   accssion2: [ (SH_pocket1, FP_pocket1, PD_pocket3, [matched residues]),
            #                (SH_pocket2, FP_pocket2, PD_pocket4, [matched residues]) ] }
    
            ##break

    return all_matches

# generates a TSV of matching pocket files with matching residues
def matches_tsv(pock_dir):
    output_file=os.path.join(pock_dir,'pocket_matches.tsv')

    count=0
    while os.path.exists(output_file):
        count+=1
        output_file=os.path.join(pock_dir,'pocket_matches'+str(count)+'.tsv')

    pocket_dir = pocket_matches(pock_dir)

    with open(output_file,'w') as output:
        for protein, matching in pocket_dir.items():
            for sett in matching:
                output.write(protein)
                for element in sett:
                    output.write('\t'+str(element))
                output.write('\n')
            print(f'Pockets added for {protein}')

    print("POCKET INFO TEXT FILE GENERATED")
    print("------------------------------------------------------")
    print()

# extracts ATOM lines from a .pdb file for a given set of residues
def extract_atom_lines(pdb_file, residues):
    atom_lines = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                residue_name = line[17:27].strip()
                if residue_name in residues:
                    atom_lines.append(line)
    return atom_lines

# generates consensus pocket files for each set of matching pdb files
def matches_pdb(input_dir,output_dir):

    # dictionary of matched residues
    pocket_dir=pocket_matches(output_dir)

    for protein_file in pocket_dir:
        ##print(f'{protein_file}:{items}')

        new=os.path.join(output_dir,protein_file,'all_match')
        if not os.path.exists(new):
            os.mkdir(new)

        i=0
        pdb_file = os.path.join(input_dir, protein_file+'.pdb')
            
        for matched in pocket_dir[protein_file]:
            i+=1
            residues = matched[-1]  # list of matching residues

            # Extract atom lines for each residue
            atom_lines = extract_atom_lines(pdb_file, residues)

            # Write the extracted atom lines to a new pdb file
            pocket_name=protein_file+'-p'+str(i)+'.pdb'
            output_file = os.path.join(new,pocket_name)
            with open(output_file, 'w') as outfile:
                outfile.writelines(atom_lines)
            ##print(f'{pocket_name} generated')
        print(f'{i} files generated for {protein_file}')
        

# get a set of confidence scores for Alphafold 
def check_numeric_values(file_path):
    numeric_values = []
    with open(file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM'):
                try:
                    numeric_value = float(line[60:66].strip())
                    numeric_values.append(numeric_value)
                except ValueError:
                    pass
    unique_numeric_values = set(numeric_values)

    threshold = 70.0
    count_below_threshold = sum(1 for value in unique_numeric_values if value < threshold)
    percentage_below_threshold = (count_below_threshold / len(unique_numeric_values)) * 100

    # returns number of residues with
    return percentage_below_threshold > 30

def move_lowconf_pdb(folder_path):
    for dirs in os.listdir(folder_path):
        print(f'FOR {dirs}')
        pocket_dirs= os.path.join(folder_path, dirs,'all_match')

        for file in os.listdir(pocket_dirs):
            ##print(file)
            if file.endswith('.pdb'):
                pdb_path = os.path.join(pocket_dirs, file)
                if check_numeric_values(pdb_path):

                    low_conf_folder=os.path.join(pocket_dirs,'Low_conf')

                    if not os.path.exists(low_conf_folder):
                        os.mkdir(low_conf_folder)

                    shutil.move(pdb_path,low_conf_folder)
                    print(f"Moved: {file}")
        print()
    print("All Low Confidence files have been moved")

#---------------------------------------------------------------------------------------------------
# MAIN PROGRAM STARTS HERE

if __name__ == "__main__":
    if len(sys.argv) ==3 and sys.argv[1] == '--tsv':
            pockets = sys.argv[2]
            matches_tsv(pockets)

    elif len(sys.argv) ==3 and sys.argv[1] == '--lcf':
            pocket_dir=sys.argv[2]
            move_lowconf_pdb(pocket_dir)

    elif len(sys.argv) == 4 and sys.argv[1] == '--pdb':
        indir = sys.argv[2] # directory with the input PDB
        outdir = sys.argv[3] # directory where the pocket files are
        matches_pdb(indir, outdir)
      
    else:
        print("PockComp.py → USAGE:")
        print ("For TSV of matched pockets: \n\tpython PockComp.py --tsv <path to dir with pocket files>")
        print ("For PDB of consensus pockets: \n\tpython PockComp.py --pdb <path to protein files dir> <path to pocket files dir>")
        print ("For moving pockets with low conf:\n\tpython PockComp.py --lcf <path to dir with pocket files>")
