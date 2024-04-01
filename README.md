# AF Annotation Code Set
A collection of scripts used to generate and obtain information from Human protein files retrieved from the Alphafold repository
(Note: Incomplete readme file. Details on 2/5 files still pending)

## Table of Contents
1. PDB_modules.py** 
2. Superfamily.py**
3. [PGet.py](#PGet.py)
4. [PComp.py](#PComp.py)

---

## PGet.py
This code runs Fpocket, Sitehound and PocketDepth for a single protein or a set of proteins when given a filename or a list of file names.

### To Consider 
- Fpocket, Sitehound and PocketDepth executable directories and the python file must be in the same parent folder
- PocketDepth is a 32bit executable and requires installation of appropriate architectures
  
- Note that for a batch of several proteins of large file size, Sitehound and PocketDepth take a considerable amount of time- hence it is more efficient to run the them separately rather than using this code.

### Usage
**1. For a single file:**
```
python PGet.py <path to input dir> <protein accession/filename>
```

**2. For a set of files:**
```
python PGet.py <path to input directory> <path to filen list.txt>
```

The script generates a directory called `result_pockets` in the working directory with subfolders for each .PDB file.

---

## PComp.py
This code compares the outputs between pocket files generated by all three pipelines. It can  generate a CSV with file matching information and consensus pockets as well 

### Usage
**1. To get CSV of matched pockets:**
```
python PComp.py --txt <path to dir with pocket files>
```
This command generates a `.tsv` file with paths of matching pockets from all three pipelines

**2. To get .PDB files of consensus pockets:** 
```
python AF_residues.py --pdb <path to protein-files-dir> <path to pocket-files-dir>
```
This command creates a folder called `matches` for each protein folder and generates pocket files containing residues common between the output pocket files from all three pipelines

**Note:**
- this script is written as the next step to [previous code](#PGet.py) and file-handling is also considered accordingly; hence `pocket-files-dir` refers to the path of `result_pockets` 
- `pocket-files-dir` refers to the directory containing the original .PDB file of the proteins
- this code generates protein pockets for all protein folders in `result_pockets` and the output  files for each protein are stored within the folder itself

**3. To move out pockets with low confidence residues**:
```
python AF_residues.py --lcf <path to pocket-files-dir>
```

This creates a folder called `low_conf` within the `matches` subfolder for each protein. All consensus pocket files that contain more than 30% residues with lower than 70% confidence score are moved into `low_conf`

---
