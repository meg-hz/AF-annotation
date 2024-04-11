# AF Annotation Code Set
A collection of scripts used to generate and obtain information from Human protein files retrieved from the Alphafold repository 

## Table of Contents
-  Basic Protein related
	1. PDB_modules.py
	2. UPGet.py
- Superfamily-related
	1. ScopGet.py
- Pocket-related
	1. PockGet.py 
	2. PockComp.py
- FLAPP-related
	1. FAuto.py
	2. FOut.py

---
## UPGet.py
This code searches through the Uniprot database for a protein given a particular input query (or a list of queries in the form of a .txt file) and generates an output file containing the query, protein accession number, protein name and gene ID
### Usage
**1. For a single query**
```
python UPget.py <query>   
```
or
```
python UPGet.py <refseqID> <output dir>
```

**2. For a list of queries**
```
python UPget.py <path to .txt file with queries>
```
or
```
python UPget.py <path to .txt file with queries> <output dir>
```

**Points To Note:**
- This code was initially written to obtain protein details for a given RefSeq ID, hence only single word query input is accepted
---
## ScopGet.py
This code generates the input fasta file for `Superfamily` webtool and analyzes the machine readable output file
### Usage
**1. To get Fasta file input for all files in a folder**
```
python ScopGet.py --all <input dir>
```
or
```
python ScopGet.py --all <input dir> <output dir>
```

**2. To get Fasta file input for a list of files in a folder**
```
python ScopGet.py --list <input dir> <text file containing list>
```
or
```
python ScopGet.py --list <input dir> <text file containing list> <output dir>
```

**3. To analyze output from Superfamily Web Tool**
```
python ScopGet.py --spf <input file>
```
or
```
python ScopGet.py --spf <input file> <output dir>
```

**Points To Note:**
- The Superfamily Web Tool is accessible [here](https://supfam.mrc-lmb.cam.ac.uk/SUPERFAMILY/hmm.html)
---
## PockGet.py
- This code runs `Fpocket`, `Sitehound` and `PocketDepth` for a single protein or a set of proteins when given a filename or a list of file names.
### Usage
**1. For a single file:**
```
python PockGet.py --list <path to input dir> <protein accession/filename>
```

**2. For a set of files:**
```
python PockGet.py --list <path to input directory> <path to .txt file containing file names>
```

**3. For all files in the folder**
```
python PockGet.py --all <path to input directory> 
```

The output pocket files from each pipeline are added to a new/existing directory- `protein_results` directory present in the working directory.

**Points To Note:**
- `Fpocket`, `Sitehound` and `PocketDepth` executable directories and the python file must be in the same parent folder
- `PocketDepth` is a 32bit executable and requires installation of appropriate architectures
- `Sitehound` runs on `python2.7` and will generate errors on higher versions of the compiler
  
- Additionally, for a batch of several proteins of large file size, `Sitehound` and `PocketDepth` take a considerable amount of time- hence it is more efficient to run the them separately rather than using this code.
---
## PockComp.py
- This code compares the outputs between pocket files generated by all three pipelines. It can  generate a TSV with file matching information and consensus pockets as well.
### Usage
**1. To get TSV of matched pockets:**
```
python PockComp.py --txt <path to dir with pocket files>
```
This command generates a .tsv file with paths of matching pockets from all three pipelines<br> 
**2. To get .PDB files of consensus pockets:** 
```
python PockComp.py --pdb <path to protein-files-dir> <path to pocket-files-dir>
```
This command creates a folder called `matches` for each protein folder and generates pocket files containing residues common between the output pocket files from all three pipelines

**3. To move out pockets with low confidence residues**:
```
python PockComp.py --lcf <path to pocket-files-dir>
```
This creates a folder called `low_conf` within the `matches` subfolder for each protein. All consensus pocket files that contain more than 30% residues with lower than 70% confidence score are moved into `low_conf`

**Points To Note:**
- This script is written as the next step to `PockGet.py`, hence the file-handling is also considered accordingly- `pocket-files-dir` refers to the path of `protein_results` directory
- The `pocket-files-dir` refers to path of the directory containing the original .PDB file of the proteins
- This code generates protein pockets for all protein folders in `protein_results` directory and the output  files for each protein are stored within the folder itself
---
## FAuto.py
- This code generates the pairs.txt file and to automate the FLAPP command for a set of files.
### Usage

**1. To run FLAPP between pocket files and template files**
```
FAuto.py --f2temp <path to pocket-file-dir> <path to template directory> <no of cores>
```
This generates an .txt alignment file in the `FLAPP` directory with information on alignment between the pocket files and template files.

**2.. To run FLAPP only among the pocket files**
```
FAuto.py --f2self <path to pocket-file-dir> <no of cores>
```
This generates an .txt alignment file in the `FLAPP` directory with information on alignment between all pocket files

**Points To Note:**
- The environment file within the `FLAPP` directory needs to be setup first. (No need to activate the environment before running the code)
- This script is written as the next step to `PockGet.py`, hence the file-handling is also considered accordingly- `pocket-files-dir`  refers to the path of `protein_results` directory
---
## FOut.py
- This code generates threshold files from the given alignment file by only extracting alignments where the Fmin is within a given cutoff value
- From the threshold output files generated above, it can also provide information on the number of instances of the aligned proteins as well as their pockets

**1. To generate threshold file from the outfile for a given Fmin threshold**
```
FAuto.py --fmin <path to alignment file> <cutoff>
```
This generates a .txt file in the `FLAPP` directory with entries only for alignments that are within the inputted cutoff.

**2. To generate threshold files for intervals of 40,50,60,70,80 and 90**
```
FAuto.py --fmin <path to alignment file> -all
```
This generates a folder containing .txt threshold files in the `FLAPP` directory for a given alignment file for threshold intervals 40 to 90 at intervals of 10 each

**3. To generate reports on protein and pocket instances for a given threshold file**
```
python FAuto.py --fset <path to threshold file>
```
This generates a .tsv file in the same folder as threshold file with information on the number of alignments for each pocket file for a given accession number

**4. To generate reports on protein and pocket instances for all threshold files in the folder**
```
python FAuto.py --fset <path to folder containing threshold files> -all
```
This generates a collection of .tsv files for each threshold file in the folder with information on the number of alignments for each pocket file for a given accession number

**Points To Note:**
- This script is written as the next step to `FAuto.py`, hence the file-handling is also considered accordingly- the `alignment file` referred to here is the output file generated in `FAuto.py`
---
