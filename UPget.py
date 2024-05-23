# UPGet.py
# to run every refseq id through uniprot to obtain protein information
# single hashes for info. double hashes for debugging/alternate code

import requests, json, sys, os

def uniprot(ncbi):

    geneID=''
    acc=''
    name=''

    response = requests.get(f'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&query=%28{ncbi}%29')
    # gives a very complicated dictionary but in the form of a json object

    data = json.loads(response.text)['results'] # convert object to dictionary

    info=[]
    for dict in data:
        if dict['entryType']=='UniProtKB reviewed (Swiss-Prot)': # only want the swissprot proteins

            acc=dict['primaryAccession']
            name= dict['proteinDescription']['recommendedName']['fullName']['value']
            organism=dict['organism']['scientificName']
        
            for i in dict['genes']:
                geneID=(i['geneName']['value']+' ')

        info.append(f'{acc}\t{name}\t{organism}\t{geneID}')
    return info

def write_uniprot(refseq_list,out_path):

    count=0
    filename='UP_List.txt'

    while os.path.exists(os.path.join(out_path,filename)):
        count=+1
        filename='UP_List'+str(count)+'.txt'
    
    file_path=os.path.join(out_path,filename) 

    with open(file_path,'w') as out_file:
        out_file.writelines(f"Query\tAccession Number\tProtein Name\tGene ID\n")
        for query in refseq_list:
            for result in uniprot(query):
                out_file.writelines(f"{query}\t{result}\n")
    
    print(f"{filename} generated.")
    
    return


if __name__=="__main__":
    if len(sys.argv)in (3,4) and sys.argv[1] in ('--single', '--list'):
        if sys.argv[1]=='--single':
            if os.path.isfile(sys.argv[2]):
                print("Filepath has been entered in place of query. Retry using '--list' command")
                exit()
            querylist = [sys.argv[2]]

        elif sys.argv[1]== '--list':
            list_path=sys.argv[2]
            #print(out_path)
            if not os.path.isfile(list_path):
                print("Invalid output path. Exiting Code.")
                exit()
            
            with open(list_path,'r') as in_file:
                reference=in_file.readlines()
            querylist=[line.split()[0] for line in reference]

        if len(sys.argv)==4:
            out_path = sys.argv[3]
            if not os.path.isdir(out_path):
                print("Invalid output path. Exiting Code.")
                exit()
            else:
                write_uniprot(querylist,out_path)

        else:
            write_uniprot(querylist,os.getcwd())
   

    else:
        print("UPGet.py → USAGE:")
        print("For Single Query:")
        print("\tpython UPGet.py --single <refseqID>")
        print("\tpython UPGet.py --single <refseqID> <output dir>")
        print("For a Text File with a List of Queries:")
        print("\tpython UPGet.py --list <path to .txt file with list of refseq IDs>")
        print("\tpython UPGet.py --list <path to .txt file with list of refseq IDs> <output dir>")
