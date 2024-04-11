# UPGet.py
# single hashes for info. double hashes for debugging/alternate code
# run every refseq id through uniprot to obtain protein information

import requests, json, sys, os

root=os.getcwd()

def uniprot(ncbi):

    geneID=''
    acc=''
    name=''

    response = requests.get(f'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&query=%28{ncbi}%29')
    # gives a very complicated dictionary but in the form of a json object

    data = json.loads(response.text)['results'] # convert object to dictionary

    for dict in data:
        if dict['entryType']=='UniProtKB reviewed (Swiss-Prot)': # only want the swissprot protein

            acc=dict['primaryAccession']
            name= dict['proteinDescription']['recommendedName']['fullName']['value']
        
            for i in dict['genes']:
                geneID+=(i['geneName']['value']+' ')

            break
    
    return f'{acc}\t{name}\t{geneID}'

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
            out_file.writelines(f"{query}\t{uniprot(query)}\n")
    
    return


if __name__=="__main__":
    if len(sys.argv)in (2,3):
        list_path=sys.argv[1]
        
        if len(sys.argv)==3:
            out_path=sys.argv[2]
            print(out_path)
            if not os.path.isdir(out_path):
                print("Invalid output path. Exiting Code.")
        else:
            out_path=root

        if os.path.exists(list_path):

            if not (os.path.isfile(list_path) and list_path.endswith('.txt')):
                print("Error in inputted file path. Exiting code.")
                exit(1)
            else:
                querylist=[]
                with open(list_path,'r') as in_file:
                    reference=in_file.readlines()

                    for line in reference:
                        if line.split()[0] not in querylist:
                            querylist.append(line.split()[0])

                    write_uniprot(querylist,out_path)

        else:
            write_uniprot([list_path],out_path)

        

    else:
        print("UPGet.py → USAGE:")
        print("For Single Query:")
        print("\tpython UPGet.py <refseqID>")
        print("\tpython UPGet.py <refseqID> <output dir>")
        print("For a Text File with a List of Queries:")
        print("\tpython UPGet.py <path to .txt file with list of refseq IDs>")
        print("\tpython UPGet.py <path to .txt file with list of refseq IDs> <output dir>")
