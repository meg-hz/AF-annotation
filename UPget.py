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

def write_uniprot(refseq_list):

    count=0
    filename='UP_List.txt'
    while os.path.exists(os.path.join(root,filename)):
        count=+1
        filename='UP_List'+str(count)+'.txt'

    with open(filename,'w') as out_file:
        for query in refseq_list:
            out_file.writelines(f"{query}\t{uniprot(query)}")
    
    return

if __name__=="__main__":
    if len(sys.argv)==2:
        list_path=sys.argv[1]

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

                    write_uniprot(querylist)

        else:
            write_uniprot([list_path])

    else:
        print("UPGet.py → USAGE:")
        print("\tpython UPGet.py <path to .txt file with list of refseq IDs>")
        print("\tpython UPGet.py <refseqID>")
