# ScopGet.py
# single hashes for info. double hashes for debugging/alternate code
# generates a tsv file with file name and corresponding scop ids

from bs4 import BeautifulSoup
import os, sys

root=os.getcwd()

def get_scop_tsv(in_dir):

    # check if output file already exists
    out_dir= os.path.join(root,'/protein_results') 
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    output_file_path=os.path.join(out_dir,'ScopList.tsv')

    if os.path.exists(output_file_path):
        response= input("ScopList.tsv already exists. Do you want to rewrite this file?(Y/N) ").strip()
        if  response in ['N','n']:
            print("Rewriting Cancelled.")
            return
        
        elif response not in ['Y', 'y']:
            print("Invalid response. Exiting Code")
            sys.exit(1)

    print("Obtaining Info from model.tab...")
    # generate spf number IDs ⇆ scop ids reference dictionary
    with open(os.path.join(root,'/superfamily/model.tab'),'r') as main:
        model=main.readlines()

    ref={int(i.split()[1]): i.split()[2] for i in model}

    print("Writing OutFile")
    with (output_file_path, 'w') as out_file:

        for file_name in os.listdir(in_dir):
            HTML_path = os.path.join(in_dir, file_name, file_name+'.html')
        
            if os.path.exits(HTML_path):

                with open(os.path.join(in_dir,file_name), 'r') as file:
                    html_read=file.read()

                # finding all the superfamily hyperlinks
                a_tag = BeautifulSoup(html_read, 'html.parser').find_all('a')
                links=[i.get('href') for i in a_tag if 'sunid' in i.get('href')][0:-1:3]

            # finding all SPF id numbers in the file
            id_no=[]
            for i in links:
                scop_id_no= i.split("=")[-1].strip()
                ##print(scop_id_no)

                if scop_id_no !='' and int(scop_id_no) not in id_no:
                    id_no.append(int(scop_id_no))

            # finding all corresponding scop ids from 
            id_ac=[]    
            for element in id_no:
                if element in ref:
                    id_ac.append(ref[element]) 

            # write the output into .tsv file
            out_file.write(f"{file_name[:5]}\t{id_ac}")

if __name__== "__main__":
    if len(sys.argv)==2:
        result_dir=sys.argv[1]
        get_scop_tsv(result_dir)
    else:
        print("ScopGet.py → USAGE:")
        print("\tpython ScopGet.py <path to html file dir>")
