# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:37:17 2019

@author: fbalzerani

Create the final sink file and the reduced_reaction file to calculate the rules
"""

# First we load libraries and modules.
import requests
import pandas as pd
import numpy as np
import scipy.io
from rdkit import Chem
import sqlite3
# import urllib
# from bs4 import BeautifulSoup
import re
from bioservices import KEGG

# =============================================================================
# The path should be:
# pwd = '.\AGREDA1_1\01_FilesGenerationAndRetroPathRLApplication'
# =============================================================================

#Connession to the DB
conn = sqlite3.connect('./input/retrorules_dump/mvc.db')
c = conn.cursor()

#Extract information about the metabolites. Useful to compare the cpd_ID if the metabolites in our list not have a related SMILES

c.execute("SELECT * FROM chemical_species")

metabolites_db = c.fetchall()

mets_df = pd.DataFrame(index = range(len(metabolites_db)), columns = ['ID', 'Name', 'cpd_ID', 'HMDB_ID', 'KEGG_ID'])
mets_df.loc[:,'ID'] = [x[0] for x in metabolites_db]
mets_df.loc[:,'Name'] = [x[1] for x in metabolites_db]
mets_df.loc[:,'cpd_ID'] = [x[12] for x in metabolites_db]
mets_df.loc[:,'KEGG_ID'] = [x[7] for x in metabolites_db]
mets_df.loc[:,'HMDB_ID'] = [x[6] for x in metabolites_db]

##### THIS PART DEVOTED TO DATA MINING FROM PUBCHEM #####
# We specify the root URL.
pubchem_server = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'
kegg_server = KEGG(verbose = False)
# kegg_server = "https://www.kegg.jp/dbget-bin/www_bget?-f+m+compound+"
hmdb_server = 'http://www.hmdb.ca/metabolites/'

# We create a function that is going to perform the web request.
def do_request(server, service, *args, **kwargs):
    url_params = ''
    for a in args:
        if a is not None:
            url_params += '/' + a
    req = requests.get('%s/%s%s' % (server, service, url_params),
                       params=kwargs,
                       headers={'Content-Type':'application/json'},
                       timeout = 5)
    if not req.ok:
        req.raise_for_status()
    return req.json()

###### IMPORTANT WARNING START #####
#
# You can use this work-around if you want to, but need to save the Matlab file
# in the proper format. Otherwise, this WON'T WORK. Save it in Matlab with the
# following command:
# 
# save('test.mat', '-v7')
# 
# Obviously, give the file whatever name you want, but keep the -v7 flag.
# Otherwise, the file gets saved as a v7.3 file, and we can't easily open that.
#
##### IMPORTANT WARNING ENDED #####

##### THE ACTUAL LOOP STARTS HERE #####
# We load the data.
mol_library = pd.read_excel('./output/library_mets.xlsx')
rxns_library = pd.read_excel('./output/library_rxns.xlsx')
mat = scipy.io.loadmat('./output/stoichiometric_matrix.mat')
S = mat['S']

sink = pd.DataFrame(columns = ['name','inchi'], index = mol_library.index)

# We now search for each of the molecules in the list.
not_done_1 = []
pubchem_count = 0
hmdb_count = 0
kegg_count = 0

for i in mol_library.index:

    print(str(i)+'/'+str(len(mol_library)-1))
    
    # Use the name to go through pubchem
    try:
        id = mol_library.Names[i]
        pubchem = do_request(pubchem_server, 'name',
                             (id + '/property/InChI/JSON'))
        inchi = pubchem['PropertyTable']['Properties'][0]['InChI']
        sink.loc[i,'name'] = mol_library.Names[i]
        sink.loc[i,'inchi'] = inchi
        pubchem_count += 1
    except:
        # extract the corresponding name from the DB and go through pubchem
        try:
            tmp_trial = mets_df[mets_df['cpd_ID']==mol_library.ID[i]]
            id = tmp_trial.Name
            id = id.reset_index(drop = True)
            id = id[0]
            pubchem = do_request(pubchem_server, 'name',
                             (id + '/property/InChI/JSON'))
            inchi = pubchem['PropertyTable']['Properties'][0]['InChI']
            sink.loc[i,'name'] = mol_library.Names[i]
            sink.loc[i,'inchi'] = inchi
            pubchem_count += 1
        except:
            # extract the corresponding HMDB ID and go through HMDB
            try:
                tmp_trial = mets_df[mets_df['cpd_ID']==mol_library.ID[i]]
                id = tmp_trial.HMDB_ID
                id = id.reset_index(drop = True)
                id = id[0]
                response = requests.request('GET',hmdb_server+id+'.xml')
                tmp_response = response.text
                # check if we can reach directly the inchi
                inchi = re.findall('InChI Identifier</th><td><div class="wrap">(\S+)</div>',tmp_response)[0]
                sink.loc[i,'name'] = mol_library.Names[i]
                sink.loc[i,'inchi'] = inchi     
                hmdb_count += 1
            except:
                try:
                    tmp_trial = mets_df[mets_df['cpd_ID']==mol_library.ID[i]]
                    id = tmp_trial.HMDB_ID
                    id = id.reset_index(drop = True)
                    id = id[0]
                    response = requests.request('GET',hmdb_server+id+'.xml')
                    tmp_response = response.text
                    # check if we can reach directly the inchi
                    inchi = re.findall('InChI Identifier</th><td><div class="word-break-all">(\S+)</div>',tmp_response)[0]
                    sink.loc[i,'name'] = mol_library.Names[i]
                    sink.loc[i,'inchi'] = inchi
                    hmdb_count += 1
                except:
                    # extract the corresponding KEGG ID and go through KEGG
                    try:
                        tmp_kegg = mets_df[mets_df['cpd_ID']==mol_library.ID[i]]
                        id = tmp_kegg.KEGG_ID
                        id = id.reset_index(drop = True)
                        id=id[0]
                        # url = kegg_server+id
                        # html = urllib.request.urlopen(url)
                        # soup = BeautifulSoup(html)
                        # text = soup.get_text()
                        text = kegg_server.get(id,option = 'mol')
                        mol = Chem.MolFromMolBlock(text)
                        inchi = Chem.MolToInchi(mol)
                        sink.loc[i,'name'] = mol_library.Names[i]
                        sink.loc[i,'inchi'] = inchi
                        kegg_count += 1
                    except:
                        sink.loc[i,'name'] = mol_library.Names[i]
                        not_done_1.append(i)
            
tmp = list(np.where(sink['inchi'].isin([''])))
tmp = np.array(tmp).tolist()

j=0
for i in tmp[0]:
    print(str(j)+'/'+str(len(tmp[0])-1))
    j+=1
    try:
        tmp_trial = mets_df[mets_df['cpd_ID']==mol_library.ID[i]]
        id = tmp_trial.HMDB_ID
        id = id.reset_index(drop = True)
        id = id[0]
        response = requests.request('GET',hmdb_server+id+'.xml')
        tmp_response = response.text
        # check if we can reach directly the inchi
        inchi = re.findall('InChI Identifier</th><td><div class="wrap">(\S+)</div>',tmp_response)[0]
        sink.loc[i,'name'] = mol_library.Names[i]
        sink.loc[i,'inchi'] = inchi
    except:
        try:
            tmp_trial = mets_df[mets_df['cpd_ID']==mol_library.ID[i]]
            id = tmp_trial.HMDB_ID
            id = id.reset_index(drop = True)
            id = id[0]
            response = requests.request('GET',hmdb_server+id+'.xml')
            tmp_response = response.text
            # check if we can reach directly the inchi
            inchi = re.findall('InChI Identifier</th><td><div class="word-break-all">(\S+)</div>',tmp_response)[0]
            sink.loc[i,'name'] = mol_library.Names[i]
            sink.loc[i,'inchi'] = inchi
        except:
            # extract the corresponding KEGG ID and go through KEGG
            try:
                tmp_kegg = mets_df[mets_df['cpd_ID']==mol_library.ID[i]]
                id = tmp_kegg.KEGG_ID
                id = id.reset_index(drop = True)
                id=id[0]
                # url = kegg_server+id
                # html = urllib.request.urlopen(url)
                # soup = BeautifulSoup(html)
                # text = soup.get_text()
                text = kegg_server.get(id,option = 'mol')
                mol = Chem.MolFromMolBlock(text)
                inchi = Chem.MolToInchi(mol)
                sink.loc[i,'name'] = mol_library.Names[i]
                sink.loc[i,'inchi'] = inchi
            except:
                sink.loc[i,'name'] = mol_library.Names[i]
                not_done_1.append(i)

tmp = list(np.where(sink['inchi'].isin([''])))
tmp = np.array(tmp).tolist()

j=0
for i in tmp[0]:
    print(str(j)+'/'+str(len(tmp[0])-1))
    j+=1
    try:
        id = mol_library.Names[i]
        pubchem = do_request(pubchem_server, 'name',
                             (id + '/property/InChI/JSON'))
        inchi = pubchem['PropertyTable']['Properties'][0]['InChI']
        sink.loc[i,'name'] = mol_library.Names[i]
        sink.loc[i,'inchi'] = inchi
    except:
        sink.loc[i,'name'] = mol_library.Names[i]
        not_done_1.append(i)

tmp = list(np.where(sink['inchi'].isin([''])))
tmp = np.array(tmp).tolist()
not_done_1 = not_done_1 + tmp[0]
not_done_1 = np.unique(not_done_1)
not_done_1 = np.array(not_done_1).tolist()
not_done_1.sort()

# We delete the reactions that contain the metabolites that had no SMARTS first,
# then we delete the rows of those metabolites.
to_del = np.unique(np.nonzero(S[not_done_1,])[1])
Snew = S.toarray()
Snew = np.delete(Snew, to_del, axis = 1)
Snew = np.delete(Snew, not_done_1, axis = 0)
# np.save('./output/S_reduced.npy',Snew)

rxns_library2 = rxns_library.drop(to_del)
rxns_library2.index = range(0,Snew.shape[1])
# rxns_library2.to_csv('./output/reduced_reactions_bounds.csv', sep = ',', index = False)

sink2 = sink.drop(not_done_1)

# Removing cofactor compounds, previously written manually

cofactor = pd.read_csv('./input/cofactor.csv')
cofactors_inchi =  sink2.loc[sink2['name'].isin(cofactor.name),]
# cofactors_inchi.to_csv('./output/cofactors_inchi.csv')

sink3 = sink2.loc[~sink2['name'].isin(pd.unique(cofactor['name'].values)),]

sink3.index = range(len(sink3))

# Load the mixture file, previously written manually
mixture = pd.read_csv('./input/mixture.csv')

sink4 = sink3.loc[~sink3['name'].isin(pd.unique(mixture['name'].values)),]
sink4.drop_duplicates(inplace=True)

# Modify manually some metabolites: those that using the name are related to chloride
# form or "mixture" form that we know are generated wrongly from the network annotation
# and then considering those compound with ions, removing the ions from the inchi
# maintaining the compounds in the sink. The ions should not change the chemical similarity
sink4.index = range(len(sink4))
sink4.loc[sink4['name'].isin(['IdB 1027']),'inchi'] = 'InChI=1S/C15H10O6/c16-8-4-11(18)9-6-13(20)15(21-14(9)5-8)7-1-2-10(17)12(19)3-7/h1-6H,(H4-,16,17,18,19,20)/p+1'
sink4.loc[sink4['name'].isin(['Ephdine']),'inchi'] = 'InChI=1S/C15H10O7/c16-7-3-9(17)8-5-12(20)15(22-13(8)4-7)6-1-10(18)14(21)11(19)2-6/h1-5H,(H5-,16,17,18,19,20,21)/p+1' 
sink4.loc[sink4['name'].isin(['Mirtillin']),'inchi'] = 'InChI=1S/C21H20O12/c22-6-15-17(28)18(29)19(30)21(33-15)32-14-5-9-10(24)3-8(23)4-13(9)31-20(14)7-1-11(25)16(27)12(26)2-7/h1-5,15,17-19,21-22,28-30H,6H2,(H4-,23,24,25,26,27)/p+1/t15-,17-,18+,19-,21-/m1/s1'
sink4.loc[sink4['name'].isin(['pelargonidin-3,5-diglucoside']),'inchi'] = 'InChI=1S/C27H30O15/c28-8-17-19(32)21(34)23(36)26(41-17)39-15-6-12(31)5-14-13(15)7-16(25(38-14)10-1-3-11(30)4-2-10)40-27-24(37)22(35)20(33)18(9-29)42-27/h1-7,17-24,26-29,32-37H,8-9H2,(H-,30,31)/p+1/t17-,18-,19-,20-,21+,22+,23-,24-,26-,27-/m1/s1'

for x in sink4.index:
    sink4.loc[x,'inchi'] = sink4.loc[x,'inchi'].replace('.Co','')
    sink4.loc[x,'inchi'] = sink4.loc[x,'inchi'].replace('.Fe','')
    sink4.loc[x,'inchi'] = sink4.loc[x,'inchi'].replace('.Na','')
    sink4.loc[x,'inchi'] = sink4.loc[x,'inchi'].replace('.Mg','')
    sink4.loc[x,'inchi'] = sink4.loc[x,'inchi'].replace('.Ni','')
    sink4.loc[x,'inchi'] = sink4.loc[x,'inchi'].replace('.2Na','')
    sink4.loc[x,'inchi'] = sink4.loc[x,'inchi'].replace('.4Na','')
    sink4.loc[x,'inchi'] = sink4.loc[x,'inchi'].replace('.Mo','')
    sink4.loc[x,'inchi'] = sink4.loc[x,'inchi'].replace(';','')
    sink4.loc[x,'inchi'] = sink4.loc[x,'inchi'].replace('.','')
    
#sink4.to_csv('./output/general_sink.csv', sep=',', index=False)

sink5 = sink4.copy()
sink5.loc[:,'InChIKey'] = [Chem.InchiToInchiKey(x) for x in sink5.inchi]
#sink5.to_excel("./output/name_inchikey.xlsx", index = False)
