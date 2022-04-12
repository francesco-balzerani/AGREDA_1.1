# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 08:45:55 2020

@author: fbalzerani
"""

from rdkit import Chem
import pandas as pd
import requests
import urllib
from bs4 import BeautifulSoup
import re

# We specify the root URL.
pubchem_server = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'

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

phenol_explorer = pd.read_excel('./input/metabolites_phenol_explorer.xls')

final_df = pd.DataFrame(index = phenol_explorer.index, columns = ['Name','InChI','ID','Smiles'])

not_done = []
for i in phenol_explorer.index:
    print(str(i)+'/'+str(len(phenol_explorer)-1))
    try:
        id = phenol_explorer.Name[i]
        pubchem = do_request(pubchem_server, 'name',
                             (id + '/property/CanonicalSMILES/JSON'))
        final_df.Smiles[i] = pubchem['PropertyTable']['Properties'][0]['CanonicalSMILES']
        try:
            mol = Chem.MolFromSmiles(final_df.Smiles[i])
            inchi = Chem.MolToInchi(mol)
            final_df.loc[i,'Name'] = phenol_explorer.Name[i]
            final_df.loc[i,'InChI'] = inchi
            final_df.loc[i,'ID'] = phenol_explorer.ID[i]
        except:
            final_df.loc[i,'Name'] = phenol_explorer.Name[i]
            not_done.append(i)
    except:
        try:
            id = int(phenol_explorer.loc[i,'PubChem ID'])
            pubchem = do_request(pubchem_server, 'CID',
                             (str(id) + '/property/CanonicalSMILES/JSON'))
            final_df.Smiles[i] = pubchem['PropertyTable']['Properties'][0]['CanonicalSMILES']
            try:
                mol = Chem.MolFromSmiles(final_df.Smiles[i])
                inchi = Chem.MolToInchi(mol)
                final_df.loc[i,'Name'] = phenol_explorer.Name[i]
                final_df.loc[i,'InChI'] = inchi
                final_df.loc[i,'ID'] = phenol_explorer.ID[i]
            except:
                final_df.loc[i,'Name'] = phenol_explorer.Name[i]
                not_done.append(i)
        except:
            try:
                id = phenol_explorer.ID[i]
                if id < 100:
                    url = "http://moldb.wishartlab.com/molecules/PE0000"+str(id)+".sdf"
                    html = urllib.request.urlopen(url)
                    soup = BeautifulSoup(html)
                    text = soup.get_text()
                    m = re.search('InChI.*', text)
                    inchi = m.group(0)
                    final_df.loc[i,'Name'] = phenol_explorer.Name[i]
                    final_df.loc[i,'InChI'] = inchi
                    final_df.loc[i,'ID'] = phenol_explorer.ID[i]
                elif id < 1000:
                    url = "http://moldb.wishartlab.com/molecules/PE000"+str(id)+".sdf"
                    html = urllib.request.urlopen(url)
                    soup = BeautifulSoup(html)
                    text = soup.get_text()
                    m = re.search('InChI.*', text)
                    inchi = m.group(0)
                    final_df.loc[i,'Name'] = phenol_explorer.Name[i]
                    final_df.loc[i,'InChI'] = inchi
                    final_df.loc[i,'ID'] = phenol_explorer.ID[i]
                else:
                    url = "http://moldb.wishartlab.com/molecules/PE00"+str(id)+".sdf"
                    html = urllib.request.urlopen(url)
                    soup = BeautifulSoup(html)
                    text = soup.get_text()
                    m = re.search('InChI.*', text)
                    inchi = m.group(0)
                    final_df.loc[i,'Name'] = phenol_explorer.Name[i]
                    final_df.loc[i,'InChI'] = inchi
                    final_df.loc[i,'ID'] = phenol_explorer.ID[i]
            except:
                final_df.loc[i,'Name'] = phenol_explorer.Name[i]
                not_done.append(i)
                
new_source = final_df.drop(not_done)
new_source.index = range(len(new_source))
new_source.drop(['ID','Smiles'],axis=1,inplace=True)
# new_source.to_csv('./output/general_source_phen_ex.csv', sep='\t', index=False)

name_inchikey_phen_exp = pd.DataFrame(data = {'name':new_source.Name, 
                                              'inchikey':[Chem.InchiToInchiKey(x) for x in new_source.InChI]})
# new_source.to_csv('./output/name_inchikey_phen_exp.csv', index=False)


