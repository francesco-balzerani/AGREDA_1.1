# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 12:30:37 2020

@author: fbalzerani

Create the rules (with and without bounds. Without for RetroPath RL application
                  with for the assignation of the bounds to the results in the 
                  optimization)
"""

# First we load libraries and modules.
import sqlite3
import pandas as pd
from rdkit import Chem
import numpy as np
import os
import re
import urllib
from bs4 import BeautifulSoup
from bioservices import KEGG
import requests
import pubchempy as pcp
import pickle

# =============================================================================
# The path should be:
# pwd = '.\AGREDA1_1\01_FilesGenerationAndRetroPathRLApplication'
# =============================================================================

# Create the rules folder

if not os.path.isdir("./output/rules"):
    os.mkdir("./output/rules")

#function to convert a tuple to dictionary

def Convert(tup,di):
    di = dict(tup)
    return di

metanetx_server = 'https://www.metanetx.org/chem_info/'

#create kegg html
kegg = 'http://rest.kegg.jp/find/reaction/'
k = KEGG(verbose = False)

#Connession to the DB
conn = sqlite3.connect('./input/retrorules_dump/mvc.db')
c = conn.cursor()

#Create dictionary to assign EC numbers to reactions where they have

c.execute("SELECT * FROM ec_reactions")

tup_ecnumbers_reactions = c.fetchall()

dictionary_ec = {}

dictionary_ec = Convert(tup_ecnumbers_reactions,dictionary_ec)

#Create dictionary to assign MXN id to reactions

c.execute("SELECT * FROM reactions")

reactions = c.fetchall()

dictionary_id_reactions = {}
dictionary_seed_reactions = {}

tup_id_mnx = [(x[0], x[1]) for x in reactions]
tup_id_seed = [(x[0], x[2]) for x in reactions]

dictionary_id_reactions = Convert(tup_id_mnx,dictionary_id_reactions)
dictionary_seed_reactions = Convert(tup_id_seed,dictionary_seed_reactions)

#Create dictionary to assign the smiles to reactions

c.execute("SELECT * FROM smiles")

smiles = c.fetchall()

dictionary_smiles = {}

tup_smiles = [(x[0], x[1]) for x in smiles]

dictionary_smiles = Convert(tup_smiles,dictionary_smiles)

#Create dictionary to assign the smarts to reactions

c.execute("SELECT * FROM smarts")

smarts = c.fetchall()

dictionary_smarts = {}

tup_smarts = [(x[0], x[1]) for x in smarts]

dictionary_smarts = Convert(tup_smarts,dictionary_smarts)

#Obtain the main table with informations about rxnsID, substrateID, diameters, isStereo, smartsID, smilesID, direction, score

c.execute("SELECT * FROM rules")

entire_rules = c.fetchall()

#Create dictionary to assign the mnxm id to the substrate compounds for the Rule ID

c.execute("SELECT * FROM chemical_species")

metabolites_db = c.fetchall()

dictionary_compounds = {}

tup_compounds = [(x[0], x[2]) for x in metabolites_db]

dictionary_compounds = Convert(tup_compounds,dictionary_compounds)

#Create the Dataframe containing all the informations about the reactions in the DB to compare with our reactions of interest

rules = pd.DataFrame(index = range(len(entire_rules)), columns = ['Rule_ID','Reaction_ID','Diameter','Direction','Rule_SMILES','Substrate_ID','EC number','Score','ID','smiles_id','seed id','smarts_id','Rule_SMARTS_EXP'])

rules.loc[:,'ID'] = [x[0] for x in entire_rules]
rules.loc[:,'EC number'] = rules['ID'].apply(lambda x: dictionary_ec.get(x))
rules.loc[:,'Diameter'] = [x[2] for x in entire_rules]
rules.loc[:,'Score'] = [x[7] for x in entire_rules]
rules.loc[:,'Substrate_ID'] = [x[1] for x in entire_rules]
rules.loc[:,'Direction'] = [x[6] for x in entire_rules]

reactions_names = rules['ID'].apply(lambda x: dictionary_id_reactions.get(x))
substrate_names = rules['Substrate_ID'].apply(lambda x: dictionary_compounds.get(x))

rules.loc[:,'Reaction_ID'] = reactions_names
rules.loc[:,'Substrate_ID'] = substrate_names

rules.loc[:,'Rule_ID'] = reactions_names+'_'+substrate_names

rules.loc[:,'smarts_id'] = [x[4] for x in entire_rules]
rules.loc[:,'Rule_SMARTS_EXP'] = rules['smarts_id'].apply(lambda x: dictionary_smarts.get(x))

rules.loc[:,'smiles_id'] = [x[5] for x in entire_rules]
rules.loc[:,'Rule_SMILES'] = rules['smiles_id'].apply(lambda x: dictionary_smiles.get(x))

rules.loc[:,'seed id'] = rules['ID'].apply(lambda x: dictionary_seed_reactions.get(x))

##############################################################################
#Create dictionary to assign the mnxm id to the substrate compounds to remove the 
# cofactor as substrate from the rules

c.execute("SELECT * FROM chemical_species")

metabolites_db = c.fetchall()

dictionary_names = {}

tup_names = [(x[1], x[2]) for x in metabolites_db]

dictionary_names = Convert(tup_names,dictionary_names)

# Removing rules that have cofactor compounds as substrate

cofactor = pd.read_csv('./input/cofactor.csv')

cofactor_id = cofactor['name'].apply(lambda x: dictionary_names.get(x))

cofactor_id.dropna(inplace=True)

manual_id_to_add = pd.Series(['MNXM1','MNXM01','MNXM2','MNXM2006'])

cofactor_id = cofactor_id.append(manual_id_to_add, ignore_index = True)

rules = rules.loc[~rules['Substrate_ID'].isin(pd.unique(cofactor_id.values)),]
rules.index = range(len(rules))
##############################################################################

#Create the dictionaries for the reaction's substrates and reaction's products
#to generate properly the product's SMILES

if os.path.isfile("./output/dictionaries.pkl"):
    # Getting back the objects:
    with open('./output/dictionaries.pkl', 'rb') as f: 
        dictionary_reac_subs, dictionary_reac_prods = pickle.load(f)
else:
    c.execute("SELECT * FROM reaction_products")
    
    reaction_products = c.fetchall()
    
    idx_reaction = [x[0] for x in reaction_products]
    
    dictionary_reac_prods = {}
    count = 0
    for ID in np.unique(idx_reaction):
        print(str(count)+'/'+str(len(np.unique(idx_reaction))))
        count += 1
        dictionary_reac_prods[dictionary_id_reactions[ID]] = [dictionary_compounds[reaction_products[x][1]] for x in range(len(reaction_products)) if reaction_products[x][0] == ID]
    
    c.execute("SELECT * FROM reaction_substrates")
    
    reaction_substrates = c.fetchall()
    
    idx_reaction = [x[0] for x in reaction_substrates]
    
    dictionary_reac_subs = {}
    count = 0
    for ID in np.unique(idx_reaction):
        print(str(count)+'/'+str(len(np.unique(idx_reaction))))
        count += 1
        dictionary_reac_subs[dictionary_id_reactions[ID]] = [dictionary_compounds[reaction_substrates[x][1]] for x in range(len(reaction_substrates)) if reaction_substrates[x][0] == ID]
    
    with open('./output/dictionaries.pkl', 'wb') as f:
        pickle.dump([dictionary_reac_subs,dictionary_reac_prods], f)

#Create a dataframe to have all the information about the compounds, needed for
#smiles' generation of substrates and products

c.execute("SELECT * FROM chemical_species")

metabolites_db = c.fetchall()

mets_df = pd.DataFrame(index = range(len(metabolites_db)), columns = ['ID', 'Name', 'METANETX_ID', 'cpd_ID', 'HMDB_ID', 'KEGG_ID'])
mets_df.loc[:,'ID'] = [x[0] for x in metabolites_db]
mets_df.loc[:,'Name'] = [x[1] for x in metabolites_db]
mets_df.loc[:,'METANETX_ID'] = [x[2] for x in metabolites_db]
mets_df.loc[:,'cpd_ID'] = [x[12] for x in metabolites_db]
mets_df.loc[:,'KEGG_ID'] = [x[7] for x in metabolites_db]
mets_df.loc[:,'HMDB_ID'] = [x[6] for x in metabolites_db]

if os.path.isfile("./output/mets_dictionary.pkl"):
    # Getting back the objects:
    with open('./output/mets_dictionary.pkl', 'rb') as f: 
        [mets_df,mets_dictionary] = pickle.load(f)
else:
    for idx in mets_df.index:
        print(str(idx)+'/'+str(mets_df.index[-1]))
        if mets_df.loc[idx,'METANETX_ID'] == 'MNXM10258':
            mets_df.loc[idx, 'Smiles'] = "CC(CCC=C(C)C)C1CCC2(C1(CCC34C2CCC5C3(C4)CCC(C5(C)CO)O)C)C"
            continue
        elif mets_df.loc[idx,'METANETX_ID'] == 'MNXM57425':
            mets_df.loc[idx, 'Smiles'] = "OP([O-])([O-])=O"
            continue
        elif mets_df.loc[idx,'METANETX_ID'] == 'MNXM861':
            mets_df.loc[idx, 'Smiles'] = 'electron_not_needed_smiles'
            continue
        p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
        pp = p.text
        try:
            mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
            if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                id = mets_df.loc[idx,'Name']
                compound = pcp.get_compounds(id, 'name')
                if len(compound) == 0:
                    try:
                        id = mets_df.loc[idx,'KEGG_ID']
                        mets_df.loc[idx, 'Smiles'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
                    except:
                        mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                        continue
                else:
                    mets_df.loc[idx, 'Smiles'] = compound[0].canonical_smiles
            try:
                sanit = Chem.SanitizeMol(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles']))
                check_inchi = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles'])))
                if check_inchi == None:
                    id = mets_df.loc[idx,'Name']
                    compound = pcp.get_compounds(id, 'name')
                    if len(compound) == 0:
                        try:
                            id = mets_df.loc[idx,'KEGG_ID']
                            mets_df.loc[idx, 'Smiles'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
                        except:
                            # old
                            # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                            # continue
                            # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                            # can be sanitizable and pass through the inchi, thus, if present, put the 
                            # SMILES of metanetx
                            p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                            pp = p.text
                            try:
                                mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                    mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                    continue
                            except:
                                mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                continue
                    else:
                        mets_df.loc[idx, 'Smiles'] = compound[0].canonical_smiles
                    try:
                        sanit = Chem.SanitizeMol(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles']))
                        check_inchi = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles'])))
                        if check_inchi == None:
                            try:
                                id = mets_df.loc[idx,'KEGG_ID']
                                mets_df.loc[idx, 'Smiles'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
                            except:
                                # old
                                # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                # continue
                                # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                                # can be sanitizable and pass through the inchi, thus, if present, put the 
                                # SMILES of metanetx
                                p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                                pp = p.text
                                try:
                                    mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                    if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                        mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                        continue
                                except:
                                    mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                    continue
                            try:
                                sanit = Chem.SanitizeMol(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles']))
                                check_inchi = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles'])))
                                if check_inchi == None:
                                    # old
                                    # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                    # continue
                                    # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                                    # can be sanitizable and pass through the inchi, thus, if present, put the 
                                    # SMILES of metanetx
                                    p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                                    pp = p.text
                                    try:
                                        mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                        if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                            mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                            continue
                                    except:
                                        mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                        continue
                            except:
                                # old
                                # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                # continue
                                # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                                # can be sanitizable and pass through the inchi, thus, if present, put the 
                                # SMILES of metanetx
                                p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                                pp = p.text
                                try:
                                    mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                    if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                        mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                        continue
                                except:
                                    mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                    continue
                    except:
                        try:
                            id = mets_df.loc[idx,'KEGG_ID']
                            mets_df.loc[idx, 'Smiles'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
                        except:
                            # old
                            # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                            # continue
                            # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                            # can be sanitizable and pass through the inchi, thus, if present, put the 
                            # SMILES of metanetx
                            p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                            pp = p.text
                            try:
                                mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                    mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                    continue
                            except:
                                mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                continue
                        try:
                            sanit = Chem.SanitizeMol(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles']))
                            check_inchi = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles'])))
                            if check_inchi == None:
                                # old
                                # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                # continue
                                # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                                # can be sanitizable and pass through the inchi, thus, if present, put the 
                                # SMILES of metanetx
                                p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                                pp = p.text
                                try:
                                    mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                    if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                        mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                        continue
                                except:
                                    mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                    continue
                        except:
                            # old
                            # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                            # continue
                            # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                            # can be sanitizable and pass through the inchi, thus, if present, put the 
                            # SMILES of metanetx
                            p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                            pp = p.text
                            try:
                                mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                    mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                    continue
                            except:
                                mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                continue
            except:
                id = mets_df.loc[idx,'Name']
                compound = pcp.get_compounds(id, 'name')
                if len(compound) == 0:
                    try:
                        id = mets_df.loc[idx,'KEGG_ID']
                        mets_df.loc[idx, 'Smiles'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
                    except:
                        # old
                        # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                        # continue
                        # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                        # can be sanitizable and pass through the inchi, thus, if present, put the 
                        # SMILES of metanetx
                        p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                        pp = p.text
                        try:
                            mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                            if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                continue
                        except:
                            mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                            continue
                else:
                    mets_df.loc[idx, 'Smiles'] = compound[0].canonical_smiles
                try:
                    sanit = Chem.SanitizeMol(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles']))
                    check_inchi = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles'])))
                    if check_inchi == None:
                        try:
                            id = mets_df.loc[idx,'KEGG_ID']
                            mets_df.loc[idx, 'Smiles'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
                        except:
                            # old
                            # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                            # continue
                            # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                            # can be sanitizable and pass through the inchi, thus, if present, put the 
                            # SMILES of metanetx
                            p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                            pp = p.text
                            try:
                                mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                    mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                    continue
                            except:
                                mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                continue 
                        try:
                            sanit = Chem.SanitizeMol(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles']))
                            check_inchi = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles'])))
                            if check_inchi == None:
                                # old
                                # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                # continue
                                # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                                # can be sanitizable and pass through the inchi, thus, if present, put the 
                                # SMILES of metanetx
                                p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                                pp = p.text
                                try:
                                    mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                    if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                        mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                        continue
                                except:
                                    mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                    continue
                        except:
                            # old
                            # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                            # continue
                            # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                            # can be sanitizable and pass through the inchi, thus, if present, put the 
                            # SMILES of metanetx
                            p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                            pp = p.text
                            try:
                                mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                    mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                    continue
                            except:
                                mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                continue
                except:
                    try:
                        id = mets_df.loc[idx,'KEGG_ID']
                        mets_df.loc[idx, 'Smiles'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
                    except:
                        # old
                        # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                        # continue
                        # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                        # can be sanitizable and pass through the inchi, thus, if present, put the 
                        # SMILES of metanetx
                        p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                        pp = p.text
                        try:
                            mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                            if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                continue
                        except:
                            mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                            continue
                    try:
                        sanit = Chem.SanitizeMol(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles']))
                        check_inchi = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles'])))
                        if check_inchi == None:
                            # old
                            # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                            # continue
                            # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                            # can be sanitizable and pass through the inchi, thus, if present, put the 
                            # SMILES of metanetx
                            p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                            pp = p.text
                            try:
                                mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                    mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                    continue
                            except:
                                mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                continue
                    except:
                        # old
                        # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                        # continue
                        # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                        # can be sanitizable and pass through the inchi, thus, if present, put the 
                        # SMILES of metanetx
                        p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                        pp = p.text
                        try:
                            mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                            if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                continue
                        except:
                            mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                            continue
        except:
            id = mets_df.loc[idx,'Name']
            compound = pcp.get_compounds(id, 'name')
            if len(compound) == 0:
                try:
                    id = mets_df.loc[idx,'KEGG_ID']
                    mets_df.loc[idx, 'Smiles'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
                except:
                    # old
                    # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                    # continue
                    # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                    # can be sanitizable and pass through the inchi, thus, if present, put the 
                    # SMILES of metanetx
                    p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                    pp = p.text
                    try:
                        mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                        if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                            mets_df.loc[idx,'Smiles'] == 'no_smiles'
                            continue
                    except:
                        mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                        continue
            else:
                mets_df.loc[idx, 'Smiles'] = compound[0].canonical_smiles
            try:
                sanit = Chem.SanitizeMol(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles']))
                check_inchi = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles'])))
                if check_inchi == None:
                    try:
                        id = mets_df.loc[idx,'KEGG_ID']
                        mets_df.loc[idx, 'Smiles'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
                    except:
                        # old
                        # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                        # continue
                        # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                        # can be sanitizable and pass through the inchi, thus, if present, put the 
                        # SMILES of metanetx
                        p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                        pp = p.text
                        try:
                            mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                            if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                continue
                        except:
                            mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                            continue 
                    try:
                        sanit = Chem.SanitizeMol(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles']))
                        check_inchi = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles'])))
                        if check_inchi == None:
                            # old
                            # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                            # continue
                            # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                            # can be sanitizable and pass through the inchi, thus, if present, put the 
                            # SMILES of metanetx
                            p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                            pp = p.text
                            try:
                                mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                                if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                    mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                    continue
                            except:
                                mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                                continue
                    except:
                        # old
                        # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                        # continue
                        # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                        # can be sanitizable and pass through the inchi, thus, if present, put the 
                        # SMILES of metanetx
                        p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                        pp = p.text
                        try:
                            mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                            if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                continue
                        except:
                            mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                            continue
            except:
                try:
                    id = mets_df.loc[idx,'KEGG_ID']
                    mets_df.loc[idx, 'Smiles'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
                except:
                    # old
                    # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                    # continue
                    # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                    # can be sanitizable and pass through the inchi, thus, if present, put the 
                    # SMILES of metanetx
                    p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                    pp = p.text
                    try:
                        mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                        if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                            mets_df.loc[idx,'Smiles'] == 'no_smiles'
                            continue
                    except:
                        mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                        continue
                try:
                    sanit = Chem.SanitizeMol(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles']))
                    check_inchi = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(mets_df.loc[idx, 'Smiles'])))
                    if check_inchi == None:
                        # old
                        # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                        # continue
                        # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                        # can be sanitizable and pass through the inchi, thus, if present, put the 
                        # SMILES of metanetx
                        p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                        pp = p.text
                        try:
                            mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                            if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                                mets_df.loc[idx,'Smiles'] == 'no_smiles'
                                continue
                        except:
                            mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                            continue
                except:
                    # old
                    # mets_df.loc[idx,'Smiles'] == 'no_smiles'
                    # continue
                    # CORRECT: at this point, none of the IDs are useful to extract a SMILES that 
                    # can be sanitizable and pass through the inchi, thus, if present, put the 
                    # SMILES of metanetx
                    p = requests.request('GET',metanetx_server+mets_df.loc[idx,'METANETX_ID'])
                    pp = p.text
                    try:
                        mets_df.loc[idx, 'Smiles'] = re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0]
                        if mets_df.loc[idx, 'Smiles'] == '&nbsp;':
                            mets_df.loc[idx,'Smiles'] == 'no_smiles'
                            continue
                    except:
                        mets_df.loc[idx, 'Smiles'] = 'no_smiles'
                        continue
        mets_df.loc[idx, 'Smiles'] = mets_df.loc[idx, 'Smiles'].replace('.[Cl-]','')
        mets_df.loc[idx, 'Smiles'] = mets_df.loc[idx, 'Smiles'].replace('.[Co]','')
        mets_df.loc[idx, 'Smiles'] = mets_df.loc[idx, 'Smiles'].replace('[Fe++].','')
    
    mets_dictionary = mets_df.set_index('METANETX_ID')['Smiles'].to_dict()
    with open('./output/mets_dictionary.pkl', 'wb') as f:
        pickle.dump([mets_df,mets_dictionary], f)
            
##############################################################################

#Modify the EC numbers where we have less than 4 digits (- included) adding 
#the -

complete_digits = rules['EC number'].str.contains('\d+.\d+.\d+.\d+')

dictionary_tmp = {False: True, True: False, None: False}

not_complete_digits = complete_digits.apply(lambda x: dictionary_tmp.get(x))

rules['EC number'] = rules['EC number'][not_complete_digits] +'.-'

##############################################################################

#Load the reactions of interest and compare with the dataframe just created

rxns_library = pd.read_csv('./output/reduced_reactions_bounds.csv')

prova_lb = rxns_library.copy()
prova_lb.drop(['Names','EC','Formula','Taxonomy','ub'],axis=1,inplace=True)
dictionary_id_lb = prova_lb.set_index('ID').T.to_dict('records')

prova_ub = rxns_library.copy()
prova_ub.drop(['Names','EC','Formula','Taxonomy','lb'],axis=1,inplace=True)
dictionary_id_ub = prova_ub.set_index('ID').T.to_dict('records')

##############################################################################
#Add EC number to reactions of agora
rules_explicit = rules.loc[:,['Rule_ID','Rule_SMARTS_EXP','Reaction_ID','Diameter','Direction','Substrate_ID','EC number','Score','seed id']]
rules_explicit.rename(columns = {'Rule_SMARTS_EXP': 'Rule_SMARTS', 'Score': 'Score_normalized', 'EC number':'Reaction_EC_number'}, inplace = True)

if os.path.isfile("./output/rxns_library_with_EC.pkl"):
    # Getting back the objects:
    with open('./output/rxns_library_with_EC.pkl', 'rb') as f: 
        rxns_library = pickle.load(f)
else:
    for i in range(len(rxns_library['EC'])):
        if type(rxns_library['EC'][i]) is float:
            print(str(i)+'/'+str(len(rxns_library)-1))
            #enzyme = name of the rules file where the space is = \s
            try:
                enzyme_name = str.replace(rxns_library['Names'][i], ' ', '%20')
                #Extract the EC numbers related to the reactions with urllib
                url = kegg+enzyme_name
                html = urllib.request.urlopen(url)
                soup = BeautifulSoup(html)
                text = soup.get_text()
                rxnID = re.findall('\w+:(R\d+)\s',text)[0]
                text = k.get(rxnID)
                ec_extracted = re.findall('EC:(\d+.\d+.\d+.\d+)',text)
                rxns_library['EC'][i] = '|'.join(ec_extracted)
            except:
                continue
    with open('./output/rxns_library_with_EC.pkl', 'wb') as f:
        pickle.dump(rxns_library, f)
        
#In the EC number of rxns library we can have more EC for a reaction
#They are write in the same line separated with a "|", so for the comparison
#Create a list with each EC in a line
    
# EC_in_model = []
# for i in rxns_library['EC']:
#     try:
#         EC_in_model += i.split('|')
#     except:
#         continue
    
#In EC_in_model variable we have all the EC numbers present in rxns_library

##############################################################################

EC_x_tax = pd.read_excel('./input/EC_x_tax.xlsx')

df1_explicit = rules_explicit.loc[rules_explicit['Reaction_EC_number'].isin(pd.unique(EC_x_tax.EC)),]

df2_explicit = rules_explicit.loc[rules_explicit['seed id'].isin(pd.unique(rxns_library['ID'].values)),]

df3_explicit = df1_explicit.merge(df2_explicit,how='outer')

#Load the information about annotated EC numbers from NCBI genomes and KEGG
#Collect them and extract the rule with those EC numbers

path_peg = './input/new_taxonomy_peg'
path_kegg = './input/species_level_extract/result'
ec_list = []

for filename in os.listdir(path_peg):
    with open(path_peg+'/'+filename,'r') as f:
        text = f.read()
        temp = re.findall('[ET]C[ :][\w\.\-]+',text)
        tmp = [re.findall('[\w.\-]+',x)[1] for x in temp]
        ec_list += tmp
        
for filename in os.listdir(path_kegg):
    with open(path_kegg+'/'+filename,'r') as f:
        text = f.read()
        tmp = re.findall('[\w.\-]+',text)
        ec_list += tmp
        
final_list_ec = np.unique(ec_list)

rules_explicit_tmp = rules_explicit.loc[rules_explicit['Reaction_EC_number'].isin(final_list_ec),]

df4_explicit = df3_explicit.merge(rules_explicit_tmp,how='outer')   
df4_explicit.loc[:,'Substrate_SMILES'] = None
df4_explicit.loc[:,'Product_SMILES'] = None 

# load the rules for sergio's reactions
sergio_rules = pd.read_excel('./output/sergio_explicit_rules.xlsx')

df5_explicit = df4_explicit.merge(sergio_rules,how='outer')

df5_explicit.drop(df5_explicit.index[df5_explicit.loc[:,'Diameter'] < 6], inplace = True)

df5_explicit.drop_duplicates(subset = ['Rule_ID','Rule_SMARTS','Reaction_ID','Direction','Substrate_ID', 'Score_normalized'], inplace = True)

df5_explicit.drop_duplicates(subset = ['Rule_ID','Diameter','Reaction_ID','Direction','Substrate_ID'], inplace = True)

df5_explicit.Rule_ID = [df5_explicit.Rule_ID[x] + '_' + str(df5_explicit.Diameter[x]) for x in df5_explicit.index]
df5_explicit.reset_index(drop = True, inplace=True)

################### Add the information about substrate's and product's SMILES

index_rxnAdd1 = df5_explicit.index[df5_explicit['Reaction_ID'].isin(['rxnAdd1'])][0]

index_of_interest = range(0,index_rxnAdd1)
    
    
for idx in index_of_interest:
    
    print(str(idx)+'/'+str(len(index_of_interest)))
    df5_explicit.loc[idx,'Substrate_SMILES'] = mets_dictionary[df5_explicit.loc[idx,'Substrate_ID']]
    if df5_explicit.loc[idx,'Substrate_ID'] in dictionary_reac_subs[df5_explicit.loc[idx,'Reaction_ID']]:
        products_id = dictionary_reac_prods[df5_explicit.loc[idx,'Reaction_ID']]
        df5_explicit.loc[idx,'Product_SMILES'] = '.'.join([mets_dictionary[x] for x in products_id])
        df5_explicit.loc[idx,'Product_IDs'] = '.'.join(products_id)
    elif df5_explicit.loc[idx,'Substrate_ID'] in dictionary_reac_prods[df5_explicit.loc[idx,'Reaction_ID']]:
        products_id = dictionary_reac_subs[df5_explicit.loc[idx,'Reaction_ID']]
        df5_explicit.loc[idx,'Product_SMILES'] = '.'.join([mets_dictionary[x] for x in products_id])
        df5_explicit.loc[idx,'Product_IDs'] = '.'.join(products_id)
        
# Check if we have problem to sanitize the structure of molecules. If yes, take
# the smiles from other databases NEEDED for the algorithm

no_san = []
for x in df5_explicit.index:
    try:
        Chem.SanitizeMol(Chem.MolFromSmiles(df5_explicit.loc[x,'Substrate_SMILES']))
    except:
        no_san.append(x)
        
no_san_1 = []
for x in no_san:
    tmp_id = df5_explicit.loc[x,'Substrate_ID']
    try:
        id = mets_df.loc[mets_df['METANETX_ID'].isin([tmp_id]),'KEGG_ID'].values[0]
        df5_explicit.loc[x,'Substrate_SMILES'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
    except:
        no_san_1.append(x)
        
# Check if we have problem to pass through the inchi of molecules. If yes, take
# the smiles from other databases NEEDED for the algorithm      
 
no_inchi = []
for x in df5_explicit.index:
    a = Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(df5_explicit.loc[x,'Substrate_SMILES'])))
    if a == None:
        no_inchi.append(x)

no_inchi_1 = []
for x in no_inchi:
    tmp_id = df5_explicit.loc[x,'Substrate_ID']
    id = mets_df.loc[mets_df['METANETX_ID'].isin([tmp_id]),'Name'].values[0]
    compound = pcp.get_compounds(id, 'name')
    if len(compound) == 0:
        try:
            id = mets_df.loc[mets_df['METANETX_ID'].isin([tmp_id]),'KEGG_ID'].values[0]
            df5_explicit.loc[x,'Substrate_SMILES'] = Chem.MolToSmiles(Chem.MolFromMolBlock(k.get(id, option = 'mol')))
        except:
            no_inchi_1.append(x)
    else:
        df5_explicit.loc[x,'Substrate_SMILES'] = compound[0].canonical_smiles
        
        
prob_prod = []
for x in df5_explicit.index:
    ids = df5_explicit.loc[x,'Product_IDs'].split('.')
    smiles = df5_explicit.loc[x,'Product_SMILES'].split('.')
    if 'no_smiles' in smiles:
        prob_prod.append(x)
  
res_prod = []
for x in np.unique(prob_prod):
    p = requests.request('GET',metanetx_server+x)
    pp = p.text
    try:
        res_prod.append(re.findall('<tr><td>SMILES</td><td class="smiles" colspan="3">(\S+)</td></tr>',pp)[0])
    except:
        res_prod.append(str(x)+'_no_smiles')
        
Rule_x_ID = pd.DataFrame({'RuleID': ['_'.join(x.split('_')[0:2]) for x in df5_explicit['Rule_ID']], 'rxnID': df5_explicit['seed id']})
Rule_x_ID.drop_duplicates(inplace = True)
# Rule_x_ID.to_excel('./output/Rule_x_ID.xlsx', index=False)
    
df6_explicit = df5_explicit.drop(prob_prod)

# Add the information about the bounds for those reactions we know that
lb = [dictionary_id_lb[0][x] if x and x in dictionary_id_lb[0].keys() else None for x in df6_explicit['seed id'].values]
ub = [dictionary_id_ub[0][x] if x and x in dictionary_id_ub[0].keys() else None for x in df6_explicit['seed id'].values]

df6_explicit['lb'] = lb
df6_explicit['ub'] = ub

# df6_explicit.to_csv('./output/rules/rules_with_H_products_bounds.csv', sep= '\t', index= False)    
# df6_explicit.to_csv('./output/rules/rules_without_H_products_bounds.csv', sep= '\t', index= False) 

df6_explicit.drop(['seed id', 'lb', 'ub'],axis=1,inplace=True)

# Check manually the None values or the presence of &nbsp;
# Remove the [Fe++] and [Cl-] from Substrate_SMILES
# df6_explicit.to_csv('./output/rules/rules_with_H_products.csv', sep= '\t', index= False)    
# df6_explicit.to_csv('./output/rules/rules_without_H_products.csv', sep= '\t', index= False) 
   
# Add the information about the bounds for those reactions we know that
lb = [dictionary_id_lb[0][x] if x and x in dictionary_id_lb[0].keys() else None for x in df5_explicit['seed id'].values]
ub = [dictionary_id_ub[0][x] if x and x in dictionary_id_ub[0].keys() else None for x in df5_explicit['seed id'].values]

df5_explicit['lb'] = lb
df5_explicit['ub'] = ub

# df5_explicit.to_csv('./output/rules/rules_with_H_only_substrate_bounds.csv', sep= '\t', index= False)    
# df5_explicit.to_csv('./output/rules/rules_without_H_only_substrate_bounds.csv', sep= '\t', index= False) 

# Remove the information about the product IDs and Smiles
df5_explicit.drop(['seed id','lb','ub','Product_IDs', 'Product_SMILES'],axis=1,inplace=True)

# Check manually the None values or the presence of &nbsp;
# Remove the [Fe++] and [Cl-] from Substrate_SMILES
# df5_explicit.to_csv('./output/rules/rules_with_H_only_substrate.csv', sep= '\t', index= False)    
# df5_explicit.to_csv('./output/rules/rules_without_H_only_substrate.csv', sep= '\t', index= False)    

##############################################################################

# Generate the formula-rules information

df4_explicit = pd.read_csv('./output/rules/rules_with_H_only_substrate_bounds.csv',sep='\t')

import requests
metanetx_server = 'https://www.metanetx.org/equa_info/'

index_rxnAdd1 = df5_explicit.index[df5_explicit['Reaction_ID'].isin(['rxnAdd1'])][0]

index_of_interest = range(0,index_rxnAdd1)

names = df4_explicit['Rule_ID'][index_of_interest].values
names_reactions = [x.split('_')[0] for x in names]
names_reactions = np.unique(names_reactions)

formulas = pd.DataFrame(index = range(len(names_reactions)),columns = ['id','formulas','lb','ub'])
formulas.loc[:,'id'] = names_reactions

for i in range(len(formulas)):
    
    print(str(i)+'/'+str(len(formulas)-1))
    
    id = formulas['id'][i]

    p = requests.request('GET',metanetx_server+id)
    pp = p.text
    try:
        text_equation = re.findall("<tr><td>equation</td><td><div class='comp_full'><span class='mono'>([\S\s]+)<em>generic compartment 1</em>",pp)[0]
        ppp = re.findall(r"<a href='/chem_info/MNXM\d+' class='chem'>(.*)</a>",text_equation)[0]
        
        pppp = ppp.split('&nbsp;&lt;=&gt;&nbsp;')
        
        substrates = pppp[0].split(' + ')
        products = pppp[1].split(' + ')
        
        substrates[0] = substrates[0].split('<')[0]
        
        for j in range(1,len(substrates)):
            substrates[j] = re.findall(r"class='chem'>(.*?)</a>",substrates[j])[0]
            
        products[-1] = re.findall("class='chem'>(.*?)$",products[-1])[0]
        
        for j in range(0,len(products)-1):
            products[j] = re.findall(r"class='chem'>(.*?)</a>",products[j])[0]
    except:
        try:
            text_equation = re.findall("<tr><td>equation</td>([\S\s]+)<em>generic compartment \d</em>",pp)[0]
            ppp = re.findall(r"<a href='/chem_info/MNXM\d+' class='chem \S+'>(.*)</a>",text_equation)[0]
            
            pppp = ppp.split('&nbsp;&lt;=&gt;&nbsp;')
                
            substrates = pppp[0].split(' + ')
            products = pppp[1].split(' + ')
            
            substrates[0] = substrates[0].split('<')[0]
            
            for j in range(1,len(substrates)):
                substrates[j] = re.findall(r"class='chem'>(.*?)</a>",substrates[j])[0]
            
            products[-1] = re.findall(r"'>(\S+)$",products[-1])[0]
            
            for j in range(0,len(products)-1):
                products[j] = re.findall(r"class='chem'>(.*?)</a>",products[j])[0]
        except:
            try:
                text_equation = re.findall("<tr><td>equation</td>([\S\s]+)<em>generic compartment \d</em>",pp)[0]
                ppp = re.findall(r"<a href='/chem_info/MNXM\d+' class='chem'>(.*)</a>",text_equation)[0]
                
                pppp = ppp.split('&nbsp;&lt;=&gt;&nbsp;')
                
                substrates = pppp[0].split(' + ')
                products = pppp[1].split(' + ')
                
                substrates[0] = substrates[0].split('<')[0]
                
                for j in range(1,len(substrates)):
                    substrates[j] = re.findall(r"class='chem'>(.*?)</a>",substrates[j])[0]
                    
                products[-1] = re.findall("class='chem'>(.*?)$",products[-1])[0]
                
                for j in range(0,len(products)-1):
                    products[j] = re.findall(r"class='chem'>(.*?)</a>",products[j])[0]
            except:
                try:
                    text_equation = re.findall("<tr><td>equation</td>([\S\s]+)<em>generic compartment \d</em>",pp)[0]
                    ppp = re.findall(r"<a href='/chem_info/MNXM\d+' class='chem'>(.*)</a>",text_equation)[0]
                    
                    pppp = ppp.split('&nbsp;&lt;=&gt;&nbsp;')
                    
                    substrates = pppp[0].split(' + ')
                    products = pppp[1].split(' + ')
                    
                    substrates[0] = substrates[0].split('<')[0]
                    
                    for j in range(1,len(substrates)):
                        try:
                            substrates[j] = re.findall(r"class='chem \S+'>(.*?)</a>",substrates[j])[0]
                        except:
                            substrates[j] = re.findall(r"class='chem'>(.*?)</a>",substrates[j])[0]
                            
                    products[-1] = re.findall("class='chem'>(.*?)$",products[-1])[0]
                    
                    for j in range(0,len(products)-1):
                        try:
                            products[j] = re.findall(r"class='chem'>(.*?)</a>",products[j])[0]
                        except:
                            products[j] = re.findall(r"class='chem \S+'>(.*?)</a>",products[j])[0]
                except:
                    continue
    first = ' + '.join(substrates)
    second = ' + '.join(products)
    
    full = ' <=> '.join([first,second])
    full = full.replace('&#39;',"'")
    formulas.loc[i,'formulas'] = full
    if np.isnan(np.unique(df4_explicit.loc[df4_explicit['Reaction_ID'].isin([id]),'lb'])[0]):
        formulas.loc[i,'lb'] = -1000
    else:
        formulas.loc[i,'lb'] = np.unique(df4_explicit.loc[df4_explicit['Reaction_ID'].isin([id]),'lb'])[0]
        
    if np.isnan(np.unique(df4_explicit.loc[df4_explicit['Reaction_ID'].isin([id]),'ub'])[0]):
        formulas.loc[i,'ub'] = 1000
    else:
        formulas.loc[i,'ub'] = np.unique(df4_explicit.loc[df4_explicit['Reaction_ID'].isin([id]),'ub'])[0]

sergio_reactions = pd.read_excel('./output/sergio_reactions_with_smiles.xlsx')
sergio_reactions = pd.DataFrame({'id': sergio_reactions.rxnID, 'formulas': sergio_reactions.Formula, 'lb': 0, 'ub': 1000})
formulas = formulas.merge(sergio_reactions, how="outer")
# formulas.to_excel('./output/formulas_rules.xlsx', index=False)

