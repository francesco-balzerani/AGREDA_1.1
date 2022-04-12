# -*- coding: utf-8 -*-
"""
Created on Fri May 21 10:39:37 2021

@author: fbalzerani
"""
       
import os
import pandas as pd
import json

path = './input/sergio_json/'

list_mets = pd.read_excel('./output/cpd_name_sergio_reactions.xlsx')
list_mets = list_mets.set_index('Name')['cpdID'].to_dict()

list_reactions = pd.read_excel('./output/sergio_reactions_with_smiles.xlsx')

new_rules_explicit = pd.DataFrame(columns = ['Rule_ID', 'Rule_SMARTS','Reaction_ID','Diameter','Direction','Substrate_ID','Reaction_EC_number','Score_normalized', 'Substrate_SMILES', 'Product_IDs','Product_SMILES'])
files = os.listdir(path)

pos = 0

for file in files:
    rxnID = file.split('.')[0]
    tmp = list_reactions.loc[list_reactions['rxnID'].isin([rxnID]),'Formula'].values[0].split(' -> ')[0]
    metID = tmp.replace(tmp,list_mets[tmp])
    Rule_ID = '_'.join([rxnID,metID])
    prodID = list_reactions.loc[list_reactions['rxnID'].isin([rxnID]),'Formula'].values[0].split(' -> ')[1].split(' + ')
    prodID = '.'.join([x.replace(x,list_mets[x]) for x in prodID])
    
    # Load the json
    with open(path+file, "r") as jsonFile:
        data = json.load(jsonFile)
        
    counts = 0
    
    for rule in data:
        if rule['direction'] == 'Reverse':
            continue
        counts += 1
        new_rules_explicit.loc[pos,'Rule_ID'] = Rule_ID + '_' + str(counts)
        new_rules_explicit.loc[pos,'Rule_SMARTS'] = rule['smarts_string']
        new_rules_explicit.loc[pos,'Reaction_ID'] = rxnID
        new_rules_explicit.loc[pos,'Diameter'] = rule['diameter']
        new_rules_explicit.loc[pos, 'Direction'] = 1
        new_rules_explicit.loc[pos,'Substrate_ID'] = metID
        new_rules_explicit.loc[pos,'Reaction_EC_number'] = rule['ec_numbers']
        new_rules_explicit.loc[pos,'Score_normalized'] = rule['score']
        new_rules_explicit.loc[pos,'Substrate_SMILES'] = list_reactions.loc[list_reactions['rxnID'].isin([rxnID]),'Smiles'].values[0].split('>>')[0]
        new_rules_explicit.loc[pos,'Product_IDs'] = prodID
        new_rules_explicit.loc[pos,'Product_SMILES'] = list_reactions.loc[list_reactions['rxnID'].isin([rxnID]),'Smiles'].values[0].split('>>')[1]
        pos += 1

# new_rules_explicit.to_excel('./output/sergio_explicit_rules.xlsx', index = False)


