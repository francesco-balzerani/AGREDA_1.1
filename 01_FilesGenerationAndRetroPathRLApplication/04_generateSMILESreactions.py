# -*- coding: utf-8 -*-
"""
Created on Wed May 19 13:52:51 2021

@author: fbalzerani
"""

import pandas as pd
import pubchempy as pcp
import re

# Code to generate the formulas by the SMILES

reactions = pd.read_excel("./output/sergio_reactions_to_generate_rules.xlsx")
cofactors = pd.read_csv('./output/cofactors_inchi.csv', sep= ',')
cofactors = cofactors.loc[:,'name']
cofactors[cofactors.index[-1]+1] = 'Methane'
cofactors[cofactors.index[-1]+1] = 'Chloride'
cofactors[cofactors.index[-1]+1] = 'Nicotinamide adenine dinucleotide'
cofactors[cofactors.index[-1]+1] = 'Nicotinamide adenine dinucleotide - reduced'

new_reactions = pd.DataFrame(columns = ['rxnID','Formula'])
count= 0
for i in reactions.index:
    print(str(i)+'/'+str(reactions.index[-1]))
    substrate = reactions.loc[i,'Formula'].split(' -> ')[0].rstrip().lstrip().split(' + ')
    products = reactions.loc[i,'Formula'].split(' -> ')[1].rstrip().lstrip().split(' + ')
    
    new_substrates = [x for x in substrate if x not in cofactors.values]
    for sub in new_substrates:
        new_reactions.loc[count,'rxnID'] = reactions.loc[i,'rxnID']
        new_reactions.loc[count,'Formula'] = ' -> '.join([sub,' + '.join(products)])
        count +=1
    
    new_products = [x for x in products if x not in cofactors.values]
    n_prod = 0
    for prod in new_products:
        n_prod += 1
        new_reactions.loc[count,'rxnID'] = reactions.loc[i,'rxnID']+'_r_'+str(n_prod)
        new_reactions.loc[count,'Formula'] = ' -> '.join([prod,' + '.join(substrate)])
        count +=1    
    
not_done = []

for i in new_reactions.index:
    print(str(i)+'/'+str(new_reactions.index[-1]))
    
    substrate = new_reactions.loc[i,'Formula'].split(' -> ')[0].rstrip().lstrip().split(' + ')
    for j in range(len(substrate)):
        if re.findall('[\d\.]+ ', substrate[j]):
            substrate[j] = substrate[j].replace(re.findall('[\d\.]+ ',substrate[j])[0],"")
        if 'IdB 1027' in substrate[j]:
            substrate[j] = substrate[j].replace('IdB 1027',"Cyanidin")
    
    products = new_reactions.loc[i,'Formula'].split(' -> ')[1].rstrip().lstrip().split(' + ')
    
    for j in range(len(products)):
        if re.findall('[\d\.]+ ',products[j]):
            products[j] = products[j].replace(re.findall('[\d\.]+ ',products[j])[0],"")
        if 'IdB 1027' in products[j]:
            products[j] = products[j].replace('IdB 1027',"Cyanidin")
    try:        
        compound_s = [pcp.get_compounds(x, 'name')[0].canonical_smiles for x in substrate]
        formula_s = '.'.join(compound_s)
    except:
        not_done.append(i)
        next
    try:
        compound_p = [pcp.get_compounds(x, 'name')[0].canonical_smiles for x in products]
        formula_p = '.'.join(compound_p)
    except:
        not_done.append(i)
        next
        
    new_reactions.loc[i,'Smiles'] = '>>'.join([formula_s,formula_p])
    
for i in not_done:
    print(str(i)+'/'+str(new_reactions.index[-1]))
    
    substrate = new_reactions.loc[i,'Formula'].split(' -> ')[0].rstrip().lstrip().split(' + ')
    for j in range(len(substrate)):
        if re.findall('[\d\.]+ ', substrate[j]):
            substrate[j] = substrate[j].replace(re.findall('[\d\.]+ ',substrate[j])[0],"")
        if 'IdB 1027' in substrate[j]:
            substrate[j] = substrate[j].replace('IdB 1027',"Cyanidin")
    
    products = new_reactions.loc[i,'Formula'].split(' -> ')[1].rstrip().lstrip().split(' + ')
    
    for j in range(len(products)):
        if re.findall('[\d\.]+ ',products[j]):
            products[j] = products[j].replace(re.findall('[\d\.]+ ',products[j])[0],"")
        if 'IdB 1027' in products[j]:
            products[j] = products[j].replace('IdB 1027',"Cyanidin")
    try:        
        compound_s = [pcp.get_compounds(x, 'name')[0].canonical_smiles for x in substrate]
        formula_s = '.'.join(compound_s)
        new_reactions.loc[i,'Smiles'] = formula_s+'>>'
    except:
        compound_p = [pcp.get_compounds(x, 'name')[0].canonical_smiles for x in products]
        formula_p = '.'.join(compound_p)
        new_reactions.loc[i,'Smiles'] = '>>'+formula_p

# Check the chloride .[Cl-] manually
# new_reactions.to_excel('./output/sergio_reactions_with_smiles.xlsx', index = False)