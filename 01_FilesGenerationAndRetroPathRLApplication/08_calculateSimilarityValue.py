# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 14:13:11 2021

@author: fbalzerani
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np

phenol_explorer = pd.read_csv('./output/general_source_phen_ex.csv', sep='\t',header = None, names=['Name','InChI'])

m_sources = [Chem.MolFromInchi(x) for x in phenol_explorer.InChI]
f_sources = [AllChem.GetMorganFingerprint(x,2) for x in m_sources]

rules = pd.read_csv("./output/rules/rules_with_H_only_substrate.csv", sep = '\t')
m_rules = [Chem.MolFromSmiles(x) for x in rules.Substrate_SMILES]
f_rules = [AllChem.GetMorganFingerprint(x,2) for x in m_rules]

similarity = []
for source in f_sources:
    for rule in f_rules:
        similarity.append(DataStructs.DiceSimilarity(source,rule))

similarity_df = pd.DataFrame(data = {'Values' : similarity})
# similarity_df.to_csv("./output/similarity_values.csv", index = False)
        
# p = np.array(similarity)
# len(p[p>0.6])/len(p)*100
