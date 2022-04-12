# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:07:46 2020

@author: fbalzerani
"""

"""
Generate the excel files to add the phenols and reactions to the network:
    - an excel for the metabolites, generating the metID and extracting the 
        metNames and metInChI from the information used to apply RetroPath.
        The mets are written manually at the end of this script.
        
    - an excel for the reactions, where the lb and ub values are defined based
        on our previous knowledge (equal to the lb and ub of the template 
        reactions if we know these informations or otherwise equal to -1000 
        and 1000) the EC number, taxonomy, the equation based on the cpdID and 
        the stoichiometry are extracted from the table results of RetroPath.
        For the taxonomy we generate both trules and trRules fields.
        The stoichiometry is generated using the chempy package.

"""

# packages

from chempy import balance_stoichiometry
import pandas as pd
import numpy as np
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import itertools

# load the table where are the equations and the list of polyphenols
table_results_known = pd.read_excel('./output/final_table_CS.xlsx', engine = "openpyxl")
phenol_explorer = pd.read_csv('../01_FilesGenerationAndRetroPathRLApplication/output/general_source_phen_ex.csv', sep='\t',header = None, names=['Name','InChI'])

names_cpd = pd.read_excel('../01_FilesGenerationAndRetroPathRLApplication/output/names_cpd_model_merged_sergio.xlsx',
                          engine = "openpyxl")
names_cpd.rename(columns={'ID':'metID','Names':'metNames'},inplace=True)

table_results_known.dropna(subset=['Formula_mets'], inplace=True)
table_results_known.index = range(len(table_results_known))

# generate the excel to add the molecules to the model
unique_phenols = table_results_known['Name_metabolite'].dropna()
unique_phenols.drop_duplicates(inplace=True)
unique_phenols.index = range(len(unique_phenols))
unique_phenols = [x.replace('_',' ') for x in unique_phenols]

new_phenols = pd.DataFrame(index = range(len(unique_phenols)), data={'metID':"", 'metNames': unique_phenols, 'metInChI': ""})
new_phenols.metInChI = [phenol_explorer.loc[phenol_explorer['Name'].isin([x]),'InChI'].values[0] for x in new_phenols.metNames]

##############################################################################
##############################################################################
############################################################################## 
# Run this part manually in order to extract those compounds that can be already
# in the universal database but not in AGREDA, extracting those cpd and modifying
# manually. Needed to do manually because of the similarity between sugars that
# can bring us to errors.

m_sources = [Chem.MolFromInchi(x) for x in new_phenols.metInChI]
f_sources = [AllChem.GetMorganFingerprint(x,2) for x in m_sources]

orig_sink = pd.read_csv("../01_FilesGenerationAndRetroPathRLApplication/output/general_sink.csv")
sink = orig_sink.copy()
m_sink = []
none_v = []
for x in sink.index:
    if Chem.MolFromInchi(sink.inchi[x]) == None:
        none_v.append(x)
    else:
        m_sink.append(Chem.MolFromInchi(sink.inchi[x]))
        
sink.drop(none_v, inplace=True)
sink.index = range(len(sink))

f_sink = [AllChem.GetMorganFingerprint(x,2) for x in m_sink]

similarity_sink = pd.DataFrame(index = range(len(new_phenols)), data = {'source':new_phenols.metNames})
for x in range(len(f_sources)):
    for y in range(len(f_sink)):
        if DataStructs.DiceSimilarity(f_sources[x],f_sink[y]) == 1: 
            similarity_sink.loc[x, 'sink'] = sink.name[y]

# After the manual run, annotate here the changes to make
new_phenols.loc[new_phenols['metNames'].isin(["Cyanidin 3-O-(3'',6''-O-dimalonyl-glucoside)"]),'metID'] = 'cpd15005'
new_phenols.loc[new_phenols['metNames'].isin(["5,6,7,4'-Tetrahydroxyisoflavone"]),'metID'] = 'cpd28833'
new_phenols.loc[new_phenols['metNames'].isin(["6,7,4'-Trihydroxyisoflavone"]),'metID'] = 'cpd10013'
new_phenols.loc[new_phenols['metNames'].isin(["Butin"]),'metID'] = 'cpd06509'
new_phenols.loc[new_phenols['metNames'].isin(["Calycosin"]),'metID'] = 'cpd01100'
new_phenols.loc[new_phenols['metNames'].isin(["Protocatechuic aldehyde"]),'metID'] = 'cpd16507'

# Doing the comparison in matlab with the ismember between the metabolites to add
# and all the metNames in model_merged_sergio, jumped out that the next 3
# compounds are already present within the model, but not present in the sink
# because they are related just to reactions without taxonomy, therefore they
# are discarded during the extraction of the universal database.
new_phenols.loc[new_phenols['metNames'].isin(["Isosakuranetin"]),'metID'] = 'cpd03162'
new_phenols.loc[new_phenols['metNames'].isin(["Homoeriodictyol"]),'metID'] = 'cpd06648'
new_phenols.loc[new_phenols['metNames'].isin(["Prunetin"]),'metID'] = 'cpd07407'

##############################################################################
##############################################################################
##############################################################################

whole_mets_info = pd.read_excel('./output/whole_metabolites_information.xlsx', engine = "openpyxl")
a = table_results_known.Formula.apply(lambda x: re.findall('(\w{14}-\w{10}-\w{1})',x))
b = np.unique(list(itertools.chain.from_iterable(a)))
predicted_mets = pd.DataFrame(index = range(len(b)), data={'metID':"", 'metNames': b, 'metInChI': ""})
predicted_mets.metInChI = [whole_mets_info.loc[whole_mets_info['InChIKey'].isin([x]),'InChI'].values[0] for x in predicted_mets.metNames]
new_mets = new_phenols.merge(predicted_mets, how = "outer")

count = 0
for x in range(len(new_mets)):
    if not new_mets.loc[x,'metID']:
        if 22+count < 100:
            new_mets.loc[x,'metID'] = 'cpd900'+ str(22+count)
        else:
            new_mets.loc[x,'metID'] = 'cpd90'+ str(22+count)
        count += 1

names_cpd_entire = names_cpd.append(new_mets[['metNames','metID']],ignore_index=True)
# Removing the duplicates we remove the new phenols that are already present in the
# universal database with the same name, and we maintain the proper cpd
names_cpd_entire.drop_duplicates(subset = 'metNames',inplace = True)
names_cpd_entire.index = range(len(names_cpd_entire))
    
# create a dictionary to obtain the cpd formulas
dictionary_names_cpd = names_cpd_entire.set_index('metNames').T.to_dict('records')

# extract the formulas
a = [re.findall('./([\S.]+)/c',new_mets['metInChI'][x]) for x in range(len(new_mets))]
new_mets['metFormulas'] = [x[0].replace('.ClH', '') for x in a]
# remove the phenols already present in the universal database
new_phenols_to_introduce_as_mets = new_mets.drop(new_mets.loc[new_mets['metID'].isin(names_cpd.metID),].index)

new_phenols_to_exchange = new_mets.loc[new_mets['metNames'].isin(new_phenols.metNames),]
new_phenols_to_exchange.drop(new_phenols_to_exchange.loc[new_phenols_to_exchange['metID'].isin(names_cpd.metID),].index, inplace=True)
new_phenols = new_mets.loc[new_mets['metNames'].isin(new_phenols.metNames),]

# new_mets.to_excel('./output/whole_list_of_mets.xlsx', index = False)
# new_phenols.to_excel('./output/list_phenols_predicted.xlsx', index = False)
# new_phenols_to_exchange.to_excel('./output/list_phenols_to_exchange.xlsx', index = False)
# new_phenols_to_introduce_as_mets.to_excel('./output/list_phenols_to_introduce_as_mets.xlsx', index = False)

##############################################################################
##############################################################################
##############################################################################

# Extract the phenols already present in the universal database but not in AGREDA
# and extract the related reactions in order to filter those already present in
# the universal database. STOP here to do the checking manually before continuing

new_phenols_already_present = new_mets.loc[new_mets['metID'].isin(names_cpd.metID),]
prediction_phenols_already_present = table_results_known.loc[table_results_known['Name_metabolite'].isin(new_phenols_already_present.metNames),]
# prediction_phenols_already_present.to_excel('./output/reactions_for_already_present_mets_to_comment.xlsx')

##############################################################################
##############################################################################
##############################################################################

already_present = [0, 48, 116, 125, 127, 129, 190, 202] # Extracted from '../output/reactions_for_already_present_mets_to_comment.xlsx'
prediction_phenols_already_present = prediction_phenols_already_present.loc[already_present,]
table_results_known.drop(already_present, inplace=True)
table_results_known.index = range(len(table_results_known))

# generate the excel to add the new reactions to the model
new_reactions = pd.DataFrame(index= range(len(table_results_known)), columns = ['Equation','Stoichiometry','lb','ub','EC','trules','trRules', 'source'])
new_reactions['lb'] = [np.unique(str(x).split(' || '))[0] for x in table_results_known['lb']]
new_reactions['ub'] = [np.unique(str(x).split(' || '))[0] for x in table_results_known['ub']]
new_reactions['EC'] = table_results_known['EC']
new_reactions['trules'] = table_results_known['Taxonomy']
new_reactions['source'] = [dictionary_names_cpd[0][x.rstrip()] for x in table_results_known.Name_metabolite]
new_reactions['ChemicalScore'] = table_results_known['ChemicalScore']

cpd_reactions = table_results_known['Formula'].copy()

for i in range(len(cpd_reactions)):
    tmp_1 = cpd_reactions[i].split(' -> ')
    tmp_2 = tmp_1[0].split(' + ')
    tmp_2 = [x.replace('_',' ') for x in tmp_2]
    
    tmp_2 = [dictionary_names_cpd[0][x.rstrip()] for x in tmp_2]
    tmp_1[0] = ' + '.join(tmp_2)
    
    tmp_3 = tmp_1[1].split(' + ')
    tmp_3 = [x.replace('_',' ') for x in tmp_3]
    tmp_3 = [dictionary_names_cpd[0][x.rstrip()] for x in tmp_3]
    tmp_1[1] = ' + '.join(tmp_3)
    
    cpd_reactions[i] = ' -> '.join(tmp_1)    
    
new_reactions['Equation'] = cpd_reactions

# use the balance_stoichiometry package to obtain the stoichiometric values
errors_tmp = []
for i in range(len(table_results_known)):
    try:
        tmp = table_results_known['Formula_mets'][i].split(sep=' -> ')
        substrates = re.findall('\w+',tmp[0])
        products = re.findall('\w+',tmp[1])
        reac, prod = balance_stoichiometry(substrates,products,underdetermined=None)
        tmp_substrates = [str(reac[x.rstrip()]) for x in substrates]
        tmp[0] = ';'.join(tmp_substrates)
        tmp_products = [str(prod[x.rstrip()]) for x in products]
        tmp[1] = ';'.join(tmp_products)
        new_reactions['Stoichiometry'][i] = ';'.join(tmp)
    except:
        tmp = table_results_known['Formula_mets'][i].split(sep=' -> ')
        substrates = re.findall('\w+',tmp[0])
        products = re.findall('\w+',tmp[1])
        if substrates == products:
            new_reactions['Stoichiometry'][i] = '1;1'
        else:
            print(table_results_known.loc[i,'Name_metabolite'])
            errors_tmp.append(i+2)

species = pd.read_excel('./input/species.xlsx', engine = "openpyxl")
temp = ['t('+str(x)+')' for x in range(1,len(species)+1)]
species = pd.DataFrame(index = range(len(species)), data = {'trules_format':temp, 'trRules_format': species.Taxonomy})
dictionary_species = species.set_index('trules_format').T.to_dict('record')

trRules = new_reactions['trules'].copy()
for i in range(len(trRules)):
    for key in dictionary_species[0].keys():
        trRules[i] = trRules[i].replace(key,dictionary_species[0][key])
    trRules[i] = trRules[i].replace(' | ',' or ')
        
new_reactions['trRules'] = trRules

trRules = prediction_phenols_already_present['Taxonomy'].copy()
for i in trRules.index:
    for key in dictionary_species[0].keys():
        trRules[i] = trRules[i].replace(key,dictionary_species[0][key])
    trRules[i] = trRules[i].replace(' | ',' or ')
    
prediction_phenols_already_present['trRules'] = trRules.values

# new_reactions.to_excel('./output/new_phenols_reactions_CS.xlsx', index = False)
# prediction_phenols_already_present.to_excel('../files_for_model/already_present_reactions.xlsx', index = False)
'''
add manually the rxnID to the already present reactions to help their modification
and the cpdID of the mets of interest to check if they have an exchange
'''