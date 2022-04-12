# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 13:24:08 2020

@author: fbalzerani
"""

# Import packages
import json
import os
import re
import pandas as pd
import numpy as np
from rdkit import Chem

# Analyze which metabolites of the source are present in the sink (and in the metabolic network)
sink = pd.read_csv('../01_FilesGenerationAndRetroPathRLApplication/output/general_sink.csv', sep=',')
phenol_explorer = pd.read_csv('../01_FilesGenerationAndRetroPathRLApplication/output/general_source_phen_ex.csv', sep='\t',header = None, names=['Name','InChI'])

# Load the cofactors
cofactors = pd.read_csv('../01_FilesGenerationAndRetroPathRLApplication/output/cofactors_inchi.csv', sep= ',')
cofactors.loc[:,'InChIKey'] = [Chem.InchiToInchiKey(x) for x in cofactors.inchi]
cofactors.loc[cofactors.name.isin(['ADP']),'InChIKey'] = "XTWYTFMLZFPYCI-UHFFFAOYSA-N"
# Load the relation name-inchikey
name_inchi = pd.read_excel("../01_FilesGenerationAndRetroPathRLApplication/output/name_inchikey.xlsx",
                           engine = "openpyxl")
name_inchi = name_inchi.merge(cofactors[['name','inchi','InChIKey']], how = 'outer')

# Load the names-cpd relation for the metabolites in model_merged_sergio
cpd_names = pd.read_excel("../01_FilesGenerationAndRetroPathRLApplication/output/names_cpd_model_merged_sergio.xlsx",
                          engine = "openpyxl")
cpd_names.drop_duplicates(inplace=True)   
# Define the pathway to the results folder
path = '../01_FilesGenerationAndRetroPathRLApplication/09_RetroPathRL/data/results_biosensor/'

source_in_sink = []
number_total_results = 0
number_total_known_molecule = 0
whole_metabolites_information = pd.DataFrame(columns = ['Name', 'InChI', 'InChIKey', 'SMILES', 'inSink','Exc'])

for j in range(len(os.listdir(path))):
    
    # If the metabolite doesn't have any result, go to the next iteration
    if os.path.isfile(path+str(j)+"/0_results") or not [x for x in os.listdir(path+str(j)) if re.search(r'full_scope.json', x)]:
        continue
    
    # Select only the full scope file to analyze the result
    number_total_results += 1
    temp_file = [x for x in os.listdir(path+str(j)) if re.search(r'full_scope.json', x)]
    
    # Load the json of the full scope and extract the information related to nodes and edges
    with open(path+str(j)+"/"+temp_file[0], "r") as jsonFile:
        data = json.load(jsonFile)
        
    elements = data['elements']
    
    nodes = elements['nodes']
    edges = elements['edges']
    
    # Extract the information about the compounds present in the network or in the sink or as source
    id_compounds = []
    number_metabolites = 0
    metabolites_table = pd.DataFrame(columns = ['Name', 'InChI', 'InChIKey', 'SMILES', 'inSink','Exc'])
    
    for i in nodes:
        if i['data']['type']=='reaction':
            continue
        number_metabolites += 1
        if i['data']['isSource']==1 or i['data']['inSink']==1:
            id_compounds.append(i['data']['id'])
            
            tmp_name = [x for x in i['data']['Names'] if not x == i['data']['id']]
            if not tmp_name:
                if i['data']['id'] in name_inchi.InChIKey.values:
                    metabolites_table = metabolites_table.append({'Name': name_inchi.loc[name_inchi['InChIKey'].isin([i['data']['id']]),'name'].values[0], 'InChI' : i['data']['InChI'], 'InChIKey' : i['data']['id'],'SMILES' : i['data']['SMILES'], 'inSink' : 1 ,'Exc' : 0}, ignore_index=True)
                else: 
                    metabolites_table = metabolites_table.append({'Name': '', 'InChI' : i['data']['InChI'], 'InChIKey' : i['data']['id'],'SMILES' : i['data']['SMILES'], 'inSink' : 1 ,'Exc' : 0}, ignore_index=True)
            else:    
                metabolites_table = metabolites_table.append({'Name': tmp_name[0], 'InChI' : i['data']['InChI'], 'InChIKey' : i['data']['id'], 'SMILES' : i['data']['SMILES'], 'inSink' : 1 ,'Exc' : 0}, ignore_index=True)
            
        else:
            id_compounds.append(i['data']['id'])
            
            tmp_name = [x for x in i['data']['Names'] if not x == i['data']['id']]
            
            if not tmp_name:
                if i['data']['id'] in name_inchi.InChIKey.values:
                    metabolites_table = metabolites_table.append({'Name': name_inchi.loc[name_inchi['InChIKey'].isin([i['data']['id']]),'name'].values[0], 'InChI' : i['data']['InChI'], 'InChIKey' : i['data']['id'],'SMILES' : i['data']['SMILES'], 'inSink' : 0 ,'Exc' : 0}, ignore_index=True)
                else: 
                    metabolites_table = metabolites_table.append({'Name': '', 'InChI' : i['data']['InChI'], 'InChIKey' : i['data']['id'], 'SMILES' : i['data']['SMILES'], 'inSink' : 0 ,'Exc' : 0}, ignore_index=True)
            else:    
                metabolites_table = metabolites_table.append({'Name': tmp_name[0], 'InChI' : i['data']['InChI'], 'InChIKey' : i['data']['id'], 'SMILES' : i['data']['SMILES'], 'inSink' : 0 ,'Exc' : 0}, ignore_index=True)
    
    # Extract the information about the reactions that contain the compounds present in the sink or as source
    # id_compounds_rules = []
    # true_reactions = []
    number_reactions = len(nodes)-number_metabolites
    reactions_table = pd.DataFrame(columns = ['Rule_ID', 'EC_number', 'Diameter'])
    S_matrix = np.zeros((number_metabolites,number_reactions))
    index = 0
    
    for i in nodes:
        if i['data']['type']=='compound':
            continue
        tmp = i['data']['id'].split("-")
        substrate = "-".join(tmp[0:3])
        tmp_comparison = [x in id_compounds for x in i['data']['Stoechiometry']]
        products = [x for x in i['data']['Stoechiometry']]
        S_matrix[np.where([x in substrate for x in id_compounds]),index] = -1
        S_matrix[np.where([x in products for x in id_compounds]),index] = +1
        rule_id = []
        for x in i['data']['Rule ID']:
            if '_r' in x:
                rule_id.append('_'.join(x.split('_')[0:4]))
            else:
                rule_id.append('_'.join(x.split('_')[0:2]))
                               
        reactions_table = reactions_table.append({'Rule_ID': np.unique(rule_id), 'EC_number' : i['data']['EC number'], 'Diameter' : i['data']['Diameter'],'ChemicalScore':i['data']['ChemicalScore']}, ignore_index=True)
        
        index += 1
        
    # New assignation of the exchanges
    count_exc = 0
    for i in range(number_metabolites):
        if metabolites_table.loc[i,'inSink'] == 1:
            count_exc += 1
            exchange = np.zeros((number_metabolites,1))
            if 1 in S_matrix[i,:]:
                exchange[i] = -1
            elif -1 in S_matrix[i,:]:
                exchange[i] = +1
            S_matrix = np.append(S_matrix, exchange, axis=1)
            reactions_table = reactions_table.append({'Rule_ID': 'Exchange_reactions_'+str(count_exc), 'EC_number' : '', 'Diameter' : ''}, ignore_index=True)
            metabolites_table.loc[i,'Exc'] = 1
        elif metabolites_table.loc[i,'InChI'] in cofactors['inchi'].values:
            metabolites_table.loc[i,'inSink'] = 1
            count_exc += 1
            exchange = np.zeros((number_metabolites,1))
            exchange[i] = -1
            S_matrix = np.append(S_matrix, exchange, axis=1)
            reactions_table = reactions_table.append({'Rule_ID': 'Exchange_reactions_cofactor'+str(count_exc), 'EC_number' : '', 'Diameter' : ''}, ignore_index=True)
            metabolites_table.loc[i,'Exc'] = 1           
    
    if len(S_matrix) == 1:
        continue
    
    metabolites_table.Name = [x.replace('_',' ') for x in metabolites_table.Name]
    # Create a variable where store the information about all the metabolites                    
    whole_metabolites_information = whole_metabolites_information.merge(metabolites_table, how='outer')
                        
    # Save
    results_path = './output/results_to_extract_CS/'   

    if os.path.exists(results_path):
        if os.path.exists(results_path+'Metabolite_'+str(j+1)):
            np.savetxt(results_path+'Metabolite_'+str(j+1)+"/S_matrix.csv", S_matrix, delimiter=",")
            metabolites_table.to_excel(results_path+'Metabolite_'+str(j+1)+'/metabolites_table.xlsx', index=False, header = metabolites_table.columns)
            reactions_table.to_excel(results_path+'Metabolite_'+str(j+1)+'/reactions_table.xlsx', index=False, header = reactions_table.columns)
        else:
            os.mkdir(results_path+'Metabolite_'+str(j+1))
            np.savetxt(results_path+'Metabolite_'+str(j+1)+"/S_matrix.csv", S_matrix, delimiter=",")
            metabolites_table.to_excel(results_path+'Metabolite_'+str(j+1)+'/metabolites_table.xlsx', index=False, header = metabolites_table.columns)
            reactions_table.to_excel(results_path+'Metabolite_'+str(j+1)+'/reactions_table.xlsx', index=False, header = reactions_table.columns)
    else:
        os.mkdir(results_path)
        os.mkdir(results_path+'Metabolite_'+str(j+1))
        np.savetxt(results_path+'Metabolite_'+str(j+1)+"/S_matrix.csv", S_matrix, delimiter=",")
        metabolites_table.to_excel(results_path+'Metabolite_'+str(j+1)+'/metabolites_table.xlsx', index=False, header = metabolites_table.columns)
        reactions_table.to_excel(results_path+'Metabolite_'+str(j+1)+'/reactions_table.xlsx', index=False, header = reactions_table.columns)

countRetro = 0
countTemp = 0
to_del = []
ccc = 0
for i in whole_metabolites_information.index:
    if len(whole_metabolites_information.loc[whole_metabolites_information['InChIKey'].isin([whole_metabolites_information.loc[i,'InChIKey']]),]) > 1:
        tmp = [x for x in whole_metabolites_information.loc[whole_metabolites_information['InChIKey'].isin([whole_metabolites_information.loc[i,'InChIKey']]),].index if whole_metabolites_information.loc[whole_metabolites_information['InChIKey'].isin([whole_metabolites_information.loc[i,'InChIKey']]),].Name[x] == '']
        if not tmp == []:
            to_del.append([x for x in whole_metabolites_information.loc[whole_metabolites_information['InChIKey'].isin([whole_metabolites_information.loc[i,'InChIKey']]),].index if whole_metabolites_information.loc[whole_metabolites_information['InChIKey'].isin([whole_metabolites_information.loc[i,'InChIKey']]),].Name[x] == ''][0])
        
    if whole_metabolites_information.inSink[i] == 1:
        try:
            whole_metabolites_information.loc[i,'metID'] = cpd_names.loc[cpd_names['metNames'].isin([whole_metabolites_information.loc[i,'Name']]),'metID'].values[0]
        except:
            if whole_metabolites_information.loc[i,'Name'].replace('_',' ') in phenol_explorer.Name.values:
                countRetro += 1
                whole_metabolites_information.loc[i,'metID'] = 'cpdRetro'+ str(countRetro)
            else:
                whole_metabolites_information.loc[i,'metID'] = 'xxxx'
    else:
        countTemp += 1
        whole_metabolites_information.loc[i,'metID'] = 'cpdTemp'+ str(countTemp)

whole_metabolites_information.drop(to_del, inplace=True)
# whole_metabolites_information.to_excel('./output/whole_metabolites_information.xlsx', index = False)
