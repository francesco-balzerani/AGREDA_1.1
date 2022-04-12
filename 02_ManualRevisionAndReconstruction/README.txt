The following are the steps necessary to extract the results of RetroPath RL and reconstruct the metabolic network.
The numeration of the steps represents the related scripts within the folder.

01. Generate the S matrix and collect the information of the metabolites and reactions involved in the predicted pathway
m02. Extract all the pathways directly connected to the sink for the new phenols and the ones already present
03. Manual reaction's balance
04. Calculate the stoichiometry of the reactions, the list of metabolites to introduced in the model, the ones to apply as core in the fastcore, and the files for the reactions and metabolites already present that need to be modified
m05. Apply the whole reconstruction process and gap filling analysis.