The following are the steps necessary to generate all the files to run the RetroPath RL algorithm.
The numeration of the steps represents the related scripts within the folder.

m01. Extract the stoichiometry matrix, the list of reactions and metabolites of the universal database. Extract the information about the reactions from the literature.
02. Create the sink file
03. Create the source file
04. Generate the SMILES version of the reaction from the literature
05. Generate SMARTS rules manually
06. Generate a files with all the rules of the reactions from the literature
07. Generate the rules files
08. Calculate the similarity between each source and each substrate of the rules, to extract the cut off value
09. Apply the RetroPath RL algorithm

To apply this pipeline, download the RetroRules database (https://retrorules.org/dl/retrorules_dump) and extract in the input folder.