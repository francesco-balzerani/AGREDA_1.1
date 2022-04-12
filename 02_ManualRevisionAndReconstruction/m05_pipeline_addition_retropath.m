% pipeline_addition_retropath.m 
% 
% Author: Francesco Balzerani
% Email: fbalzerani@tecnun.es
% Date: 21/10/2020

close all
clear
clc

%%%%%%%%%%%%%% file names

mets_to_add_name_file = 'list_phenols_to_introduce_as_mets';
mets_to_exchange_name_file = 'list_phenols_to_exchange';
mets_to_fastcore_name_file = 'list_phenols_predicted';
new_reactions_name_file = 'new_phenols_reactions_CS';
alredy_present_reactions_name_file = 'already_present_reactions';
phenol_explorer_name_file = 'general_source_phen_ex';
metsPredictedNotKnown_name_file = 'InformationPredictedMetabolites';
plants_name_file = 'no_bacteria';
fastFVA_path = fullfile('.', 'output', 'optionsFastcorePhenols', 'removeNotAnnotatedRxnsFastFVAPhenols');

%%%%%%%%%%%%%% model names

merg_model_name = 'agora_seed_merged_sergio';
agreda_orig_model_name = 'AGREDA_orig';
minim_model_name = 'AGREDA';
phenols_model_name = 'phenols_retropath_model';
final_model_name = 'finalModelfastcore';
definitive_model_name = 'AGREDA11';

%%%%%%%%%%%%%% load

mets_to_add = readtable(fullfile('.', 'output', [mets_to_add_name_file '.xlsx']));
mets_to_exchange = readtable(fullfile('.', 'output', [mets_to_exchange_name_file '.xlsx']));
mets_to_fastcore = readtable(fullfile('.', 'output', [mets_to_fastcore_name_file '.xlsx']));
new_reactions = readtable(fullfile('.', 'output', [new_reactions_name_file '.xlsx']));
already_present_reactions = readtable(fullfile('.', 'output', [alredy_present_reactions_name_file '.xlsx']));
metsPredictedNotKnown = readtable(fullfile('.', 'input', [metsPredictedNotKnown_name_file '.xlsx']));
phenol_explorer = readtable(fullfile('..','01_FilesGenerationAndRetroPathRLApplication','output',[phenol_explorer_name_file '.csv']),'ReadVariableNames',false);
plants = readtable(fullfile('.', 'input',[plants_name_file '.xlsx']));

%%%%%%%%%%%%%% define the options for the FastcoreWeighted

options.modelName = phenols_model_name;
options.fastFVAName = fullfile('.', 'output', 'optionsFastcorePhenols', ['fastFVAPhenols_', phenols_model_name]);
options.saving = 1;
if exist(fullfile('.','output', 'optionsFastcorePhenols', ['blockedPhenolsRxnsUp_', phenols_model_name, '.mat']),'file')
    load(fullfile('.','output', 'optionsFastcorePhenols', ['blockedPhenolsRxnsUp_', phenols_model_name , '.mat']))
    options.blockedPhenolsRxnsUp = blockedPhenolsRxnsUp;
end
if exist(fullfile('.','output', 'optionsFastcorePhenols', ['blockedPhenolsRxnsOut_', phenols_model_name, '.mat']),'file')
    load(fullfile('.','output', 'optionsFastcorePhenols', ['blockedPhenolsRxnsOut_', phenols_model_name, '.mat']))
    options.blockedPhenolsRxnsOut = blockedPhenolsRxnsOut;
end
options.plantRxns = table2cell(plants(:,1));

%% Add Paths

addpath(genpath(fullfile('..','..','AGREDA1_1')));

%% Add the metabolites and reactions and connect them to AGREDA
clc
disp('Loading the Universal Database')
load(fullfile('..','01_FilesGenerationAndRetroPathRLApplication','input', [merg_model_name '.mat']));

clc
disp('Loading the AGREDA model')
load(fullfile('.', 'input', [agreda_orig_model_name '.mat']));
agreda = model;

% Change the ID of some exchange because they don't correspond between
% model_merged_sergio and the final version of agreda, maybe due to manual
% curation. The change regards just the IDs, nothing structural

agreda_addEX = [agreda.rxnID(~cellfun(@isempty,regexp(agreda.rxnID,'rxnAddEX'))),agreda.rxns(~cellfun(@isempty,regexp(agreda.rxnID,'rxnAddEX')))];
model_merged_addEX = [model_merged_sergio.rxnID(~cellfun(@isempty,regexp(model_merged_sergio.rxnID,'rxnAddEX'))),model_merged_sergio.rxns(~cellfun(@isempty,regexp(model_merged_sergio.rxnID,'rxnAddEX')))];
n_ids = length(agreda_addEX);

for i = 1 : n_ids
    pos = find(ismember(model_merged_addEX(:,2),agreda_addEX(i,2)));
    agreda.rxnID(ismember(agreda.rxnID,agreda_addEX(i,1))) = model_merged_sergio.rxnID(ismember(model_merged_sergio.rxnID,model_merged_addEX(pos,1))); 
end

agreda_addNewEX = [agreda.rxnID(~cellfun(@isempty,regexp(agreda.rxnID,'rxnAddNewEX'))),agreda.rxns(~cellfun(@isempty,regexp(agreda.rxnID,'rxnAddNewEX')))];
model_merged_addNewEX = [model_merged_sergio.rxnID(~cellfun(@isempty,regexp(model_merged_sergio.rxnID,'rxnAddNewEX'))),model_merged_sergio.rxns(~cellfun(@isempty,regexp(model_merged_sergio.rxnID,'rxnAddNewEX')))];
n_ids = length(agreda_addNewEX);

for i = 1 : n_ids
    pos = find(ismember(model_merged_addNewEX(:,2),agreda_addNewEX(i,2)));
    agreda.rxnID(ismember(agreda.rxnNames,agreda_addNewEX(i,2))) = model_merged_sergio.rxnID(ismember(model_merged_sergio.rxnID,model_merged_addNewEX(pos,1))); 
end

% save(fullfile('.', 'output', 'outputReconstruction', [minim_model_name '.mat']), 'agreda', '-v7');

% Add information about phenols to the network
model_with_phenols = addPhenols(model_merged_sergio, mets_to_add, mets_to_fastcore);

% Add the predicted reactions for the phenols
[model_with_retropath, cpdID_rxnID] = addPredictedReactions(model_with_phenols, new_reactions, already_present_reactions);

model_with_retropath_corrected = removeNotKnownMets(model_with_retropath, metsPredictedNotKnown);
% save(fullfile('.', 'output', 'outputReconstruction', 'model_with_mets_and_rxns_retropath.mat'), 'model_with_retropath_corrected', '-v7');

initCobraToolbox(false)

core = find(ismember(model_merged_sergio.rxnID,agreda.rxnID));
sense = 'both';

[finalModel, PhenolsAddedRxns] = applyFastcorePhenols(model_with_retropath_corrected, core, cpdID_rxnID, sense, options);
% save(fullfile('.', 'output', 'outputReconstruction', [final_model_name '.mat']), 'finalModel', '-v7');
% save(fullfile('.', 'output', 'outputReconstruction', 'addedrxns.mat'), 'PhenolsAddedRxns');

% Add manually a reaction (rxn20708) for the 2,4-Dihydroxybenzoic acid (cpd90073),
% because it's not selected by the fastcore, but has the information we
% need for the species analysis

reaction_to_add = ismember(model_with_retropath_corrected.rxnID,'rxn20708');

comp = find(model_with_retropath_corrected.S(:,reaction_to_add));

tmp_S = zeros(length(finalModel.metNames),1);

for i = 1 : length(comp)
    
    pos_comp = find(ismember(finalModel.metID,model_with_retropath_corrected.metID(comp(i))));
    
    if length(pos_comp) > 1 
        pos_cyt = ~cellfun(@isempty,regexp(finalModel.mets(pos_comp),'\[c\]'));
        pos_comp = pos_comp(pos_cyt);
    end
    
    tmp_S(pos_comp) = model_with_retropath_corrected.S(comp(i),reaction_to_add);    
end

finalModel.S = [finalModel.S,tmp_S];
finalModel.rxns = [finalModel.rxns;model_with_retropath_corrected.rxns(reaction_to_add)];
finalModel.lb = [finalModel.lb;model_with_retropath_corrected.lb(reaction_to_add)];
finalModel.ub = [finalModel.ub;model_with_retropath_corrected.ub(reaction_to_add)];
finalModel.c  = [finalModel.c;model_with_retropath_corrected.c(reaction_to_add)];
finalModel.rules = [finalModel.rules;model_with_retropath_corrected.rules(reaction_to_add)];
finalModel.rxnGeneMat = [finalModel.rxnGeneMat ; model_with_retropath_corrected.rxnGeneMat(reaction_to_add,:)];
finalModel.rxnNames = [finalModel.rxnNames;model_with_retropath_corrected.rxnNames(reaction_to_add)];
finalModel.subSystems = [finalModel.subSystems;model_with_retropath_corrected.subSystems(reaction_to_add)];
finalModel.grRules = [finalModel.grRules;model_with_retropath_corrected.grRules(reaction_to_add)];
finalModel.comments = [finalModel.comments;model_with_retropath_corrected.comments(reaction_to_add)];
finalModel.citations = [finalModel.citations;model_with_retropath_corrected.citations(reaction_to_add)];
finalModel.rxnConfidenceScores = [finalModel.rxnConfidenceScores;model_with_retropath_corrected.rxnConfidenceScores(reaction_to_add)];
finalModel.rxnECNumbers = [finalModel.rxnECNumbers;model_with_retropath_corrected.rxnECNumbers(reaction_to_add)];
finalModel.rxnKEGGID = [finalModel.rxnKEGGID;model_with_retropath_corrected.rxnKEGGID(reaction_to_add)];
finalModel.trRules = [finalModel.trRules; model_with_retropath_corrected.trRules(reaction_to_add)];
finalModel.trules = [finalModel.trules; model_with_retropath_corrected.trules(reaction_to_add)];
finalModel.rxnTaxMat = [finalModel.rxnTaxMat; model_with_retropath_corrected.rxnTaxMat(reaction_to_add,:) ];
finalModel.correctedRxns = [finalModel.correctedRxns;model_with_retropath_corrected.correctedRxns(reaction_to_add)];
finalModel.rxnID = [finalModel.rxnID;model_with_retropath_corrected.rxnID(reaction_to_add)];
finalModel.rxnCode = [finalModel.rxnCode;model_with_retropath_corrected.rxnCode(reaction_to_add)];
finalModel.rxnStoichiometry = [finalModel.rxnStoichiometry;model_with_retropath_corrected.rxnStoichiometry(reaction_to_add)];
finalModel.rxnIsTransport = [finalModel.rxnIsTransport;model_with_retropath_corrected.rxnIsTransport(reaction_to_add)];
finalModel.rxnEquation = [finalModel.rxnEquation;model_with_retropath_corrected.rxnEquation(reaction_to_add)];
finalModel.rxnDefinition = [finalModel.rxnDefinition;model_with_retropath_corrected.rxnDefinition(reaction_to_add)];
finalModel.rxnReversibility = [finalModel.rxnReversibility;model_with_retropath_corrected.rxnReversibility(reaction_to_add)];
finalModel.rxnDirection = [finalModel.rxnDirection;model_with_retropath_corrected.rxnDirection(reaction_to_add)];
finalModel.rxnAbstractReaction = [finalModel.rxnAbstractReaction;model_with_retropath_corrected.rxnAbstractReaction(reaction_to_add)];
finalModel.rxnPathways = [finalModel.rxnPathways;model_with_retropath_corrected.rxnPathways(reaction_to_add)];
finalModel.rxnAliases = [finalModel.rxnAliases;model_with_retropath_corrected.rxnAliases(reaction_to_add)];
finalModel.rxnDeltaG = [finalModel.rxnDeltaG;model_with_retropath_corrected.rxnDeltaG(reaction_to_add)];
finalModel.rxnDeltaGErr = [finalModel.rxnDeltaGErr;model_with_retropath_corrected.rxnDeltaGErr(reaction_to_add)];
finalModel.rxnCompoundIDs = [finalModel.rxnCompoundIDs;model_with_retropath_corrected.rxnCompoundIDs(reaction_to_add)];
finalModel.rxnStatus = [finalModel.rxnStatus;model_with_retropath_corrected.rxnStatus(reaction_to_add)];
finalModel.rxnIsObsolete = [finalModel.rxnIsObsolete;model_with_retropath_corrected.rxnIsObsolete(reaction_to_add)];
finalModel.rxnLinkedReaction = [finalModel.rxnLinkedReaction;model_with_retropath_corrected.rxnLinkedReaction(reaction_to_add)];
finalModel.rxnNotes = [finalModel.rxnNotes;model_with_retropath_corrected.rxnNotes(reaction_to_add)];
finalModel.rxnSource = [finalModel.rxnSource;model_with_retropath_corrected.rxnSource(reaction_to_add)];
finalModel.rxnOntology = [finalModel.rxnOntology;model_with_retropath_corrected.rxnOntology(reaction_to_add)];

% Add the transport reactions already present in AGREDA    
    
transport_rxns = find(~cellfun(@isempty,regexp(agreda.rxnID, 'rxnTransport')));
n_transport = length(transport_rxns);

for i = 1 : n_transport
    tmp_met = agreda.metID(find(agreda.S(:,transport_rxns(i))));
    if isequal(tmp_met,{'cpd01202'})
        tmp_pos = find(ismember(finalModel.metNames,'Phloretate'));
    else
        tmp_pos = find(ismember(finalModel.metID,tmp_met));
    end
    
    tmp_S = zeros(length(finalModel.metNames),1);
    tmp_S(tmp_pos) = 1;
    
    finalModel.S = [finalModel.S,tmp_S];
    finalModel.rxns = [finalModel.rxns;agreda.rxns(transport_rxns(i))];
    finalModel.lb = [finalModel.lb;agreda.lb(transport_rxns(i))];
    finalModel.ub = [finalModel.ub;agreda.ub(transport_rxns(i))];
    finalModel.c  = [finalModel.c;agreda.c(transport_rxns(i))];
    finalModel.rules = [finalModel.rules;agreda.rules(transport_rxns(i))];
    finalModel.rxnGeneMat = [finalModel.rxnGeneMat ; agreda.rxnGeneMat(transport_rxns(i),:)];
    finalModel.rxnNames = [finalModel.rxnNames;agreda.rxnNames(transport_rxns(i))];
    finalModel.subSystems = [finalModel.subSystems;agreda.subSystems(transport_rxns(i))];
    finalModel.grRules = [finalModel.grRules;agreda.grRules(transport_rxns(i))];
    finalModel.comments = [finalModel.comments;agreda.comments(transport_rxns(i))];
    finalModel.citations = [finalModel.citations;agreda.citations(transport_rxns(i))];
    finalModel.rxnConfidenceScores = [finalModel.rxnConfidenceScores;agreda.rxnConfidenceScores(transport_rxns(i))];
    finalModel.rxnECNumbers = [finalModel.rxnECNumbers;agreda.rxnECNumbers(transport_rxns(i))];
    finalModel.rxnKEGGID = [finalModel.rxnKEGGID;agreda.rxnKEGGID(transport_rxns(i))];
    finalModel.trRules = [finalModel.trRules; agreda.trRules(transport_rxns(i))];
    finalModel.trules = [finalModel.trules; agreda.trules(transport_rxns(i))];
    finalModel.rxnTaxMat = [finalModel.rxnTaxMat; agreda.rxnTaxMat(transport_rxns(i),:) ];
    finalModel.correctedRxns = [finalModel.correctedRxns;agreda.correctedRxns(transport_rxns(i))];
    finalModel.rxnID = [finalModel.rxnID;agreda.rxnID(transport_rxns(i))];
    finalModel.rxnCode = [finalModel.rxnCode;agreda.rxnCode(transport_rxns(i))];
    finalModel.rxnStoichiometry = [finalModel.rxnStoichiometry;agreda.rxnStoichiometry(transport_rxns(i))];
    finalModel.rxnIsTransport = [finalModel.rxnIsTransport;agreda.rxnIsTransport(transport_rxns(i))];
    finalModel.rxnEquation = [finalModel.rxnEquation;agreda.rxnEquation(transport_rxns(i))];
    finalModel.rxnDefinition = [finalModel.rxnDefinition;agreda.rxnDefinition(transport_rxns(i))];
    finalModel.rxnReversibility = [finalModel.rxnReversibility;agreda.rxnReversibility(transport_rxns(i))];
    finalModel.rxnDirection = [finalModel.rxnDirection;agreda.rxnDirection(transport_rxns(i))];
    finalModel.rxnAbstractReaction = [finalModel.rxnAbstractReaction;agreda.rxnAbstractReaction(transport_rxns(i))];
    finalModel.rxnPathways = [finalModel.rxnPathways;agreda.rxnPathways(transport_rxns(i))];
    finalModel.rxnAliases = [finalModel.rxnAliases;agreda.rxnAliases(transport_rxns(i))];
    finalModel.rxnDeltaG = [finalModel.rxnDeltaG;agreda.rxnDeltaG(transport_rxns(i))];
    finalModel.rxnDeltaGErr = [finalModel.rxnDeltaGErr;agreda.rxnDeltaGErr(transport_rxns(i))];
    finalModel.rxnCompoundIDs = [finalModel.rxnCompoundIDs;agreda.rxnCompoundIDs(transport_rxns(i))];
    finalModel.rxnStatus = [finalModel.rxnStatus;agreda.rxnStatus(transport_rxns(i))];
    finalModel.rxnIsObsolete = [finalModel.rxnIsObsolete;agreda.rxnIsObsolete(transport_rxns(i))];
    finalModel.rxnLinkedReaction = [finalModel.rxnLinkedReaction;agreda.rxnLinkedReaction(transport_rxns(i))];
    finalModel.rxnNotes = [finalModel.rxnNotes;agreda.rxnNotes(transport_rxns(i))];
    finalModel.rxnSource = [finalModel.rxnSource;agreda.rxnSource(transport_rxns(i))];
    finalModel.rxnOntology = [finalModel.rxnOntology;agreda.rxnOntology(transport_rxns(i))];
end

model = removeNotAnnotatedRxnsPhenols(finalModel, fastFVA_path);

% Merge the tRules of AGREDA with those that derive from
% model_merged_sergio
n_rxnid = length(agreda.rxnID);

for i = 1 : n_rxnid
    tmp_pos_model = find(ismember(model.rxnID,agreda.rxnID(i)));
    if isempty(model.trRules{tmp_pos_model})
        continue
    end
    a = strsplit(agreda.trRules{i},' or ');
    b = strsplit(model.trRules{tmp_pos_model},' or ');
    c = unique([a,b]);
    model.trRules{tmp_pos_model} = strjoin(c, ' or ');
end

rxnsWithoutTax = find(cellfun(@isempty, model.trules));
for i = 1:length(rxnsWithoutTax)
    met = find(model.S(:,rxnsWithoutTax(i)));
    rxns = find(model.S(met,:));
    rxns = rxns(~ismember(rxns,rxnsWithoutTax(i)));
    if length(rxns) > 1
        tmp = cellfun(@strsplit, model.trules(rxns), repmat({' | '}, length(rxns), 1), ...
            'UniformOutput', false);
        tmp = cellfun(@(x) x', tmp, 'UniformOutput', false);
        model.trules{rxnsWithoutTax(i)} = strjoin(unique(cat(1,tmp{:})),' | ');
        tmp = cellfun(@strsplit, model.trRules(rxns), repmat({' or '}, length(rxns), 1), ...
            'UniformOutput', false);
        tmp = cellfun(@(x) x', tmp, 'UniformOutput', false);
        model.trRules{rxnsWithoutTax(i)} = strjoin(unique(cat(1,tmp{:})),' or ');
    else
        model.trules(rxnsWithoutTax(i)) = model.trules(rxns);
        model.trRules(rxnsWithoutTax(i)) = model.trRules(rxns);
    end
end

mergeMetsSEED = readtable(fullfile('.','input','mergeMetabolitesSEED.xlsx'));
removeMets = readtable(fullfile('.','input','removeMetabolitesSEED.xlsx'));

model = removeProblematicMets(model,table2cell(removeMets(:,'ID')));
model = mergeMets(model,table2cell(mergeMetsSEED(:,'manID')),table2cell(mergeMetsSEED(:,'remID')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Level out the bounds of the AGREDA's reactions in the new reconstruction
% to the value in AGREDA
for i = 1 : n_rxnid
    tmp_pos_model = find(ismember(model.rxnID,agreda.rxnID(i)));
    
    if model.lb(tmp_pos_model) ~= agreda.lb(i)
        disp([num2str(i),': ',num2str(model.lb(tmp_pos_model)),' / ',num2str(agreda.lb(i))])
        model.lb(tmp_pos_model) = agreda.lb(i);
        model.ub(tmp_pos_model) = agreda.ub(i);
    end
end

% Add the metFormulas to all the AGREDA's mets

n_mets = length(agreda.metID);

for i = 1 : n_mets
    tmp_pos_model = find(ismember(model.metID,agreda.metID(i)));
    model.metFormulas(tmp_pos_model,1) = agreda.metFormulas(i);
end

% save(fullfile('.', 'output', 'outputReconstruction', 'model_before_single_species.mat'), 'model', '-v7');

path_new = fullfile('.','output','new_fva');
model_new = balanceSpecies(model,path_new);

% save(fullfile('.', 'output', 'outputReconstruction', 'model_after_single_species.mat'), 'model_new', '-v7');

%%%%%%%%%%%%%%%%%%%%%%% remove reactions related to predicted mets that 
%%%%%%%%%%%%%%%%%%%%%%% generate cycles, applyting fastFVA when they are
%%%%%%%%%%%%%%%%%%%%%%% irreversible

predicted_mets_info = table2cell(metsPredictedNotKnown(:,[1,14]));
rxns_to_change = [];

for i = 1 : length(predicted_mets_info)
    if predicted_mets_info{i,2} == 1
        tmp_pos_predicted = find(ismember(model_new.metID,predicted_mets_info(i,1)));
        related_rxns = find(model_new.S(tmp_pos_predicted,:));
        rxns_to_change = [rxns_to_change, related_rxns];
    end
end

model_new_tmp = model_new;
model_new_tmp.lb(rxns_to_change) = 0;

[vmin_new, vmax_new] = fastFVA(model_new_tmp,0);
epsilon = 1e-06;
rxns_to_delete = find(abs(vmin_new)<=epsilon & vmax_new<=epsilon);

model_new.S(:,rxns_to_delete) = [];
model_new.rxns(rxns_to_delete) = [];
model_new.lb(rxns_to_delete) = [];
model_new.ub(rxns_to_delete) = [];
model_new.c(rxns_to_delete) = [];
model_new.rules(rxns_to_delete) = [];
model_new.rxnNames(rxns_to_delete) = [];
model_new.subSystems(rxns_to_delete) = [];
model_new.grRules(rxns_to_delete) = [];
model_new.comments(rxns_to_delete) = [];
model_new.citations(rxns_to_delete) = [];
model_new.rxnConfidenceScores(rxns_to_delete) = [];
model_new.correctedRxns(rxns_to_delete) = [];
model_new.rxnECNumbers(rxns_to_delete) = [];
model_new.rxnKEGGID(rxns_to_delete) = [];
model_new.rxnID(rxns_to_delete) = [];
model_new.rxnCode(rxns_to_delete) = [];
model_new.rxnStoichiometry(rxns_to_delete) = [];
model_new.rxnIsTransport(rxns_to_delete) = [];
model_new.rxnEquation(rxns_to_delete) = [];
model_new.rxnDefinition(rxns_to_delete) = [];
model_new.rxnReversibility(rxns_to_delete) = [];
model_new.rxnDirection(rxns_to_delete) = [];
model_new.rxnAbstractReaction(rxns_to_delete) = [];
model_new.rxnPathways(rxns_to_delete) = [];
model_new.rxnAliases(rxns_to_delete) = [];
model_new.rxnDeltaG(rxns_to_delete) = [];
model_new.rxnDeltaGErr(rxns_to_delete) = [];
model_new.rxnCompoundIDs(rxns_to_delete) = [];
model_new.rxnStatus(rxns_to_delete) = [];
model_new.rxnIsObsolete(rxns_to_delete) = [];
model_new.rxnLinkedReaction(rxns_to_delete) = [];
model_new.rxnNotes(rxns_to_delete) = [];
model_new.rxnSource(rxns_to_delete) = [];
model_new.rxnOntology(rxns_to_delete) = [];
model_new.trules(rxns_to_delete) = [];
model_new.trRules(rxns_to_delete) = [];
model_new.rxnGeneMat(rxns_to_delete,:) = [];
model_new.rxnTaxMat(rxns_to_delete,:) = [];

metsToDelete = find(sum(model_new.S~=0, 2)==0);

model_new.S(metsToDelete,:) = [];
model_new.mets(metsToDelete) = [];
model_new.b(metsToDelete) = [];
model_new.csense(metsToDelete) = [];
model_new.metNames(metsToDelete) = [];
model_new.metCharges(metsToDelete) = [];
model_new.metFormulas(metsToDelete) = [];
model_new.metChEBIID(metsToDelete) = [];
model_new.metKEGGID(metsToDelete) = [];
model_new.metPubChemID(metsToDelete) = [];
model_new.metInChIString(metsToDelete) = [];
model_new.metHMDBID(metsToDelete) = [];
model_new.metSmiles(metsToDelete) = [];
model_new.metID(metsToDelete) = [];
model_new.metsTax(metsToDelete) = [];
model_new.metMass(metsToDelete) = [];
model_new.metSource(metsToDelete) = [];
model_new.metInchiKey(metsToDelete) = [];
model_new.metIsCore(metsToDelete) = [];
model_new.metIsObsolete(metsToDelete) = [];
model_new.metLinkedCompound(metsToDelete) = [];
model_new.metIsCofactor(metsToDelete) = [];
model_new.metDeltaG(metsToDelete) = [];
model_new.metDeltaGErr(metsToDelete) = [];
model_new.metPKA(metsToDelete) = [];
model_new.metPKB(metsToDelete) = [];
model_new.metAbstractCompound(metsToDelete) = [];
model_new.metComprisedOf(metsToDelete) = [];
model_new.metAliases(metsToDelete) = [];
model_new.metOntology(metsToDelete) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modify the Mirtillin reaction

id_rxn = 'rxnAdd18';

pos_rxn = find(ismember(model_new.rxnID,id_rxn));

substrates = {'cpd08920';'cpd00001'}; % Mirtillin + Water
stoic_sub = [-1;-4];

products = {'cpd00027';'cpd01020';'cpd01473';'cpd00067';'cpd01024';'cpd00011'}; % D-glucose + Gallate + Phloroglucinol + proton + Methane + carbon dioxide
stoic_prod = [1;1;1;1;1;1];

stoic = {'1;4;1;1;1;1;1;1'};

model_new.S(:,pos_rxn) = 0;

for i = 1 : length(stoic_sub)
    tmp_met = find(ismember(model_new.metID,substrates(i)));
    if length(tmp_met) > 1
        tmp_met = tmp_met(~cellfun(@isempty,regexp(model_new.mets(tmp_met),'\[c\]')));
    end
    model_new.S(tmp_met,pos_rxn) = stoic_sub(i);
end

for i = 1 : length(stoic_prod)
    tmp_met = find(ismember(model_new.metID,products(i)));
    if length(tmp_met) > 1
        tmp_met = tmp_met(~cellfun(@isempty,regexp(model_new.mets(tmp_met),'\[c\]')));
    end
    model_new.S(tmp_met,pos_rxn) = stoic_prod(i);
end

model_new.rxnStoichiometry(pos_rxn) = stoic;

model = model_new;

% Modify 18 reactions accordingly to AGREDA 1.0
n = length(agreda.rxnID);
count = 0;
not_agreda10 = [];
not_agreda11 = [];
count_bad = 0;
for i = 1 : n
    pos_rxn = find(ismember(model.rxnID,agreda.rxnID(i)));
    supra_compounds = model.mets(find(model.S(:,pos_rxn)));
    agreda10_compounds = agreda.mets(find(agreda.S(:,i)));
    if isequal(supra_compounds,agreda10_compounds) || isequal(model.rxnID(pos_rxn),{'rxnAdd18'})
        count = count + 1;
    else
        not_agreda10 = [not_agreda10;i];
        not_agreda11 = [not_agreda11;pos_rxn];
        model.S(:,pos_rxn) = zeros(length(model.metID),1);
        for j = 1 : length(agreda10_compounds)
            agreda_value = full(agreda.S(find(ismember(agreda.mets,agreda10_compounds(j))),i));
            agreda_value = agreda_value(agreda_value~=0);
            if length(find(ismember(model.mets,agreda10_compounds(j)))) > 1
                count_bad = count +1 ;
            end
            model.S(ismember(model.mets,agreda10_compounds(j)),pos_rxn) = agreda_value;
        end
    end
end

% Remove 2 reactions and change directionality of 1 reaction after manual
% checking of the predicted reactions

to_remove_id = [{'rxnPhenAdd195'};{'rxnPhenAdd292'}];
rxns_to_delete = find(ismember(model.rxnID,to_remove_id));

model.S(:,rxns_to_delete) = [];
model.rxns(rxns_to_delete) = [];
model.lb(rxns_to_delete) = [];
model.ub(rxns_to_delete) = [];
model.c(rxns_to_delete) = [];
model.rules(rxns_to_delete) = [];
model.rxnNames(rxns_to_delete) = [];
model.subSystems(rxns_to_delete) = [];
model.grRules(rxns_to_delete) = [];
model.comments(rxns_to_delete) = [];
model.citations(rxns_to_delete) = [];
model.rxnConfidenceScores(rxns_to_delete) = [];
model.correctedRxns(rxns_to_delete) = [];
model.rxnECNumbers(rxns_to_delete) = [];
model.rxnKEGGID(rxns_to_delete) = [];
model.rxnID(rxns_to_delete) = [];
model.rxnCode(rxns_to_delete) = [];
model.rxnStoichiometry(rxns_to_delete) = [];
model.rxnIsTransport(rxns_to_delete) = [];
model.rxnEquation(rxns_to_delete) = [];
model.rxnDefinition(rxns_to_delete) = [];
model.rxnReversibility(rxns_to_delete) = [];
model.rxnDirection(rxns_to_delete) = [];
model.rxnAbstractReaction(rxns_to_delete) = [];
model.rxnPathways(rxns_to_delete) = [];
model.rxnAliases(rxns_to_delete) = [];
model.rxnDeltaG(rxns_to_delete) = [];
model.rxnDeltaGErr(rxns_to_delete) = [];
model.rxnCompoundIDs(rxns_to_delete) = [];
model.rxnStatus(rxns_to_delete) = [];
model.rxnIsObsolete(rxns_to_delete) = [];
model.rxnLinkedReaction(rxns_to_delete) = [];
model.rxnNotes(rxns_to_delete) = [];
model.rxnSource(rxns_to_delete) = [];
model.rxnOntology(rxns_to_delete) = [];
model.trules(rxns_to_delete) = [];
model.trRules(rxns_to_delete) = [];
model.rxnGeneMat(rxns_to_delete,:) = [];
model.rxnTaxMat(rxns_to_delete,:) = [];

to_modify = {'rxnPhenAdd192'};
model.ub(ismember(model.rxnID,to_modify)) = 0;

[vmin, vmax] = fastFVA(model,0);
epsilon = 1e-06;
rxns_to_delete = find(abs(vmin)<=epsilon & vmax<=epsilon);

model.S(:,rxns_to_delete) = [];
model.rxns(rxns_to_delete) = [];
model.lb(rxns_to_delete) = [];
model.ub(rxns_to_delete) = [];
model.c(rxns_to_delete) = [];
model.rules(rxns_to_delete) = [];
model.rxnNames(rxns_to_delete) = [];
model.subSystems(rxns_to_delete) = [];
model.grRules(rxns_to_delete) = [];
model.comments(rxns_to_delete) = [];
model.citations(rxns_to_delete) = [];
model.rxnConfidenceScores(rxns_to_delete) = [];
model.correctedRxns(rxns_to_delete) = [];
model.rxnECNumbers(rxns_to_delete) = [];
model.rxnKEGGID(rxns_to_delete) = [];
model.rxnID(rxns_to_delete) = [];
model.rxnCode(rxns_to_delete) = [];
model.rxnStoichiometry(rxns_to_delete) = [];
model.rxnIsTransport(rxns_to_delete) = [];
model.rxnEquation(rxns_to_delete) = [];
model.rxnDefinition(rxns_to_delete) = [];
model.rxnReversibility(rxns_to_delete) = [];
model.rxnDirection(rxns_to_delete) = [];
model.rxnAbstractReaction(rxns_to_delete) = [];
model.rxnPathways(rxns_to_delete) = [];
model.rxnAliases(rxns_to_delete) = [];
model.rxnDeltaG(rxns_to_delete) = [];
model.rxnDeltaGErr(rxns_to_delete) = [];
model.rxnCompoundIDs(rxns_to_delete) = [];
model.rxnStatus(rxns_to_delete) = [];
model.rxnIsObsolete(rxns_to_delete) = [];
model.rxnLinkedReaction(rxns_to_delete) = [];
model.rxnNotes(rxns_to_delete) = [];
model.rxnSource(rxns_to_delete) = [];
model.rxnOntology(rxns_to_delete) = [];
model.trules(rxns_to_delete) = [];
model.trRules(rxns_to_delete) = [];
model.rxnGeneMat(rxns_to_delete,:) = [];
model.rxnTaxMat(rxns_to_delete,:) = [];

metsToDelete = find(sum(model.S~=0, 2)==0);

model.S(metsToDelete,:) = [];
model.mets(metsToDelete) = [];
model.b(metsToDelete) = [];
model.csense(metsToDelete) = [];
model.metNames(metsToDelete) = [];
model.metCharges(metsToDelete) = [];
model.metFormulas(metsToDelete) = [];
model.metChEBIID(metsToDelete) = [];
model.metKEGGID(metsToDelete) = [];
model.metPubChemID(metsToDelete) = [];
model.metInChIString(metsToDelete) = [];
model.metHMDBID(metsToDelete) = [];
model.metSmiles(metsToDelete) = [];
model.metID(metsToDelete) = [];
model.metsTax(metsToDelete) = [];
model.metMass(metsToDelete) = [];
model.metSource(metsToDelete) = [];
model.metInchiKey(metsToDelete) = [];
model.metIsCore(metsToDelete) = [];
model.metIsObsolete(metsToDelete) = [];
model.metLinkedCompound(metsToDelete) = [];
model.metIsCofactor(metsToDelete) = [];
model.metDeltaG(metsToDelete) = [];
model.metDeltaGErr(metsToDelete) = [];
model.metPKA(metsToDelete) = [];
model.metPKB(metsToDelete) = [];
model.metAbstractCompound(metsToDelete) = [];
model.metComprisedOf(metsToDelete) = [];
model.metAliases(metsToDelete) = [];
model.metOntology(metsToDelete) = [];

% save(fullfile('.', 'output', 'outputReconstruction', [definitive_model_name '.mat']), 'model', '-v7');

% extract the list of new phenols present in the final reconstruction
[r,~] = find(model.S(:,~cellfun(@isempty,regexp(model.rxnID,'rxnAddPhenEX'))));
T = table(model.metNames(r),'VariableNames',{'metNames'});
% writetable(T,'./output/outputReconstruction/list_phenols_added_met_net.xlsx') 
% modify manually the name of those that were already present in the UB:
% 6-Hydroxydaidzein -> 6,7,4'-Trihydroxyisoflavone
% 4',5,6,7-tetrahydroxyisoflavone -> 5,6,7,4'-Tetrahydroxyisoflavone
% Cyanidin 3-O-3'',6''-O-dimalonylglucoside -> Cyanidin 3-O-(3'',6''-O-dimalonyl-glucoside)
% Protocatechualdehyde -> Protocatechuic aldehyde

% Extract new reactions-species and new EX-species

t1 = [model.rxnID(~cellfun(@isempty,regexp(model.rxnID,'rxnPhenAdd'))), model.trRules(~cellfun(@isempty,regexp(model.rxnID,'rxnPhenAdd')))];
T = table(t1(:,1), t1(:,2),'VariableNames',{'ID','Taxonomy'});
% writetable(T,'./output/outputReconstruction/new_reactions_species.xlsx') 