function [model, cpdID_rxnID] = addPredictedReactions(model,new_reactions_table,already_present_reactions)

% addPredictedReactions.m
% 
% Author: Francesco Balzerani
% Email: fbalzerani@tecnun.es
% Date: 21/10/2020

% Extract the information regarding the new reactions and the list of all
% the phenols for which we generate new reactions

rxnEquation_file = table2cell(new_reactions_table(:,1));
Stoichiometry_file = table2cell(new_reactions_table(:,2));
EC_file = table2cell(new_reactions_table(:,5));
trules_file = table2cell(new_reactions_table(:,6));
trRules_file = table2cell(new_reactions_table(:,7));

n_rows = length(rxnEquation_file);

% Calculate rxnTaxMat

rxnTaxMat_file = sparse(length(trRules_file), length(model.taxonomy));
pos = cellfun(@regexp, trules_file, repmat({'\d*'}, n_rows, 1), repmat({'Match'}, n_rows, 1), 'UniformOutput', false);
pos = cellfun(@str2double, pos, 'UniformOutput', false);

for i = 1:n_rows
    rxnTaxMat_file(i, pos{i}) = 1;
end

% Calculate position compounds for S matrix

cpd_id_mets = cellfun(@regexp, rxnEquation_file, repmat({'cpd(\d)+'}, n_rows, 1), repmat({'Match'}, n_rows, 1), 'UniformOutput', false);
stoichiometry_values = cellfun(@strsplit, Stoichiometry_file, repmat({';'}, n_rows,1), 'UniformOutput', false);
S_file = zeros(length(model.mets),n_rows);
for i = 1 : n_rows
    clc
    disp('Creating the new reactions');
    disp(['Reaction: ' num2str(i) '/' num2str(n_rows)]);
    temp_cpd = strsplit(rxnEquation_file{i},' -> ');
    n_substrates = length(regexp(temp_cpd{1},'cpd(\d)+','Match'));
    for j = 1:length(cpd_id_mets{i})        
        if j <= n_substrates
            pos_mets = find(~cellfun(@isempty,regexp(model.metID,cpd_id_mets{i}{j})));
            if length(pos_mets) > 1
                pos_cyt = ~cellfun(@isempty,regexp(model.mets(pos_mets),'\[c\]'));
                pos_mets = pos_mets(pos_cyt);
            end
            S_file(pos_mets,i) = -str2double(stoichiometry_values{i}{j});
        else
            pos_mets = find(~cellfun(@isempty,regexp(model.metID,cpd_id_mets{i}{j})));
            if length(pos_mets) > 1
                pos_cyt = ~cellfun(@isempty,regexp(model.mets(pos_mets),'\[c\]'));
                pos_mets = pos_mets(pos_cyt);
            end
            S_file(pos_mets,i) = str2double(stoichiometry_values{i}{j});
        end
    end
end

% Define rxnNames, rxnID

rxnNames_file = cell(length(trules_file),1);
rxnID_file = cell(length(trules_file),1);
for i = 1 : n_rows
    rxnNames_file{i} = ['new_phenol_reaction_',num2str(i)];
    rxnID_file{i} = ['rxnPhenAdd',num2str(i)];
end

% Extract relation between cpdID and rxnID
cpdID_rxnID = [table2cell(new_reactions_table(:,'source')),rxnID_file];

% Fill in all the fields

model.S = [model.S,S_file];
model.rxns = [model.rxns;rxnNames_file];
model.lb = [model.lb;repmat(-1000,n_rows,1)];
model.ub = [model.ub;repmat(+1000,n_rows,1)];
model.c  = [model.c;zeros(length(rxnNames_file),1)];
model.rules = [model.rules;cell(length(rxnNames_file),1)];
model.rxnGeneMat = [model.rxnGeneMat ; sparse(length(rxnNames_file), length(model.genes))];
model.rxnNames = [model.rxnNames;rxnNames_file];
model.subSystems = [model.subSystems;cell(length(rxnNames_file),1)];
model.grRules = [model.grRules;cell(length(rxnNames_file),1)];
model.comments = [model.comments;cell(length(rxnNames_file),1)];
model.citations = [model.citations;cell(length(rxnNames_file),1)];
model.rxnConfidenceScores = [model.rxnConfidenceScores;zeros(length(rxnNames_file),1)];
model.rxnECNumbers = [model.rxnECNumbers;EC_file];
model.rxnKEGGID = [model.rxnKEGGID;cell(length(rxnNames_file),1)];
model.trRules = [model.trRules; trRules_file];
model.trules = [model.trules; trules_file];
model.rxnTaxMat = [model.rxnTaxMat; rxnTaxMat_file ];
model.correctedRxns = [model.correctedRxns;zeros(length(rxnNames_file),1)];
model.rxnID = [model.rxnID;rxnID_file];
model.rxnCode = [model.rxnCode;cell(length(rxnNames_file),1)];
model.rxnStoichiometry = [model.rxnStoichiometry;Stoichiometry_file];
model.rxnIsTransport = [model.rxnIsTransport;cell(length(rxnNames_file),1)];
model.rxnEquation = [model.rxnEquation;rxnEquation_file];
model.rxnDefinition = [model.rxnDefinition;cell(length(rxnNames_file),1)];
model.rxnReversibility = [model.rxnReversibility;cell(length(rxnNames_file),1)];
model.rxnDirection = [model.rxnDirection;cell(length(rxnNames_file),1)];
model.rxnAbstractReaction = [model.rxnAbstractReaction;cell(length(rxnNames_file),1)];
model.rxnPathways = [model.rxnPathways;cell(length(rxnNames_file),1)];
model.rxnAliases = [model.rxnAliases;cell(length(rxnNames_file),1)];
model.rxnDeltaG = [model.rxnDeltaG;cell(length(rxnNames_file),1)];
model.rxnDeltaGErr = [model.rxnDeltaGErr;cell(length(rxnNames_file),1)];
model.rxnCompoundIDs = [model.rxnCompoundIDs;cell(length(rxnNames_file),1)];
model.rxnStatus = [model.rxnStatus;cell(length(rxnNames_file),1)];
model.rxnIsObsolete = [model.rxnIsObsolete;cell(length(rxnNames_file),1)];
model.rxnLinkedReaction = [model.rxnLinkedReaction;cell(length(rxnNames_file),1)];
model.rxnNotes = [model.rxnNotes;cell(length(rxnNames_file),1)];
model.rxnSource = [model.rxnSource;cell(length(rxnNames_file),1)];
model.rxnOntology = [model.rxnOntology;cell(length(rxnNames_file),1)];

% apply changes

already_rxnID = table2cell(already_present_reactions(:,'rxnID'));
already_present_reactions = table2cell(already_present_reactions);
n_rxn_already = length(already_rxnID);

% loop to change the reactions
rxnPhenAdd = model.rxnID(find(~cellfun(@isempty, regexp(model.rxnID, 'rxnPhenAdd'))));
numbers = regexp(rxnPhenAdd,'\d+','match');
last_number = str2double(numbers{end}{1});

for k = 1 : n_rxn_already
    clc
    disp('Modifying the already present reactions');
    disp(['Reaction: ' num2str(k) '/' num2str(n_rxn_already)]);
    pos_rxn = find(ismember(model.rxnID, already_rxnID(k)));
    model.rxnAliases{pos_rxn,1} = [model.rxnAliases{pos_rxn,1},' || old_info_in_model_merged_sergio_corresponding_to_rxnID = ', model.rxnID{pos_rxn,1}];
    model.rxnNames{pos_rxn,1} = ['new_phenol_reaction_',num2str(last_number+k)];
    model.rxnID{pos_rxn,1} = ['rxnPhenAdd',num2str(last_number+k)];
    tmp_cpd_rxn{k,1} = already_present_reactions{k,11};
    tmp_cpd_rxn{k,2} = model.rxnID{pos_rxn,1};
    if isempty(model.trules{pos_rxn,1})
        model.trules{pos_rxn,1} = already_present_reactions{k,7};
        model.trRules{pos_rxn,1} = already_present_reactions{k,10};
    else
        % merge trules
        tmp = [model.trules{pos_rxn,1},' | ',already_present_reactions{k,7}];
        tmp = strjoin(unique(split(tmp,' | ')), ' | ');
        model.trules{pos_rxn,1} = tmp;
        % merge trRules
        tmp = [model.trRules{pos_rxn,1},' or ',already_present_reactions{k,10}];
        tmp = strjoin(unique(split(tmp,' or ')), ' or ');
        model.trRules{pos_rxn,1} = tmp;
    end
    if isempty(model.rxnECNumbers{pos_rxn,1})
        model.rxnECNumbers{pos_rxn,1} = already_present_reactions{k,6};
    else
        % merge EC numbers
        if ~isempty(already_present_reactions{k,6})            
            tmp = [model.rxnECNumbers{pos_rxn,1},'|',already_present_reactions{k,6}];
            tmp = strjoin(unique(split(tmp,'|')), '|');
            model.rxnECNumbers{pos_rxn,1} = tmp;
        end
    end
    model.lb(pos_rxn,1) = -1000;
    model.ub(pos_rxn,1) = +1000;
end

cpdID_rxnID = [cpdID_rxnID;tmp_cpd_rxn];