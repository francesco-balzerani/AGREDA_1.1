function model = addPhenols(model,phenols,mets_to_fastcore)

% addPhenols.m
% 
% Author: Francesco Balzerani
% Email: fbalzerani@tecnun.es
% Date: 19/10/2020

% Script to add a list of phenolic compounds into an AGORA-SEED fields 
% compatible model and create their corresponding exchange reactions.

n_phenols = size(phenols,1);
phenols_names = table2cell(phenols(:,'metNames'));

% Add metabolic information for all the phenols created and the predicted
% metabolites
model.metID = [model.metID; table2cell(phenols(:,'metID'))];
model.mets  = [model.mets; phenols_names];
model.metNames = [model.metNames; phenols_names];
model.metInChIString = [model.metInChIString; table2cell(phenols(:,'metInChI'))];
model.metFormulas = [model.metFormulas;table2cell(phenols(:,'metFormulas'))];

model.S = [model.S;zeros(n_phenols,length(model.rxns))];

model.b  = [model.b;zeros(n_phenols,1)];
model.csense  = [model.csense; repmat('E', n_phenols,1)];
model.metCharges = [model.metCharges;cell(n_phenols,1)];
model.metChEBIID = [model.metChEBIID;cell(n_phenols,1)];
model.metKEGGID = [model.metKEGGID;cell(n_phenols,1)];
model.metPubChemID = [model.metPubChemID;cell(n_phenols,1)];
model.metSmiles = [model.metSmiles;cell(n_phenols,1)];
model.metHMDBID = [model.metHMDBID;cell(n_phenols,1)];
model.metsTax = [model.metsTax;cell(n_phenols,1)];
model.metMass = [model.metMass;cell(n_phenols,1)];
model.metSource = [model.metSource;cell(n_phenols,1)];
model.metInchiKey = [model.metInchiKey;cell(n_phenols,1)];
model.metIsCore = [model.metIsCore;cell(n_phenols,1)];
model.metIsObsolete = [model.metIsObsolete;cell(n_phenols,1)];
model.metLinkedCompound = [model.metLinkedCompound;cell(n_phenols,1)];
model.metIsCofactor = [model.metIsCofactor;cell(n_phenols,1)];
model.metDeltaG = [model.metDeltaG;cell(n_phenols,1)];
model.metDeltaGErr = [model.metDeltaGErr;cell(n_phenols,1)];
model.metPKA = [model.metPKA;cell(n_phenols,1)];
model.metPKB = [model.metPKB;cell(n_phenols,1)];
model.metAbstractCompound = [model.metAbstractCompound;cell(n_phenols,1)];
model.metComprisedOf = [model.metComprisedOf;cell(n_phenols,1)];
model.metAliases = [model.metAliases;cell(n_phenols,1)];
model.metOntology = [model.metOntology;cell(n_phenols,1)];

% Create the exchange reactions for the new phenols, not for the predicted
% metabolites
n_mets = size(mets_to_fastcore,1);
metID_ex = table2cell(mets_to_fastcore(:,'metID'));
exchanges = find(sum(model.S ~=0,1) == 1)';

count = 0;
for i = 1 : n_mets
    pos_phen = find(ismember(model.metID, metID_ex(i,1)));
    rxn_related_met = find(model.S(pos_phen,:));
    if sum(ismember(rxn_related_met,exchanges)) == 0
        count = count+1;
        model.rxnID = [model.rxnID;['rxnAddPhenEX' num2str(i)]];
        model.S = [model.S, zeros(size(model.S,1),1)];
        model.S(pos_phen,end) = -1;
        model.rxns = [model.rxns;strcat('EX_',model.metNames(pos_phen))];
        model.rxnNames = [model.rxnNames;strcat('EX_',model.metNames(pos_phen))];
    elseif sum(ismember(rxn_related_met,exchanges)) == 0
        warning('more than one exchange')
    else
        pos_rxn = rxn_related_met(ismember(rxn_related_met,exchanges));
        model.rxnAliases{pos_rxn,1} = [model.rxnAliases{pos_rxn,1},' || old_info_in_model_merged_sergio_corresponding_to_rxnID = ', model.rxnID{pos_rxn,1}];
        model.rxnNames{pos_rxn,1} = strcat('EX_',model.metNames{pos_phen});
        model.rxns{pos_rxn,1} = strcat('EX_',model.metNames{pos_phen});
        model.rxnID{pos_rxn,1} = ['rxnAddPhenEX',num2str(i)];
        model.lb(pos_rxn,1) = -1000;
        model.ub(pos_rxn,1) = +1000;
    end
end

model.lb = [model.lb;repmat(-1000,count,1)];
model.ub = [model.ub;repmat(1000,count,1)];

model.c  = [model.c;zeros(count,1)];
model.rxnGeneMat = [model.rxnGeneMat ; sparse(count, length(model.genes))];
model.subSystems = [model.subSystems;cell(count,1)];
model.rules = [model.rules;cell(count,1)];
model.grRules = [model.grRules;cell(count,1)];
model.rxnConfidenceScores = [model.rxnConfidenceScores;zeros(count,1)];
model.rxnECNumbers = [model.rxnECNumbers;cell(count,1)];
model.rxnKEGGID = [model.rxnKEGGID;cell(count,1)];
model.correctedRxns = [model.correctedRxns;zeros(count,1)];
model.rxnCode = [model.rxnCode;cell(count,1)];
model.rxnStoichiometry = [model.rxnStoichiometry;cell(count,1)];
model.rxnIsTransport = [model.rxnIsTransport;cell(count,1)];
model.rxnEquation = [model.rxnEquation;cell(count,1)];
model.rxnDefinition = [model.rxnDefinition;cell(count,1)];
model.rxnReversibility = [model.rxnReversibility;cell(count,1)];
model.rxnDirection = [model.rxnDirection;cell(count,1)];
model.rxnAbstractReaction = [model.rxnAbstractReaction;cell(count,1)];
model.rxnPathways = [model.rxnPathways;cell(count,1)];
model.rxnAliases = [model.rxnAliases;cell(count,1)];
model.rxnDeltaG = [model.rxnDeltaG;cell(count,1)];
model.rxnDeltaGErr = [model.rxnDeltaGErr;cell(count,1)];
model.rxnCompoundIDs = [model.rxnCompoundIDs;cell(count,1)];
model.rxnStatus = [model.rxnStatus;cell(count,1)];
model.rxnIsObsolete = [model.rxnIsObsolete;cell(count,1)];
model.rxnLinkedReaction = [model.rxnLinkedReaction;cell(count,1)];
model.rxnNotes = [model.rxnNotes;cell(count,1)];
model.rxnSource = [model.rxnSource;cell(count,1)];
model.rxnOntology = [model.rxnOntology;cell(count,1)];
model.comments = [model.comments;cell(count,1)];
model.citations = [model.citations;cell(count,1)];
model.trRules = [model.trRules; cell(count,1)];
model.trules = [model.trules; cell(count,1)];
model.rxnTaxMat = [model.rxnTaxMat; sparse(count,length(model.taxonomy))];

end