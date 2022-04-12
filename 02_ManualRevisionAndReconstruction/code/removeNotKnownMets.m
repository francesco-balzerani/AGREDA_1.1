function model = removeNotKnownMets(model, information_predicted_mets)

% addPredictedReactions.m
% 
% Author: Francesco Balzerani
% Email: fbalzerani@tecnun.es
% Date: 21/10/2020

% Remove those predicted metabolites by RetroPath RL for which it's not
% possible to find information in any database and the related reactions.
% For those present in databases, update the name and aliases field

information_predicted_mets = table2cell(information_predicted_mets);
information_predicted_mets(65:67,:) = [];

metsToDelete = [];
rxns_to_delete = [];
for i = 1 : length(information_predicted_mets)
    
    % find the position of the mets
    pos_met = find(ismember(model.metNames,information_predicted_mets(i,2)));
    
    % if is known in DBs, modify the name and write down the ID, otherwise
    % add to the list of mets to remove and the related rxns
    if information_predicted_mets{i,14} == 1
        model.metNames(pos_met) = information_predicted_mets(i,15);
        model.mets(pos_met) = information_predicted_mets(i,15);
        model.metAliases(pos_met) = information_predicted_mets(i,16);
    else    
        pos_rxns = find(model.S(pos_met,:));        
        metsToDelete = [metsToDelete;pos_met];
        rxns_to_delete = [rxns_to_delete,pos_rxns];
    end
end

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

end