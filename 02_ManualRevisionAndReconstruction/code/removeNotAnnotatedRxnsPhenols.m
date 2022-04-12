function model = removeNotAnnotatedRxnsPhenols(model,fastFVAName)
    
    if nargin == 1
        fastFVAName = 'removeNotAnnotatedRxnsFastFVAPhenols';
    end

    rxns_without_tax_allowed = find(~cellfun(@isempty,regexp(model.rxnID,'rxnDiet')) | ~cellfun(@isempty,regexp(model.rxnID,'rxnAdd')) | ~cellfun(@isempty,regexp(model.rxnID,'rxnTransport')));
    rxns_without_tax = find(cellfun(@isempty,model.trules));

    rxns_to_delete = rxns_without_tax(~ismember(rxns_without_tax,rxns_without_tax_allowed));

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

%     Add an exchange reaction because we lost it during the
%     fastcoreWeighted. This exchange is added also before fastcore, but
%     had a weight greated than a reaction with taxonomy or without tax but
%     with EC. So to have it in the model, introduce it now manually.
% 
%     1-Feruloyl-D-glucose
%     
    tmp_pos = find(~cellfun(@isempty,regexp(model.metNames,'1-Feruloyl-D-glucose')));
    tmp_id = model.rxnID(~cellfun(@isempty,regexp(model.rxnID,'rxnAddNewEX')));
    tmp_id = strsplit(tmp_id{end},'EX');
    count_n= str2double(tmp_id{2});

    model.rxnID = [model.rxnID;['rxnAddNewEX' num2str(count_n+1)]];
    model.S = [model.S, zeros(size(model.S,1),1)];
    model.S(tmp_pos,end) = -1;
    model.rxns = [model.rxns;strcat('EX_',model.metNames(tmp_pos))];
    model.rxnNames = [model.rxnNames;strcat('EX_',model.metNames(tmp_pos))];
    model.lb = [model.lb;-1000];
    model.ub = [model.ub;0];
    
    n_mets = 1;

    model.c  = [model.c;zeros(n_mets,1)];
    model.rxnGeneMat = [model.rxnGeneMat ; sparse(n_mets, length(model.genes))];
    model.subSystems = [model.subSystems;cell(n_mets,1)];
    model.rules = [model.rules;cell(n_mets,1)];
    model.grRules = [model.grRules;cell(n_mets,1)];
    model.rxnConfidenceScores = [model.rxnConfidenceScores;zeros(n_mets,1)];
    model.rxnECNumbers = [model.rxnECNumbers;cell(n_mets,1)];
    model.rxnKEGGID = [model.rxnKEGGID;cell(n_mets,1)];
    model.correctedRxns = [model.correctedRxns;zeros(n_mets,1)];
    model.rxnCode = [model.rxnCode;cell(n_mets,1)];
    model.rxnStoichiometry = [model.rxnStoichiometry;cell(n_mets,1)];
    model.rxnIsTransport = [model.rxnIsTransport;cell(n_mets,1)];
    model.rxnEquation = [model.rxnEquation;cell(n_mets,1)];
    model.rxnDefinition = [model.rxnDefinition;cell(n_mets,1)];
    model.rxnReversibility = [model.rxnReversibility;cell(n_mets,1)];
    model.rxnDirection = [model.rxnDirection;cell(n_mets,1)];
    model.rxnAbstractReaction = [model.rxnAbstractReaction;cell(n_mets,1)];
    model.rxnPathways = [model.rxnPathways;cell(n_mets,1)];
    model.rxnAliases = [model.rxnAliases;cell(n_mets,1)];
    model.rxnDeltaG = [model.rxnDeltaG;cell(n_mets,1)];
    model.rxnDeltaGErr = [model.rxnDeltaGErr;cell(n_mets,1)];
    model.rxnCompoundIDs = [model.rxnCompoundIDs;cell(n_mets,1)];
    model.rxnStatus = [model.rxnStatus;cell(n_mets,1)];
    model.rxnIsObsolete = [model.rxnIsObsolete;cell(n_mets,1)];
    model.rxnLinkedReaction = [model.rxnLinkedReaction;cell(n_mets,1)];
    model.rxnNotes = [model.rxnNotes;cell(n_mets,1)];
    model.rxnSource = [model.rxnSource;cell(n_mets,1)];
    model.rxnOntology = [model.rxnOntology;cell(n_mets,1)];
    model.comments = [model.comments;cell(n_mets,1)];
    model.citations = [model.citations;cell(n_mets,1)];
    model.trRules = [model.trRules; cell(n_mets,1)];
    model.trules = [model.trules; cell(n_mets,1)];
    model.rxnTaxMat = [model.rxnTaxMat; sparse(n_mets,length(model.taxonomy))];

    if exist([fastFVAName '.mat'], 'file')
        load([fastFVAName '.mat'], 'vMin', 'vMax')
        if length(vMin)~=length(model.rxns) || length(vMax)~=length(model.rxns)
            error('The fastFVA fluxes do not correspond to the model reactions')
        end
    else
        [vMin, vMax] = fastFVA(model);
        save([fastFVAName '.mat'], 'vMin', 'vMax')
    end
    
    epsilon = 1e-08;
    rxns_to_delete = find(abs(vMax)<epsilon & abs(vMin)<epsilon);
    
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
    
    idx = find(model.ub<0);
    model.S(:,idx) = -1*model.S(:,idx);
    tmp_lb = model.lb;
    model.lb(idx) = -model.ub(idx);
    model.ub(idx) = -tmp_lb(idx);
end


