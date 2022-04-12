function [finalModel, PhenolsAddedRxns] = applyFastcorePhenols(model, core, cpdID_rxnID, sense, options)

% 
% USAGE:
%
%   [finalModel, iDietAddedRxns] = applyFastcorePhenols(model, core, sense, predicted_mets, options)
%
% INPUTS:
%   model:      Metabolic model structure (COBRA Toolbox format).
%   core:       Index of the reactions that formed the core for fastcore
%               algorithm.
%   sense:      Variable to indicate the direction of the pathway for each
%               Phenol. Can take the next values: uptake, outside
%               or both.
%
% OPTIONAL INPUT:
%   options:    Structure with fields:
%
%                   *.modelName - Name to attach to the model.
%                   *.plantRxns - Set of reaction IDs that are known to
%                    belong to plants. These reactions have a greater
%                    weight than the others.
%
% Authors: Francesco Balzerani
% Email: fbalzerani@tecnun.es
% Last modification: 26/10/2020

%Set parameters
if nargin == 4
    options.modelName = 'modelFastcorePhenols';
    options.plantRxns = [];
    options.fastFVAName = fullfile('.', 'files_for_model', 'optionsFastcorePhenols', ['fastFVAPhenols_', modelName]);
    options.blockedPhenolsRxnsUp = [];
    options.blockedPhenolsRxnsOut = [];
    options.rxnPhenolsID = 'rxnAddPhenEX';
else
    if isfield(options, 'modelName')
        modelName = options.modelName;
    else
        modelName = 'modelFastcorePhenols';
    end
    if isfield(options, 'plantRxns')
        plantRxns = options.plantRxns;
    else
        plantRxns = [];
    end
    if isfield(options, 'fastFVAName')
        fastFVAName = options.fastFVAName;
    else
        fastFVAName = fullfile('.', 'output', 'optionsFastcorePhenols', ['fastFVAPhenols_', modelName]);
    end
    if isfield(options, 'blockedPhenolsRxnsUp')
        blockedPhenolsRxnsUp = options.blockedPhenolsRxnsUp;
    else
        blockedPhenolsRxnsUp = [];
    end
    if isfield(options, 'blockedPhenolsRxnsOut')
        blockedPhenolsRxnsOut = options.blockedPhenolsRxnsOut;
    else
        blockedPhenolsRxnsOut = [];
    end
    if isfield(options, 'rxnPhenolsID')
        rxnPhenolsID = options.rxnPhenolsID;
    else
        rxnPhenolsID = 'rxnAddPhenEX';
    end
end

%Check the value of sense variable
if strcmp(sense,'uptake')
    origSense = 'uptake';
    nIter = 1;
elseif strcmp(sense,'outside')
    origSense = 'outside';
    nIter = 1;
elseif strcmp(sense,'both')
    origSense = 'uptake';
    oppositeSense = 'outside';
    nIter = 2;
else
    error('The sense is not correct. Options are uptake, outside or both.')
end
    
% modify the name for some reactions because they are defined as NaN
nan_name = find(cellfun(@(x) x,cellfun(@length,cellfun(@isnan,model.rxns,'UniformOutput',false),'UniformOutput',false))==1);
model.rxns(nan_name) = cellfun(@strcat, repmat({'Name_'},length(nan_name),1), cellfun(@(x) num2str(x), ...
    num2cell((1:length(nan_name))'),'UniformOutput',false),'UniformOutput',false);
nan_name = find(cellfun(@(x) x,cellfun(@length,cellfun(@isnan,model.rxnNames,'UniformOutput',false),'UniformOutput',false))==1);
model.rxnNames(nan_name) = cellfun(@strcat, repmat({'Name_'},length(nan_name),1), cellfun(@(x) num2str(x), ...
    num2cell((1:length(nan_name))'),'UniformOutput',false),'UniformOutput',false);

%Define core reactions
model.core = zeros(length(model.rxns),1);
model.core(core) = 1;

%Define reactions related to plants
if ~isempty(plantRxns)
    model.plantRxns = zeros(length(model.rxns),1);
    model.plantRxns(ismember(model.rxnID,plantRxns)) = 1;
end

model.modelName = modelName;
model_orig = model;

%Split the model in order to apply the fastFVA
revs_f = find(model.lb < 0 & model.ub > 0);
revRxnsRelation = [revs_f (length(model.rxns)+1:length(model.rxns)+length(revs_f))'];

model.S = [model.S, -model.S(:, revs_f)];
model.ub = [model.ub; -model.lb(revs_f)];
model.lb = zeros(size(model.ub));
model.rxns = [model.rxns; model.rxns(revs_f)];
model.c = [model.c; model.c(revs_f)];
model.rules = [model.rules; model.rules(revs_f)];

if exist([fastFVAName '.mat'], 'file')
    load([fastFVAName '.mat'], 'vMin', 'vMax')
    if length(vMin)~=length(model.rxns) || length(vMax)~=length(model.rxns)
        error('The fastFVA fluxes do not correspond to the model reactions')
    end
else
    [vMin, vMax] = fastFVA(model);
    save([fastFVAName '.mat'], 'vMin', 'vMax')
end

%Change reaction bounds and delete blocked reactions from the model
epsilon = 1e-06;
blockedRxns = find(abs(vMin)<=epsilon & vMax<=epsilon);

allRevRxns = [revRxnsRelation(:,1);revRxnsRelation(:,2)];
irrevRxnsToDelete = blockedRxns(~ismember(blockedRxns,allRevRxns));
revRxnsToChange = blockedRxns(ismember(blockedRxns,allRevRxns));

for i = 1:length(revRxnsToChange)
    if ismember(revRxnsToChange(i), revRxnsRelation(:,1))
        model_orig.ub(revRxnsToChange(i)) = 0;
    elseif ismember(revRxnsToChange(i), revRxnsRelation(:,2))
        pos = revRxnsRelation(ismember(revRxnsRelation(:,2),revRxnsToChange(i)),1);
        model_orig.lb(pos) = 0;
    end
end

model_orig.S(:,irrevRxnsToDelete) = [];
model_orig.lb(irrevRxnsToDelete) = [];
model_orig.ub(irrevRxnsToDelete) = [];
model_orig.rxns(irrevRxnsToDelete) = [];
model_orig.rxnNames(irrevRxnsToDelete) = [];
model_orig.c(irrevRxnsToDelete) = [];
model_orig.trules(irrevRxnsToDelete) = [];
model_orig.rules(irrevRxnsToDelete) = [];
model_orig.subSystems(irrevRxnsToDelete) = [];
model_orig.grRules(irrevRxnsToDelete) = [];
model_orig.comments(irrevRxnsToDelete) = [];
model_orig.citations(irrevRxnsToDelete) = [];
model_orig.rxnConfidenceScores(irrevRxnsToDelete) = [];
model_orig.correctedRxns(irrevRxnsToDelete) = [];
model_orig.rxnECNumbers(irrevRxnsToDelete) = [];
model_orig.rxnKEGGID(irrevRxnsToDelete) = [];
model_orig.rxnID(irrevRxnsToDelete) = [];
model_orig.rxnCode(irrevRxnsToDelete) = [];
model_orig.rxnStoichiometry(irrevRxnsToDelete) = [];
model_orig.rxnIsTransport(irrevRxnsToDelete) = [];
model_orig.rxnEquation(irrevRxnsToDelete) = [];
model_orig.rxnDefinition(irrevRxnsToDelete) = [];
model_orig.rxnReversibility(irrevRxnsToDelete) = [];
model_orig.rxnDirection(irrevRxnsToDelete) = [];
model_orig.rxnAbstractReaction(irrevRxnsToDelete) = [];
model_orig.rxnPathways(irrevRxnsToDelete) = [];
model_orig.rxnAliases(irrevRxnsToDelete) = [];
model_orig.rxnDeltaG(irrevRxnsToDelete) = [];
model_orig.rxnDeltaGErr(irrevRxnsToDelete) = [];
model_orig.rxnCompoundIDs(irrevRxnsToDelete) = [];
model_orig.rxnStatus(irrevRxnsToDelete) = [];
model_orig.rxnIsObsolete(irrevRxnsToDelete) = [];
model_orig.rxnLinkedReaction(irrevRxnsToDelete) = [];
model_orig.rxnNotes(irrevRxnsToDelete) = [];
model_orig.rxnSource(irrevRxnsToDelete) = [];
model_orig.rxnOntology(irrevRxnsToDelete) = [];
model_orig.trRules(irrevRxnsToDelete) = [];
model_orig.rxnGeneMat(irrevRxnsToDelete,:) = [];
model_orig.rxnTaxMat(irrevRxnsToDelete,:) = [];
model_orig.core(irrevRxnsToDelete) = [];
if isfield(model_orig, 'plantRxns')
    model_orig.plantRxns(irrevRxnsToDelete,:) = [];
end

model = model_orig;
finalModel = model_orig;

%Calculate iDiet blocked metabolites for the chosen sense
PhenolsRxnsIdx = find(~cellfun(@isempty,regexp(model.rxnID,rxnPhenolsID)));
PhenolsRxns = model.rxnID(PhenolsRxnsIdx);
for j = 1:nIter
    if j == 1
        sense = origSense;
    else
        sense = oppositeSense;
    end
    if strcmp(sense,'uptake')
        if isempty(blockedPhenolsRxnsUp)
            nPhenols = length(PhenolsRxnsIdx);
            blockedPhenolsRxnsUp = zeros(nPhenols,1);
            for i = 1:nPhenols
                modelTmp = model;
                modelTmp.c = zeros(length(modelTmp.c),1);
                modelTmp.c(PhenolsRxnsIdx(i)) = 1;
                cplex = FBA(modelTmp,'minimize');
                if abs(cplex.Solution.objval) < 1e-08
                    blockedPhenolsRxnsUp(i) = 1;
                end
            end
            save(fullfile('.','output','optionsFastcorePhenols',['blockedPhenolsRxnsUp_' modelName '.mat']), 'blockedPhenolsRxnsUp')
            PhenolsRxnsOfIntUp = PhenolsRxns(~logical(blockedPhenolsRxnsUp));
        else
            PhenolsRxnsOfIntUp = PhenolsRxns(~logical(blockedPhenolsRxnsUp));
        end
    elseif strcmp(sense,'outside')
        if isempty(blockedPhenolsRxnsOut)
            nPhenols = length(PhenolsRxnsIdx);
            blockedPhenolsRxnsOut = zeros(nPhenols,1);
            for i = 1:nPhenols
                modelTmp = model;
                modelTmp.c = zeros(length(modelTmp.c),1);
                modelTmp.c(PhenolsRxnsIdx(i)) = 1;
                cplex = FBA(modelTmp,'maximize');
                if abs(cplex.Solution.objval) < 1e-08
                    blockedPhenolsRxnsOut(i) = 1;
                end
            end
            save(fullfile('.','output','optionsFastcorePhenols',['blockedPhenolsRxnsOut_' modelName '.mat']), 'blockedPhenolsRxnsOut')
            PhenolsRxnsOfIntOut = PhenolsRxns(~logical(blockedPhenolsRxnsOut));
        else
            PhenolsRxnsOfIntOut = PhenolsRxns(~logical(blockedPhenolsRxnsOut));
        end
    end
end

%Split the model
nOrigRxns = length(model.rxns);
origCore = find(model.core);
revs_f = find(model.lb < 0 & model.ub > 0);
revRxnsRelation = [revs_f (length(model.rxns)+1:length(model.rxns)+length(revs_f))'];

model.S = [model.S, -model.S(:, revs_f)];
model.ub = [model.ub; -model.lb(revs_f)];
model.lb = zeros(size(model.ub));
model.rxns = [model.rxns; model.rxns(revs_f)];
model.c = [model.c; model.c(revs_f)];
model.rules = [model.rules; model.rules(revs_f)];
model.trules = [model.trules; model.trules(revs_f)];
model.rxnID = [model.rxnID; model.rxnID(revs_f)];
model.core = [model.core; model.core(revs_f)];
if isfield(model, 'plantRxns')
    model.plantRxns = [model.plantRxns; model.plantRxns(revs_f)];
end

%Select the index of the Phenols reactions
PhenolsRxnsIdx = find(~cellfun(@isempty,regexp(model.rxnID,rxnPhenolsID)));

%Calculate iteratively the reactions added by each Phenols
for j = 1:nIter
    if j == 1
        sense = origSense;
    else
        sense = oppositeSense;
    end
    %Select those that are not blocked
    if strcmp(sense,'uptake')
        PhenolsRxnsUpIdx = PhenolsRxnsIdx(sum(model.S(:,PhenolsRxnsIdx), 1) == 1);
        PhenolsRxnsOutIdx = PhenolsRxnsIdx(sum(model.S(:,PhenolsRxnsIdx), 1) == -1);
        PhenolsRxnsUp = model.rxnID(PhenolsRxnsUpIdx);
        PhenolsRxnsOut = model.rxnID(PhenolsRxnsOutIdx);
        PhenolsRxnsUpIdx = PhenolsRxnsUpIdx(ismember(PhenolsRxnsUp,PhenolsRxnsOfIntUp));
        PhenolsRxnsOutIdx = PhenolsRxnsOutIdx(ismember(PhenolsRxnsOut,PhenolsRxnsOfIntUp));
        nPhenols = length(PhenolsRxnsOfIntUp);
    else
        PhenolsRxnsUpIdx = PhenolsRxnsIdx(sum(model.S(:,PhenolsRxnsIdx), 1) == 1);
        PhenolsRxnsOutIdx = PhenolsRxnsIdx(sum(model.S(:,PhenolsRxnsIdx), 1) == -1);
        PhenolsRxnsUp = model.rxnID(PhenolsRxnsUpIdx);
        PhenolsRxnsOut = model.rxnID(PhenolsRxnsOutIdx);
        PhenolsRxnsUpIdx = PhenolsRxnsUpIdx(ismember(PhenolsRxnsUp,PhenolsRxnsOfIntOut));
        PhenolsRxnsOutIdx = PhenolsRxnsOutIdx(ismember(PhenolsRxnsOut,PhenolsRxnsOfIntOut));
        nPhenols = length(PhenolsRxnsOfIntOut);
    end
    PhenolsAddedRxns = cell(nPhenols,8);
    for i = 1:nPhenols
        clc
        disp([num2str(i) '/' num2str(nPhenols) ' sense:' sense]);
        
        %If the Phenol reaction was already in the core, avoid the fastcore
        if strcmp(sense,'uptake')
            if model.core(PhenolsRxnsUpIdx(i)) == 1
                PhenolsAddedRxns{i,1} = model.rxns(PhenolsRxnsUpIdx(i));
                PhenolsAddedRxns{i,8} = {'CONNECTED'};
                continue;
            end
        elseif strcmp(sense,'outside')
            if model.core(PhenolsRxnsOutIdx(i)) == 1
                PhenolsAddedRxns{i,1} = model.rxns(PhenolsRxnsOutIdx(i));
                PhenolsAddedRxns{i,8} = {'CONNECTED'};
                continue;
            end
        end

        weights = ones(length(model.rxns),1)*100;
        weights(~cellfun(@isempty,model.rxnECNumbers)) = 30;
        if isfield(model, 'plantRxns')
            weights(logical(model.plantRxns)) = 1000;
        end
        weights(~cellfun(@isempty,model.trules)) = 0.1;
        if strcmp(sense,'uptake')
            met_id = model.metID(find(model.S(:,PhenolsRxnsUpIdx(i))));
            related_rxns = find(ismember(model.rxnID,cpdID_rxnID(ismember(cpdID_rxnID(:,1),met_id),2)));
            tmp_phenol = find(model.S(:,PhenolsRxnsUpIdx(i)));
            related_rxns_up = [];
            for k = 1 : length(related_rxns)
                if model.S(tmp_phenol,related_rxns(k)) > 0
                    continue
                else
                    related_rxns_up = [related_rxns_up; related_rxns(k)];
                end
            end
            
            model.core(related_rxns_up) = 1;
            model.core(PhenolsRxnsUpIdx(i)) = 1;
            weights(PhenolsRxnsOutIdx(i),1) = 100000;
        elseif strcmp(sense,'outside')
            related_rxns = find(model.S(find(model.S(:,PhenolsRxnsOutIdx(i))),:));
            related_rxns = related_rxns(~cellfun(@isempty,regexp(model.rxnID(related_rxns),'rxnPhenAdd')));
            tmp_phenol = find(model.S(:,PhenolsRxnsOutIdx(i)));
            related_rxns_out = [];
            for k = 1 : length(related_rxns)
                if model.S(tmp_phenol,related_rxns(k)) < 0
                    continue
                else
                    related_rxns_out = [related_rxns_out; related_rxns(k)];
                end
            end

            model.core(related_rxns_out) = 1;
            model.core(PhenolsRxnsOutIdx(i)) = 1;
            weights(PhenolsRxnsUpIdx(i),1) = 100000;
        end
        weights(logical(model.core)) = 0;

        fastcoreResult = fastCoreWeighted(find(model.core),model,weights);
        coreRxnBoolSplitted = false(length(model.rxns),1);
        coreRxnBoolSplitted(fastcoreResult) = true;

        firstSenseRxnBool = coreRxnBoolSplitted(1:nOrigRxns);
        secondSenseRxnBool = false(nOrigRxns,1);
        tmp = fastcoreResult(fastcoreResult>nOrigRxns);
        tmp = revRxnsRelation(ismember(revRxnsRelation(:,2),tmp),1);
        secondSenseRxnBool(tmp) = true;
        coreRxnBoolOrig = firstSenseRxnBool | secondSenseRxnBool;

        rxnsFound = find(coreRxnBoolOrig);
        PhenolsAddedRxnsIdx = rxnsFound(~ismember(rxnsFound,origCore));

        %Define the reactions added by each Phenols
        PhenolsAddedRxns{i,1} = model.rxns(PhenolsRxnsOutIdx(i));
        PhenolsAddedRxns{i,2} = model.rxnID(PhenolsAddedRxnsIdx);
        PhenolsAddedRxns{i,3} = model.rxns(PhenolsAddedRxnsIdx);
        PhenolsAddedRxns{i,4} = model.rxnAliases(PhenolsAddedRxnsIdx);
        PhenolsAddedRxns{i,5} = printRxnFormula(model,'rxnAbbrList',model.rxns(PhenolsAddedRxnsIdx),'metNameFlag',true);
        clc
        PhenolsAddedRxns{i,6} = model.rxnECNumbers(PhenolsAddedRxnsIdx);
        PhenolsAddedRxns{i,7} = model.trules(PhenolsAddedRxnsIdx);
        rxnsWithoutTax = PhenolsAddedRxnsIdx(cellfun(@isempty,model.trules(PhenolsAddedRxnsIdx)));
        admittedRxns = (~cellfun(@isempty,regexp(model.rxnID(rxnsWithoutTax),'rxnDiet'))) | (~cellfun(@isempty,regexp(model.rxnID(rxnsWithoutTax),'rxnAdd')));
        if sum(admittedRxns) == length(rxnsWithoutTax)
            PhenolsAddedRxns{i,8} = {'CONNECTED'};
        else
            PhenolsAddedRxns{i,8} = {'NOT CONNECTED'};
        end
        
        if strcmp(sense,'uptake')
            model.core(related_rxns_up) = 0;
            model.core(PhenolsRxnsUpIdx(i)) = 0;
        elseif strcmp(sense,'outside')
            model.core(related_rxns_out) = 0;
            model.core(PhenolsRxnsOutIdx(i)) = 0;
        end      
    end
    if j == 1
        PhenolsAddedRxnsTmp = PhenolsAddedRxns;
    end
end

%Build the final model
if nIter == 2
    PhenolsAddedRxns = [PhenolsAddedRxnsTmp; PhenolsAddedRxns];
end
newRxns = unique(cat(1,PhenolsAddedRxns{:,2}));
newRxnsIdx = find(ismember(finalModel.rxnID,newRxns));
rxnsToDelete = true(length(finalModel.rxns),1);
rxnsToDelete([origCore; newRxnsIdx]) = false;

finalModel.S(:,rxnsToDelete) = [];
finalModel.rxns(rxnsToDelete) = [];
finalModel.lb(rxnsToDelete) = [];
finalModel.ub(rxnsToDelete) = [];
finalModel.c(rxnsToDelete) = [];
finalModel.rules(rxnsToDelete) = [];
finalModel.rxnNames(rxnsToDelete) = [];
finalModel.subSystems(rxnsToDelete) = [];
finalModel.grRules(rxnsToDelete) = [];
finalModel.comments(rxnsToDelete) = [];
finalModel.citations(rxnsToDelete) = [];
finalModel.rxnConfidenceScores(rxnsToDelete) = [];
finalModel.correctedRxns(rxnsToDelete) = [];
finalModel.rxnECNumbers(rxnsToDelete) = [];
finalModel.rxnKEGGID(rxnsToDelete) = [];
finalModel.rxnID(rxnsToDelete) = [];
finalModel.rxnCode(rxnsToDelete) = [];
finalModel.rxnStoichiometry(rxnsToDelete) = [];
finalModel.rxnIsTransport(rxnsToDelete) = [];
finalModel.rxnEquation(rxnsToDelete) = [];
finalModel.rxnDefinition(rxnsToDelete) = [];
finalModel.rxnReversibility(rxnsToDelete) = [];
finalModel.rxnDirection(rxnsToDelete) = [];
finalModel.rxnAbstractReaction(rxnsToDelete) = [];
finalModel.rxnPathways(rxnsToDelete) = [];
finalModel.rxnAliases(rxnsToDelete) = [];
finalModel.rxnDeltaG(rxnsToDelete) = [];
finalModel.rxnDeltaGErr(rxnsToDelete) = [];
finalModel.rxnCompoundIDs(rxnsToDelete) = [];
finalModel.rxnStatus(rxnsToDelete) = [];
finalModel.rxnIsObsolete(rxnsToDelete) = [];
finalModel.rxnLinkedReaction(rxnsToDelete) = [];
finalModel.rxnNotes(rxnsToDelete) = [];
finalModel.rxnSource(rxnsToDelete) = [];
finalModel.rxnOntology(rxnsToDelete) = [];
finalModel.trules(rxnsToDelete) = [];
finalModel.trRules(rxnsToDelete) = [];
finalModel.rxnGeneMat(rxnsToDelete,:) = [];
finalModel.rxnTaxMat(rxnsToDelete,:) = [];

metsToDelete = find(sum(finalModel.S~=0, 2)==0);

finalModel.S(metsToDelete,:) = [];
finalModel.mets(metsToDelete) = [];
finalModel.b(metsToDelete) = [];
finalModel.csense(metsToDelete) = [];
finalModel.metNames(metsToDelete) = [];
finalModel.metCharges(metsToDelete) = [];
finalModel.metFormulas(metsToDelete) = [];
finalModel.metChEBIID(metsToDelete) = [];
finalModel.metKEGGID(metsToDelete) = [];
finalModel.metPubChemID(metsToDelete) = [];
finalModel.metInChIString(metsToDelete) = [];
finalModel.metHMDBID(metsToDelete) = [];
finalModel.metSmiles(metsToDelete) = [];
finalModel.metID(metsToDelete) = [];
finalModel.metsTax(metsToDelete) = [];
finalModel.metMass(metsToDelete) = [];
finalModel.metSource(metsToDelete) = [];
finalModel.metInchiKey(metsToDelete) = [];
finalModel.metIsCore(metsToDelete) = [];
finalModel.metIsObsolete(metsToDelete) = [];
finalModel.metLinkedCompound(metsToDelete) = [];
finalModel.metIsCofactor(metsToDelete) = [];
finalModel.metDeltaG(metsToDelete) = [];
finalModel.metDeltaGErr(metsToDelete) = [];
finalModel.metPKA(metsToDelete) = [];
finalModel.metPKB(metsToDelete) = [];
finalModel.metAbstractCompound(metsToDelete) = [];
finalModel.metComprisedOf(metsToDelete) = [];
finalModel.metAliases(metsToDelete) = [];
finalModel.metOntology(metsToDelete) = [];

%Remove unnecessary fields
finalModel = rmfield(finalModel,'core');
if isfield(finalModel,'plantRxns')
    finalModel = rmfield(finalModel,'plantRxns');
end

end