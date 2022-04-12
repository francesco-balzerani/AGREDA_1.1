% extractionRxnsMetsUB.m
%
% Author: Francesco Balzerani
% Email: fbalzerani@tecnun.es
% Date: 25/11/2019
% 
% Script to extract the stoichiometric matrix and all the information
% regarding the reactions of the universal database with taxonomy
% information. Extract the information for the metabolites as well.

clear 
close all
clc

load('.\input\agora_seed_merged_sergio.mat')

model = model_merged_sergio;

ex_rxns = find(sum(model.S~=0, 1)==1);

pos_trules = find(~cellfun(@isempty,model.trules));

S = model.S(:,pos_trules(~ismember(pos_trules,ex_rxns)));

mets_to_delete = find(sum(S~=0, 2)==0);

S(mets_to_delete,:) = [];

% save('./output/stoichiometric_matrix.mat','S','-v7')

varNames = {'Names','ID','EC','Formula','Taxonomy','lb','ub'};

% modify the name for some reactions because they are defined as NaN
nan_name = find(cellfun(@(x) x,cellfun(@length,cellfun(@isnan,model.rxns,'UniformOutput',false),'UniformOutput',false))==1);
model.rxns(nan_name) = cellfun(@strcat, repmat({'Name_'},length(nan_name),1), cellfun(@(x) num2str(x), ...
    num2cell((1:length(nan_name))'),'UniformOutput',false),'UniformOutput',false);
nan_name = find(cellfun(@(x) x,cellfun(@length,cellfun(@isnan,model.rxnNames,'UniformOutput',false),'UniformOutput',false))==1);
model.rxnNames(nan_name) = cellfun(@strcat, repmat({'Name_'},length(nan_name),1), cellfun(@(x) num2str(x), ...
    num2cell((1:length(nan_name))'),'UniformOutput',false),'UniformOutput',false);


library_rxns = cell(length(model.trules(pos_trules(~ismember(pos_trules,ex_rxns)))),4);

library_rxns(:,1) = model.rxnNames(pos_trules(~ismember(pos_trules,ex_rxns)));
library_rxns(:,2) = model.rxnID(pos_trules(~ismember(pos_trules,ex_rxns)));
library_rxns(:,3) = model.rxnECNumbers(pos_trules(~ismember(pos_trules,ex_rxns)));
library_rxns(:,4) = printRxnFormula(model, 'rxnAbbrList', model.rxns(pos_trules(~ismember(pos_trules,ex_rxns))),'metNameFlag',true);
library_rxns(:,5) = model.trules(pos_trules(~ismember(pos_trules,ex_rxns)));
lb = model.lb(pos_trules(~ismember(pos_trules,ex_rxns)));
ub = model.ub(pos_trules(~ismember(pos_trules,ex_rxns)));

P = table(library_rxns(:,1),library_rxns(:,2),library_rxns(:,3),library_rxns(:,4),library_rxns(:,5),lb,ub,'VariableNames',varNames);

% writetable(P,'./output/library_rxns.xlsx')

varNames = {'Names','ID','InchiKey','Formula','Smiles'};

metNames = model.metNames;
metNames(mets_to_delete) = [];
metID = model.metID;
metID(mets_to_delete) = [];
metInchiKey = model.metInchiKey;
metInchiKey(mets_to_delete) = [];
metFormulas = model.metFormulas;
metFormulas(mets_to_delete) = [];
metSmiles = model.metSmiles;
metSmiles(mets_to_delete) = [];

library_mets = cell(length(metNames),5);

library_mets(:,1) = metNames;
library_mets(:,2) = metID;
library_mets(:,3) = metInchiKey; 
library_mets(:,4) = metFormulas;
library_mets(:,5) = metSmiles;

T = table(library_mets(:,1),library_mets(:,2),library_mets(:,3),library_mets(:,4),library_mets(:,5),'VariableNames',varNames);

% writetable(T,'./output/library_mets.xlsx')

temp = model_merged_sergio.rxnID(~cellfun(@isempty,regexp(model_merged_sergio.rxnID,'rxnAdd\d+')));
literature = table2cell(readtable('./input/literatureKnowledge.xlsx', 'Sheet', 'Reactions'));
for i = 1 : length(literature)
    tmp_cpd = strsplit(literature{i,1},';');
    tmp_stoic = strsplit(literature{i,2},';');
    substrate = cell(1,1);
    product = cell(1,1);
    count_s = 0;
    count_p = 0;
    for j = 1 : length(tmp_stoic)
        if str2double(tmp_stoic(j))<0
            count_s = count_s + 1;
            tmp_name = model_merged_sergio.metNames(find(ismember(model_merged_sergio.metID,tmp_cpd(j))));
            substrate(count_s,1) = unique(tmp_name);
        else
            count_p = count_p + 1;
            tmp_name = model_merged_sergio.metNames(find(ismember(model_merged_sergio.metID,tmp_cpd(j))));
            if ismember('Phloretate',tmp_name)
                tmp_name = {'Phloretate'};
            end
            product(count_p,1) = unique(tmp_name);
        end
    end
    substrate = strjoin(substrate,' + ');
    product = strjoin(product, ' + ');
    temp{i,2} = strjoin({substrate, product}, ' -> ');
end

% writetable(cell2table(temp, 'VariableNames', {'rxnID','Formula'}),'./output/sergio_reactions_to_generate_rules.xlsx');

list_cpd = cell(length(literature),1);
for i = 1 : length(literature)
    list_cpd{i,1} = regexp(literature{i,1},'cpd\d+','match')';
end

list_cpd = unique(cat(1,list_cpd{:}));
for i = 1 : length(list_cpd)
    tmp_name = model_merged_sergio.metNames(find(ismember(model_merged_sergio.metID,list_cpd(i))));
    if ismember('Phloretate',tmp_name)
        tmp_name = {'Phloretate'};
    end
    list_cpd{i,2} = unique(tmp_name);
end

% writetable(cell2table(list_cpd, 'VariableNames', {'cpdID','Name'}),'./output/cpd_name_sergio_reactions.xlsx')

names_id = [model_merged_sergio.metNames,model_merged_sergio.metID];
% writetable(cell2table(names_id, 'VariableNames', {'metNames','metID'}),'./output/names_cpd_model_merged_sergio.xlsx')

