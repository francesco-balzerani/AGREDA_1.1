% ExtractPathways.m
%
% Author: Francesco Balzerani
% Email: fbalzerani@tecnun.es
% Date: 08/06/2021
% 
% Calculate all the pathways that reach metabolites of the sink

clc
clear

initCobraToolbox(false)
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64') % change the path to the proper one related to the position of cplex
addpath(genpath('C:\Users\fbalzerani\cobratoolbox')) % change the path to the proper one related to the position of the cobratoolbox

common = table2cell(readtable('./input/CommonAGREDA10.xlsx'));
name_inchikey = table2cell(readtable('../01_FilesGenerationAndRetroPathRLApplication/output/name_inchikey.xlsx'));
name_inchikey_phen_exp = table2cell(readtable('../01_FilesGenerationAndRetroPathRLApplication/output/name_inchikey_phen_exp.csv'));
for i = 1 : length(name_inchikey_phen_exp)
    pos_phen = find(ismember(common(:,2),name_inchikey_phen_exp(i,1)));
    if ~isempty(pos_phen)
        name_inchikey_phen_exp{i,1} = common{pos_phen,1};
    end
end
name_inchikey = [name_inchikey(:,[1,3]);name_inchikey_phen_exp];

load('../01_FilesGenerationAndRetroPathRLApplication/input/EC_x_tax.mat')
Rule_x_ID = table2cell(readtable('../01_FilesGenerationAndRetroPathRLApplication/output/Rule_x_ID.xlsx'));
sources = table2cell(readtable('../01_FilesGenerationAndRetroPathRLApplication/output/general_source_phen_ex.csv', 'ReadVariableNames',0));
n_sources = length(sources);
whole_metabolites_information = table2cell(readtable(fullfile('.','output','whole_metabolites_information.xlsx')));

load('../01_FilesGenerationAndRetroPathRLApplication/input/agora_seed_merged_sergio.mat')
ID_x_tax = cell(length(model_merged_sergio.rxnID),2);
ID_x_tax(:,1) = model_merged_sergio.rxnID;
ID_x_tax(:,2) = model_merged_sergio.trules;

rxns_id_formula = table2cell(readtable('../01_FilesGenerationAndRetroPathRLApplication/output/formulas_rules.xlsx'));
path = fullfile('.','output','results_to_extract_CS/');
files = dir([path, 'Metabolite_*']);
files = {files(:).name}';
n_files = length(files);

count = 0;
count_new = 0;
cccc = 0;

final_table = cell(n_sources,1);

for i = 1 : n_files
    disp(['iteration ' num2str(i) ' of ' num2str(n_files)])
    disp(files{i})
    
    tmp = split(files{i},'_');
    
    % Load the files for each metabolite
    S_matrix = importdata([path,files{i},'/S_matrix.csv']);
    mets = readtable([path,files{i},'/metabolites_table.xlsx']);
    rxns = readtable([path,files{i},'/reactions_table.xlsx']);

    mets = table2cell(mets);
    rxns = table2cell(rxns);
   
    % Build the model
    model = createModel();
    model.rxns = rxns(:,1);
    for m = 1 : size(mets,1)
        if isempty(mets{m,1})
            pos_name = find(ismember(name_inchikey(:,2),mets{m,3}));
            if isempty(pos_name)
                model.mets{m,1} = mets{m,3};
            else
                model.mets{m,1} = name_inchikey{pos_name,1};
            end
        else
            model.mets{m,1} = mets{m,1};
        end
    end

    model.S = S_matrix;
    model.ub = 1000 * ones(size(model.rxns));
    model.lb = zeros(size(model.rxns));
    model.c = zeros(size(model.rxns));
    model.b = zeros(size(model.mets));
    model.csense = char(zeros(size(model.mets)));
    model.csense(:) = 'E';
    model.rules = cell(size(model.rxns));
    model.rules(:) = {''};
    model.c(end) = 1;
    model.metFormula = cellfun(@(x) x{2},cellfun(@strsplit,mets(:,2),repmat({'/'},length(mets(:,2)),1),'UniformOutput',false),'UniformOutput',false);
    model.metFormula = cellfun(@strrep,model.metFormula,repmat({'p+1'},length(model.metFormula),1),repmat({'H'},length(model.metFormula),1),'UniformOutput',false);
    model = buildRxnGeneMat(model);
    
    for posss = 1 : length(model.mets)
        tmp_formula = model_merged_sergio.metFormulas(ismember(model_merged_sergio.metNames, model.mets(posss)));
        if iscell(tmp_formula) && ~isempty(tmp_formula)
            tmp_formula = tmp_formula{1};
        end
        if ~isempty(tmp_formula) && ~isequal(tmp_formula, 'null')
            model.metFormula{posss} = tmp_formula;
        end           
    end

    [vmin,vmax] = fastFVA(model,0);
    index_active = find(vmax~=0);
    result_reactions = cell(1,9);
    count_active_not_exchange = 0;
    
    tmp_model = model;
    tmp_model.ub(end) = 0;
    [~,tmp_vmax] = fastFVA(tmp_model,0);
    tmp_index_active = find(tmp_vmax~=0);
    cycle_rxns = [];
    for j = 1 : length(tmp_index_active)
        if ~cellfun(@isempty,regexp(tmp_model.rxns(tmp_index_active(j)),'Exchange', 'match'))
            continue
        else
            cycle_rxns = [cycle_rxns; tmp_index_active(j)];
        end
    end
    index_active(ismember(index_active,cycle_rxns)) = [];
    for j = 1 : length(index_active)
        if ~cellfun(@isempty,regexp(model.rxns(index_active(j)),'Exchange', 'match'))
            continue
        end
        count_active_not_exchange = count_active_not_exchange + 1;
%       metabolite name
        result_reactions{count_active_not_exchange,1} = sources{str2double(tmp(2)),1};
%       predicted formula
        substrate = model.mets(find(model.S(:,index_active(j)) == -1));
        products = join(model.mets(find(model.S(:,index_active(j)) == 1)),' + ');
        result_reactions{count_active_not_exchange,2} = join([substrate,products],' -> ');
%       ec number
        result_reactions{count_active_not_exchange,3} = regexp(rxns{index_active(j),2},'\d+.\d+.\d+.[\d-+]','match');

%       taxonomy
        if ~isempty(result_reactions{count_active_not_exchange,3})
            tmp_EC = regexp(result_reactions{count_active_not_exchange,3},'\d+.\d+.\d+.[\d-+]','match');
            ec_tax = EC_x_tax(ismember(EC_x_tax(:,1),tmp_EC{1}),2);
        else
            ec_tax = [];
        end
        if isempty(regexp(model.rxns{index_active(j)},'M\w+','match'))
            tmp_name = regexp(model.rxns{index_active(j)},'rxnAdd\d+','match');
            for k = 1 : length(tmp_name)
                tmp_trules{1,k} = model_merged_sergio.trules(ismember(model_merged_sergio.rxnID,tmp_name{k}));
            end
        elseif isempty(regexp(model.rxns{index_active(j)},'rxnAdd\w+','match'))
            tmp_name = regexp(model.rxns{index_active(j)},'M\w+','match');
            remove = [];
            for k = 1:length(tmp_name)
                tmp_ID = Rule_x_ID(ismember(Rule_x_ID(:,1),tmp_name{k}),2);
                tmp_trules{1,k} = ID_x_tax(ismember(ID_x_tax(:,1),tmp_ID),2);
                if isempty(tmp_trules{1,k}) || isempty(tmp_trules{1,k}{1})
                    remove = [remove;k];
                end
            end
            tmp_trules(remove) = [];
        else
            tmp_name = regexp(model.rxns{index_active(j)},'rxnAdd\d+','match');
            for k = 1 : length(tmp_name)
                tmp_trules_add{1,k} = model_merged_sergio.trules(ismember(model_merged_sergio.rxnID,tmp_name{k}));
            end
            remove = [];
            tmp_name = regexp(model.rxns{index_active(j)},'M\w+','match');
            for k = 1:length(tmp_name)
                tmp_ID = Rule_x_ID(ismember(Rule_x_ID(:,1),tmp_name{k}),2);
                tmp_trules_mnxr{1,k} = ID_x_tax(ismember(ID_x_tax(:,1),tmp_ID),2);
                if isempty(tmp_trules_mnxr{1,k}) || isempty(tmp_trules_mnxr{1,k}{1})
                    remove = [remove;k];
                end
            end
            tmp_trules_mnxr(remove) = [];
            tmp_trules = [tmp_trules_add, tmp_trules_mnxr];
        end
        if ~isempty(tmp_trules)
            tmp_trules = [tmp_trules, ec_tax];
            temp = cat(1,tmp_trules{:});
            temp = strjoin(temp,' | ');
            temp = unique(split(temp,' | '));
            temp = strjoin(temp,' | ');
            result_reactions{count_active_not_exchange,4} = temp;
        else
            result_reactions{count_active_not_exchange,4} = ec_tax;
        end
        
        clear tmp_trules tmp_trules_add tmp_trules_mnxr
               
%       template_reaction, lb, ub
        if isempty(regexp(model.rxns{index_active(j)},'MNXR\d+','match'))
            if ~isempty(regexp(model.rxns{index_active(j)},'rxnAdd\d+_r_\d','match'))
                tmp_rxns = regexp(model.rxns{index_active(j)},'rxnAdd\d+_r_\d','match');
            else
                tmp_rxns = regexp(model.rxns{index_active(j)},'rxnAdd\d+','match');
            end
        else
            tmp_rxns = regexp(model.rxns{index_active(j)},'MNXR\d+','match');
        end
        template = [];
        lb = [];
        ub = [];
        for r = 1 : length(tmp_rxns)
            template = [template; unique(rxns_id_formula(ismember(rxns_id_formula(:,1),tmp_rxns{r}),2))];
            lb = [lb; rxns_id_formula(ismember(rxns_id_formula(:,1),tmp_rxns{r}),3)];
            ub = [ub; rxns_id_formula(ismember(rxns_id_formula(:,1),tmp_rxns{r}),4)];
        end
        
        result_reactions{count_active_not_exchange,6} = template;
        result_reactions{count_active_not_exchange,8} = lb;
        result_reactions{count_active_not_exchange,9} = ub;
        
        result_reactions{count_active_not_exchange,10} = tmp_rxns;
        
%       score
        result_reactions{count_active_not_exchange,11} = rxns{index_active(j),4};

%       formula_mets
        formula_sub = model.metFormula(ismember(model.mets,substrate));
        formula_prod = [];
        prods = split(products, ' + ');
        for p = 1 : length(prods)
            formula_prod = [formula_prod, model.metFormula(ismember(model.mets,prods(p)))];
        end
        prod_comp = join(formula_prod,' + ');
        result_reactions{count_active_not_exchange,7} = join([formula_sub,prod_comp],' -> ');
    end
    final_table{str2double(tmp(2)),1} = sources{str2double(tmp(2)),1};
    final_table{str2double(tmp(2)),2} = result_reactions;
    pos = find(ismember(common(:,2),sources{str2double(tmp(2)),1}));
    if isempty(pos)
        final_table{str2double(tmp(2)),3} = 'new_phenol';
    else
        final_table{str2double(tmp(2)),3} = 'already_present';
    end
end

index_new = [];
for i = 1 : length(final_table)
    if isempty(final_table{i,2}) || isempty(final_table{i,2}{1})
        continue
    else
        if isequal(final_table{i,3},'new_phenol')
            index_new = [index_new; i];
        end
    end
end

complete_results = {};
position = 0;

for i = 1 : length(index_new)
   for  j = 1 : size(final_table{index_new(i),2},1)
       position = position + 1;
       complete_results{position,1} = final_table{index_new(i),2}{j,1};
       complete_results{position,2} = final_table{index_new(i),2}{j,2}{1};
       complete_results{position,3} = strjoin(final_table{index_new(i),2}{j,3},' || ');
       if iscell(final_table{index_new(i),2}{j,4})
           complete_results{position,4} = final_table{index_new(i),2}{j,4}{1};
       else
           complete_results{position,4} = final_table{index_new(i),2}{j,4};
       end

       complete_results{position,6} = strjoin(final_table{index_new(i),2}{j,6},' || ');
       complete_results{position,7} = final_table{index_new(i),2}{j,7}{1};
       complete_results{position,10} = strjoin(final_table{index_new(i),2}{j,10},' || ');
       for k = 1 : size(final_table{index_new(i),2}{j,8},1)
           final_table{index_new(i),2}{j,8}{k} = num2str(final_table{index_new(i),2}{j,8}{k});
           final_table{index_new(i),2}{j,9}{k} = num2str(final_table{index_new(i),2}{j,9}{k});
       end
       complete_results{position,8} = strjoin(final_table{index_new(i),2}{j,8},' || ');
       complete_results{position,9} = strjoin(final_table{index_new(i),2}{j,9},' || ');
       complete_results{position,11} = final_table{index_new(i),2}{j,11};
   end 
end

varNames = {'Name_metabolite','Formula','Formula_mets','Template_reaction','ID_template','EC','Taxonomy','lb','ub','ChemicalScore'};
P = table(complete_results(:,1),complete_results(:,2),complete_results(:,7),complete_results(:,6),complete_results(:,10),complete_results(:,3),complete_results(:,4),complete_results(:,8),complete_results(:,9),complete_results(:,11),'VariableNames',varNames);

% writetable(P,fullfile('.','output','whole_pathways_CS.xlsx'))
% checked manually the reactions. For some metabolites there is the
% inchikey even though they are present in the sink, so add their name and
% check their metID. (file: mets_not_known.csv)

index_already = [];
for i = 1 : length(final_table)
    if isempty(final_table{i,2}) || isempty(final_table{i,2}{1})
        continue
    else
        if isequal(final_table{i,3},'already_present')
            index_already = [index_already; i];
        end
    end
end

already_results = {};
position = 0;

for i = 1 : length(index_already)
   for  j = 1 : size(final_table{index_already(i),2},1)
       position = position + 1;
       already_results{position,1} = final_table{index_already(i),2}{j,1};
       already_results{position,2} = final_table{index_already(i),2}{j,2}{1};
       already_results{position,3} = strjoin(final_table{index_already(i),2}{j,3},' || ');
       if iscell(final_table{index_already(i),2}{j,4})
           already_results{position,4} = final_table{index_already(i),2}{j,4}{1};
       else
           already_results{position,4} = final_table{index_already(i),2}{j,4};
       end
       already_results{position,6} = strjoin(final_table{index_already(i),2}{j,6},' || ');
       already_results{position,7} = final_table{index_already(i),2}{j,7}{1};
       already_results{position,10} = strjoin(final_table{index_already(i),2}{j,10},' || ');
       for k = 1 : size(final_table{index_already(i),2}{j,8},1)
           final_table{index_already(i),2}{j,8}{k} = num2str(final_table{index_already(i),2}{j,8}{k});
           final_table{index_already(i),2}{j,9}{k} = num2str(final_table{index_already(i),2}{j,9}{k});
       end
       already_results{position,8} = strjoin(final_table{index_already(i),2}{j,8},' || ');
       already_results{position,9} = strjoin(final_table{index_already(i),2}{j,9},' || ');
   end 
end

varNames = {'Name_metabolite','Formula','Formula_mets','Template_reaction','ID_template','EC','Taxonomy','lb','ub'};
P = table(already_results(:,1),already_results(:,2),already_results(:,7),already_results(:,6),already_results(:,10),already_results(:,3),already_results(:,4),already_results(:,8),already_results(:,9),'VariableNames',varNames);

% writetable(P,fullfile('.','output','predictions_control.xlsx'))
% checked manually the reactions. For some metabolites there is the
% inchikey even though they are present in the sink, so add their name and
% check their metID. (file: mets_not_known.csv)

% revise the file manually, balance the reaction and save a new excel "./output/final_table.xlsx"
