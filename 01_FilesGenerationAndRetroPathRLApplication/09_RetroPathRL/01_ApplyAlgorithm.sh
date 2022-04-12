#!/bin/bash

#Author: Francesco Balzerani
#Contact: fbalzerani@tecnun.es
#Date: 30/01/2020

#This script works to apply the RetroPartRL analysis to find new pathways
#for the phenolic compounds present in Phenol-Explorer.
#The first two lines are necessary to create the files we need to run the pipeline,
#the four files for the rules (with and without H) and the files related to the 
#organisms (our sink).
#The RetroRules's rules without H are avoided to have a bigger metabolic space and put 
#also the sergio's reactions (the metabolic space is a lot bigger without removing the 
#rules for which we can't create the implicit rules.
#The data were prepared with python in a csv file with the name in the first column
#and the InChI ID in the second one.
#We extract the names and the IDs in two variables in order to apply the pipeline
#iteratively to all of the metabolites.
#If the folder of the results and the folders for each metabolites are not created yet
#we create them.
#We specify the "organism_name" == none and biosensor == True in order to apply the
#biosensor study and not the retrosynthesis (the organism_name == none is in order to
#not referred to any specific organism, modifying the file calculate_organism.py in the
#"detectable_compound" variable to insert our sink, in with H as well as in without H). 
#In case of retrosynthesis, put biosensor = False.
#We decided to define the parameter add_Hs = True, because we lose too much information
#removing those rules without implicit form. However, the authors say
#that using rules without H (add_Hs = False), the algorithm perfoms better and faster.
#We parallel the pipeline for N= 20 metabolites each iteration (trying fixing semaphore
#and timeout errors).
#The value of the biological and chemical score cut off are selected as 0.1 and 0.6
#respectively.
#The time_budget has to be high, with 8hours of time budget we should obtain a proper
#analysis of the metabolic space for all the metabolites
#We define a seed for the reproducibility: seed = 1234
#######################################################################################

#The folder (PWD) must be the 09_RetroPathRL folder where are present all the scripts
#we need to apply the pipeline
#The sink file must have the header "name,inchi" and the separator=comma, instead in
#the source file we must remove the header, otherwise we'll have problem with the loop.
#Moreover, in the source file, the separator=\t and the names of the molecules mustn't 
#have an empty space in because otherwise there will be an error in the process, since 
#the algorithm read wrong the name to apply, using the two part of the name separately.
#Remember to put the source,sink and rules files compatible with the unix format.
#######################################################################################



python calculate_rule_sets_similarity.py --rule_address_with_H ./data/rules/rules_with_H_only_substrate.csv --rule_address_without_H ./data/rules/rules_without_H_only_substrate.csv

python calculate_organisms.py

awk 'BEGIN {FS = "\t"} ; {print $1}' $PWD/data/general_source_phen_ex.csv > $PWD/data/names.txt
awk 'BEGIN {FS = "\t"} ; {print $2}' $PWD/data/general_source_phen_ex.csv > $PWD/data/keys.txt

#awk 'BEGIN {FS = "\t"} ; {print $1}' $PWD/data/metabolites_2nd_step.csv > $PWD/data/names.txt
#awk 'BEGIN {FS = "\t"} ; {print $2}' $PWD/data/metabolites_2nd_step.csv > $PWD/data/keys.txt

NAMES=($(cat $PWD/data/names.txt))
KEYS=($(cat $PWD/data/keys.txt))
NUM=${#NAMES[@]}
N=20
SEQ=($(seq $(($N-1)) $N $(($NUM-1))))

if [ ! -d "$PWD/data/results_biosensor" ]; then
        mkdir $PWD/data/results_biosensor
fi

python change_config.py --add_Hs True --biosensor True 

for ((i=0;i<NUM;i++))
do
if [ ! -d "$PWD/data/results_biosensor/$i" ]; then
mkdir $PWD/data/results_biosensor/$i
fi

echo ${NAMES[i]}
echo ${KEYS[i]}
echo $i

python Tree.py \
    --log_file tree.log \
    --itermax 1000 \
    --expansion_width 10 \
    --time_budget 28800 \
    --max_depth 7 \
    --UCT_policy Biochemical_UCT_1 \
    --UCTK 20 \
    --bias_k 0 \
    --k_rave 0 \
    --Rollout_policy Rollout_policy_random_uniform_on_biochemical_multiplication_score \
    --max_rollout 3 \
    --chemical_scoring SubstrateChemicalScorer \
    --virtual_visits 0 \
    --progressive_bias_strategy max_reward \
    --diameter 6 8 10 12 14 16 \
    --organism_name "none" \
    --c_name ${NAMES[i]} \
    --c_inchi ${KEYS[i]} \
    --folder_to_save $PWD/data/results_biosensor/$i \
    --biological_score_cut_off 0.1 \
    --substrate_only_score_cut_off 0.6 \
    --chemical_score_cut_off 0.6 \
    --minimal_visit_counts 1 \
    --seed 1234 \
    --fire_timeout 5 &

for j in ${SEQ[@]}
do
if [ $i -eq $j ]; then
wait
break
fi
done
done
