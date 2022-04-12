python change_config.py --DB_CACHE True --DB_time 0  --use_cache True --add_Hs True --biosensor True

python Tree.py  \
    --log_file tree.log \
    --itermax 1000  \
    --expansion_width 20 \
    --time_budget 7200 \
    --max_depth 2 \
    --UCT_policy Biochemical_UCT_1 \
    --UCTK 20 \
    --bias_k 0 \
    --k_rave 50 \
    --Rollout_policy Rollout_policy_random_uniform_on_biochemical_multiplication_score \
    --max_rollout 3 \
    --chemical_scoring SubandprodChemicalScorer \
    --virtual_visits 0 \
    --progressive_bias_strategy max_reward  \
    --diameter 10 12 14 16 \
    --c_name pipecolate \
    --c_inchi "InChI=1S/C6H11NO2/c8-6(9)5-3-1-2-4-7-5/h5,7H,1-4H2,(H,8,9)" \
    --folder_to_save pipecolate \
    --EC_filter 1.5.3.7 1.5.3 \
    --biological_score_cut_off 0.1  \
    --substrate_only_score_cut_off 0.9 \
    --chemical_score_cut_off 0.9 \
    --minimal_visit_counts 1
