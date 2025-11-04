import pandas as pd
import os

"""
This scripts aims to look at tasks penalties that strongly co-correlate.
The input is a dataframe of tasks pairs and their correlation (built in R)
For each pair, the scripts goes in the results to look at the reactions used for sample 1 (in filtered files) 
    and computes the % of reactions in common
The script also looks at the genes in the GPR of these used reactions and computes the % of genes in common

This was done before any summary files were produced from the SWAMP pipeline and could eventually be
    modified for the summary produced
    
"""
df_corr = pd.read_excel(
    'C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 15 UMC cardio/stats_UMC cardio/results/test_co_correlation_between_taskscores.xlsx')

# for each pair, go into the result files for sample 1 and extract reactions used for both tasks
# compute % of reactions in common

folder_results_sample_1 = r"C:\Users\inapa\Desktop\results_test_run_CARIM_route_optimization\sample_1_all_tasks_edgeR_BC_data"

df_corr_copy = df_corr
df_corr_copy['perc_rxn_common_1'] = [None for x in range(len(df_corr_copy))]
df_corr_copy['perc_rxn_common_2'] = [None for x in range(len(df_corr_copy))]
df_corr_copy['perc_gene_common_1'] = [None for x in range(len(df_corr_copy))]
df_corr_copy['perc_gene_common_2'] = [None for x in range(len(df_corr_copy))]

from cobra.io import (load_json_model)
from cobra import Model

model_inhouse_v7 = load_json_model(
    "C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/model_modifications/versions/version_7/output/model_inhouse_v7.json")

for row in range(len(df_corr_copy)):
#for row in [0]:

    task_no_1 = df_corr_copy.loc[row, 'var1']
    task_no_2 = df_corr_copy.loc[row, 'var2']

    name_filtered_file_1 = f"filtered_Task_{'%03d' % task_no_1}_Sample001_irr.xlsx"
    name_filtered_file_2 = f"filtered_Task_{'%03d' % task_no_2}_Sample001_irr.xlsx"

    path_1 = os.path.join(folder_results_sample_1, name_filtered_file_1)
    path_2 = os.path.join(folder_results_sample_1, name_filtered_file_2)

    filtered_sample_1_task_no_1 = pd.read_excel(path_1)
    filtered_sample_1_task_no_2 = pd.read_excel(path_2)

    reactions_1 = [rxn.replace('_f', '').replace('_r','') for rxn in filtered_sample_1_task_no_1['ID'] if 'temporary' not in rxn]
    reactions_2 = [rxn.replace('_f', '').replace('_r','')  for rxn in filtered_sample_1_task_no_2['ID'] if 'temporary' not in rxn]

    rxn_common = set(reactions_1).intersection(set(reactions_2))

    perc_rxn_common_1 = (len(rxn_common)/len(reactions_1)) * 100
    perc_rxn_common_2 = (len(rxn_common) / len(reactions_2)) * 100

    df_corr_copy.loc[row, 'perc_rxn_common_1'] = perc_rxn_common_1
    df_corr_copy.loc[row, 'perc_rxn_common_2'] = perc_rxn_common_2

    genes_rxn_1 = [
        gene.id
        for rxn in reactions_1
        for gene in model_inhouse_v7.reactions.get_by_id(rxn).genes
    ]

    genes_rxn_2 = [
        gene.id
        for rxn in reactions_2
        for gene in model_inhouse_v7.reactions.get_by_id(rxn).genes
    ]

    genes_common = set(genes_rxn_1).intersection(genes_rxn_2)
    if len(genes_rxn_1) != 0:
        perc_genes_common_1 = (len(genes_common) / len(genes_rxn_1)) * 100
        df_corr_copy.loc[row, 'perc_gene_common_1'] = perc_genes_common_1
    if len(genes_rxn_2) != 0:
        perc_genes_common_2 = (len(genes_common) / len(genes_rxn_2)) * 100
        df_corr_copy.loc[row, 'perc_gene_common_2'] = perc_genes_common_2

df_corr_copy.to_excel(
    'C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 15 UMC cardio/stats_UMC cardio/results/test_co_correlation_between_taskscores_w_perc_common.xlsx')
