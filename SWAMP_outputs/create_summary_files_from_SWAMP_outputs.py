import os
import pandas as pd

df_tasklist = pd.read_excel(r"P:\FSE_MACSBIO\Maltais-Payette, Ina\DCM project\for_SWaPAM\tasklists\tasklist_MACSBIO_v0_6_20250904\tasklist_MACSBIO_v0_6_20250904.xlsx")

df_tasklist_clean = df_tasklist[['ID','DESCRIPTION','Metabolism Subcategory','SHOULD FAIL']].dropna(how='all')

list_task_ids = list(df_tasklist_clean.loc[:,'ID'])

folder_results = r"C:\Users\inapa\Desktop\results_test_run_CARIM_route_optimization\energy_substrate_tasks_all_samples_combat_data_fixed_fpk"

# as specified in route optimisation options
list_samples = list(range(1,324))

df_scores = pd.DataFrame(index = list_task_ids, columns = list_samples)
df_times = pd.DataFrame(index = list_task_ids, columns = list_samples)
df_reactions_id = pd.DataFrame(index = list_task_ids, columns = list_samples)
df_n_reactions = pd.DataFrame(index = list_task_ids, columns = list_samples)

for task_id in list_task_ids:
    for sample in list_samples:
        name_filtered_file = f"filtered_Task_{'%03d' % task_id}_Sample{'%03d' % sample}_irr.xlsx"
        path = os.path.join(folder_results, name_filtered_file)
        if os.path.exists(path):
            print(path) # to see progression
            df_filtered = pd.read_excel(path, usecols=["ID", "Option_values"])
            df_scores.loc[task_id, sample] = df_filtered.loc[4,'Option_values']
            df_times.loc[task_id, sample] = df_filtered.loc[3, 'Option_values']

            df_reactions_id.loc[task_id, sample] = [rxn for rxn in df_filtered['ID'] if 'temporary' not in rxn]
            df_n_reactions.loc[task_id, sample] = len(df_reactions_id.loc[task_id, sample])


df_scores.to_excel(os.path.join(folder_results, r"_penalties.xlsx"))
df_times.to_excel(os.path.join(folder_results, r"_times.xlsx"))
df_reactions_id.to_excel(os.path.join(folder_results, r"_reactions_id.xlsx"))
df_n_reactions.to_excel(os.path.join(folder_results, r"_n_reactions.xlsx"))
