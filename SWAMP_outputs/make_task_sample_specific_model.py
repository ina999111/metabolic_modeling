
from cobra.io import (save_json_model,
                      load_json_model,)

from cobra import (Model,
                   Reaction,
                   Metabolite)

import pandas as pd

# TODO: change this path for that of the model used for the run
model_inhouse_v7 = load_json_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/model_modifications/versions/version_7/output/model_inhouse_v7.json")

# TODO: change this path for the location of the summary file
df_summary = pd.read_excel(r"C:\Users\inapa\Documents\Repertoire etudiant_IMP\Projets\Projet 15 UMC cardio\MSP_project_period\summary_of_all_filtered_results.xlsx")

def create_task_sample_specific_model(
        model: Model,
        task_number: int,
        sample_number: int
):
    """
    This function can be used to create a task and sample specific model
    That is a small metabolic model containing only the reactions used by a specific sample for a specific task
    The reactions used are taken from a summary file

    :param task_number:
    :param sample_number:
    :return: a task and sample specific model, can be used for mapping in escher
    """

    specific_model = Model() # empty model
    # get list of reactions from the summary df
    row = df_summary[(df_summary['Sample'] == sample_number) & (df_summary['Task'] == task_number)]
    all_rxns_used = row['reactions_used']

    list_all_rxns_used = all_rxns_used.str.split("__").explode().tolist()

    list_non_exchange_rxns =[rxn.replace('_f', '').replace('_r', '') for rxn in list_all_rxns_used if
                   'temporary' not in rxn]

    # creating model
    for reaction_id in list_non_exchange_rxns:
        individual_reaction = model.reactions.get_by_id(reaction_id)
        specific_model.add_reactions([individual_reaction])

    return(specific_model)

# TODO: change the task and sample number
task_and_sample_specific_model = create_task_sample_specific_model(model = model_inhouse_v7,
                                                                   task_number = 18,
                                                                   sample_number = 2)
# TODO: change this path to where you want to save the model
save_json_model(task_and_sample_specific_model, r"C:\Users\inapa\Documents\Repertoire etudiant_IMP\Projets\Projet 15 UMC cardio\MSP_project_period\test_model.json")
