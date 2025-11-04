import cobra
from cobra.io import (save_json_model,
                      save_matlab_model,
                      write_sbml_model,
                      load_json_model,
                      load_matlab_model,
                      read_sbml_model)

from cobra import (Model,
                   Reaction,
                   Metabolite)
import os

import pandas as pd
import re
import mygene
import warnings

model_inhouse = load_json_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/model_modifications/versions/version_IMP/output/model_inhouse_IMP_20250307.json")

dict_ids_and_names = {metabolite.id:metabolite.name[:-3] for metabolite in model_inhouse.metabolites}

def change_ids_for_names_in_tasklist(
        cobra_model: Model,
        tasklist: pd.DataFrame
) -> pd.DataFrame:
    """
    from a standard excel tasklist, replaces the names of metabolites inputs and outputs with the corresponding id
    also verifies if the metabolites exist in that compartment and throws an error message if not
    :param tasklist_name: standard tasklist with columns "IN" and "OUT" for input and output metabolites
    :return:
    """
    list_inputs_ids = tasklist['IN'].values.tolist()
    list_inputs_names = [None for x in list_inputs_ids]
    for idx in range(len(list_inputs_ids)):
        id = list_inputs_ids[idx]
        if name == '':
            list_inputs_names[idx] = ''
        else:
            name = dict_ids_and_names[id[:-3]]
            name_with_comp = name + id[-3:]
            list_inputs_name[idx] = id_with_comp

            if id not in [metabolite.id for metabolite in cobra_model.metabolites]:
                warning_message = f"the metabolite id {id_with_comp} was not found in the model"
                warnings.warn(warning_message)

    list_output_names = tasklist['OUT'].values.tolist()
    list_output_id = [None for x in list_output_names]
    for idx in range(len(list_output_names)):
        name = list_output_names[idx]
        if name == '':
            list_output_id[idx] = ''
        else:
            id = dict_names_and_ids[name[:-3]]
            id_with_comp = id + name[-3:]
            list_output_id[idx] = id_with_comp

            if id_with_comp not in [metabolite.id for metabolite in cobra_model.metabolites]:
                warning_message = f"the metabolite id {id_with_comp} was not found in the model"
                warnings.warn(warning_message)

    tasklist['IN'] = list_inputs_id
    tasklist['OUT'] = list_output_id

    return tasklist
