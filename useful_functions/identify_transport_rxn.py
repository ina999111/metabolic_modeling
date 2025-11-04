import os
import pandas as pd
from cobra.io import (load_json_model)
from cobra import Model

def is_transport(model: Model,
                 rxn_id: str
) -> bool:
    """
    Looks at the reaction id in the specified model and evaluates if it is a transport reaction
    defined as having at least 1 metabolite be the same in substrate and product, but with different compartment
    :param model: a cobra model
    :param rxn_id: reaction id as in the model
    :return: true or false
    """

    rxn = model.reactions.get_by_id(rxn_id)
    list_substrates_names = [met.name for met in rxn.reactants]
    list_products_names = [met.name for met in rxn.products]

    list_substrates_compartment = [met.id[-3:] for met in rxn.reactants]
    list_products_compartment = [met.id[-3:] for met in rxn.products]

    # if a substrate is also a product, then look at compartment
    common_met = set(list_substrates_names).intersection(set(list_products_names))

    if bool(common_met) is True:  # if set not empty
        for met_name in common_met:
            substrate_compartment = list_substrates_compartment[list_substrates_names.index(met_name)]
            product_compartment = list_products_compartment[list_products_names.index(met_name)]

            if substrate_compartment != product_compartment:
                return True
            else:
                return False
    else:
        return False

model_inhouse_v7 = load_json_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/model_modifications/versions/version_7/output/model_inhouse_v7.json")


# test
is_transport(model_inhouse_v7, 'MAR00444')
model_inhouse_v7.reactions.get_by_id('MAR00444').reaction


is_transport(model_inhouse_v7, 'MAR06509')
model_inhouse_v7.reactions.get_by_id('MAR06509').reaction
