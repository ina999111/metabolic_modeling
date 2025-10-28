from scripts_IMP.random_code import model_inhouse


def change_formula_id_for_names(
        model: Model,
        formula: str
) -> str:
    dict_id_and_names = dict()
    for metabolite in model.metabolites:
        dict_id_and_names[metabolite.id[:-3]] = metabolite.name

    for id, name in dict_id_and_names.items():
        formula = formula.replace(id, name)

    return (formula)


change_formula_id_for_names(model_inhouse, model_inhouse.reactions.get_by_id('MAR04143').reaction)

list_rxn_ids = ['MAR05998',
                'MAR04388',
                'MAR04926',
                'MAR04143',
                'MAR04141',
                'MAR04852',
                'MAR04139',
                'MAR04101',
                'MAR04363',
                'MAR04365',
                'MAR04368',
                'MAR04373',
                'MAR04391',
                'MAR04375',
                'MAR04377',
                'MAR04381',
                'MAR04856',
                'MAR04521',
                'MAR05027',
                'MAR05029']

for rxn_id in list_rxn_ids:
    formula_id = model_inhouse.reactions.get_by_id(rxn_id).reaction
    formula_name = change_formula_id_for_names(model_inhouse,formula_id)
    print(rxn_id)
    print(formula_name)

rxn_used = ['MAR01377', 'MAR04101', 'MAR04139', 'MAR04141', 'MAR04143', 'MAR04363', 'MAR04365', 'MAR04368', 'MAR04373',
            'MAR04375', 'MAR04377', 'MAR04381', 'MAR04388', 'MAR04391', 'MAR04521', 'MAR04852', 'MAR04856', 'MAR04926',
            'MAR05998', 'MAR09034', 'MAR09079', 'MAR09135']

for rxn_id in rxn_used:
    formula = model_inhouse.reactions.get_by_id(rxn_id).reaction
    formula_name = change_formula_id_for_names(model_inhouse, formula)
    print(formula_name)
