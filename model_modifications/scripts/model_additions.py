
import sys

sys.path.append(r"C:\Users\inapa\PycharmProjects\SWAMP")

from src.cobrapy_fork._cobra import (
    load_matlab_model,
    save_json_model,
)


from src.cobrapy_fork.io import (
    load_json_model,
    # load_matlab_model, # NOT OKAY
    # save_json_model, # NOT OKAY, but doesn't lead to an issue in the case of subsystems
)


from cobra import (Model,
                   Reaction,
                   Metabolite)
import os

import pandas as pd
import re
import mygene
import warnings



# TODO add names for excel files designation __reactions __metabolites
# TODO make so that it works with all model file formats
# todo add csv and tsv support for metabolites and reactions to change files

# TODO: modify to handle not only additions, but also removal and change to mets/rxns
## the setup for this would be to have tables with a "action" flag stating if it is an addition/change/removal
## as well as columns as follows (here for mets but similar for rxns)
## before_metabolite_ID
## before_metabolite_name
## before_metabolite_compartment
## before_metabolite_formula
## before_metabolite_charge
## after_metabolite_ID
## after_metabolite_name
## after_metabolite_compartment
## after_metabolite_formula
## after_metabolite_charge
## project_notes : to make the link with the tasks used to test the rxns/mets added
## for additions, only the after_ columns are filled
## for removals, only the before_ columns are filled
## for changes, both before_ and after_ columns are filled



'''
This scripts adds metabolites and reactions to a cobra model using the following inputs:
    - a cobra model (matlab format) 
    - excel file (.xlsx) containing the metabolites to add, must include the following columns: 
        metabolite_ID: starts with MAM
        metabolite_name: no non-alphanumeric characters
        metabolite_compartment: single letter, lowercase e.g. 'c'
        metabolite_formula: chemical formula. e.g. C6H12O6
        metabolite_charge: integer
    - excel file containing the reactions to add, must include the following columns:
        reaction_ID: starts with MAR
        reactions_name: no non-alphanumeric characters
        reaction_subsystem: used to specify the project for which the 
            reactions are added: no non-alphanumeric characters 
        reaction_formula_readable: use readable metabolites names (exactly as in the 
            model or the excel file containing the metabolites to add), 
            followed by compartment between []: 
            These should always be forward '-->' or reversible '<=>'
            e.g. Glucose[c] + ATP[c] --> Glucose-6-phosphate[c] + ADP[c] + H[c]
            Reverisble reactions should have <=> instead of -->
            
        reaction_GPR_readable: Containing gene symbols, separated by 'or'/'and' (lower case!).
             If there is no GPR to add, it should contain 'no gene rule':
             e.g. (G6PD or G6PDH) and (G6PD and G6PDH)
             e.g. no gene rule
'''

def create_log_for_metabolites(
        cobra_model: Model,
        list_metabolites_added: list
) -> pd.DataFrame:
    """
    for each metabolite added to the model, finds existing metabolites with the same formula
    as well as similar formula (-1H, -2H, -3H, -4H, +1H and +2H)
    outputs a log data frame
    :param cobra_model: a cobra model BEFORE adding metabolites
    :param list_metabolites_added: list of dictionaries for each added metabolite
    the dictionaries are created with the Metabolite function
    :return: a dataframe with added metabolites in rows and the metabolites with same or similar formulas in columns
    columns(n=7) = same formula as well as formula -1H, -2H, -3H, -4H, +1H and +2H
    """
    amount_of_H_differences_negative = 4
    amount_of_H_differences_positive = 2

    metabolites_with_same_or_similar_formula = {
        metabolite: {
            "same_formula": [],
            **{f"similar_formula_+{i}H": [] for i in range(1, amount_of_H_differences_positive + 1)},
            **{f"similar_formula_-{i}H": [] for i in range(1, amount_of_H_differences_negative + 1)}
        }
        for metabolite in list_metabolites_added
    }

    model_metabolites = cobra_model.metabolites  # removed self.
    pattern_digits_after_H = re.compile(r'^(.*?H)(\d+)(.*)$')

    for metabolite in list_metabolites_added:

        formula = metabolite.formula
        id_without_compartment = metabolite.id[:metabolite.id.rfind(metabolite.compartment)]
        if formula is None:
            continue
        match = pattern_digits_after_H.match(formula)
        if match:
            if not len(match.groups()) == 3:
                continue
            else:
                number_of_H = int(match.group(2))
                base_formula = match.group(1) + match.group(3)
                for other_metabolite in model_metabolites:
                    other_formula = other_metabolite.formula
                    other_id_without_compartment = other_metabolite.id[
                                                   :other_metabolite.id.rfind(other_metabolite.compartment)]
                    if other_formula is None:
                        continue
                    elif id_without_compartment == other_id_without_compartment:
                        continue
                    other_match = pattern_digits_after_H.match(other_formula)
                    if other_match:
                        if not len(other_match.groups()) == 3:
                            continue
                        else:
                            other_number_of_H = int(other_match.group(2))
                            other_base_formula = other_match.group(1) + other_match.group(3)
                            difference_h = number_of_H - other_number_of_H
                            if difference_h > amount_of_H_differences_positive or difference_h < -amount_of_H_differences_negative:
                                continue
                            if base_formula == other_base_formula:
                                if number_of_H == other_number_of_H:
                                    metabolites_in_list = metabolites_with_same_or_similar_formula[metabolite][
                                        "same_formula"]
                                    metabolites_without_compartment_in_list = [met.id[:met.id.rfind(met.compartment)] for
                                                                        met in metabolites_in_list]
                                    if other_id_without_compartment not in metabolites_without_compartment_in_list:
                                        metabolites_with_same_or_similar_formula[metabolite]["same_formula"].append(
                                            other_metabolite)
                                elif difference_h > 0:
                                    metabolites_in_list = metabolites_with_same_or_similar_formula[metabolite][
                                        f"similar_formula_+{difference_h}H"]
                                    metabolites_without_compartment_in_list = [met.id[:met.id.rfind(met.compartment)] for
                                                                        met in metabolites_in_list]
                                    if other_id_without_compartment not in metabolites_without_compartment_in_list:
                                        metabolites_with_same_or_similar_formula[metabolite][
                                            f"similar_formula_+{difference_h}H"].append(other_metabolite)
                                elif difference_h < 0:
                                    metabolites_in_list = metabolites_with_same_or_similar_formula[metabolite][
                                        f"similar_formula_{difference_h}H"]
                                    metabolites_without_compartment_in_list = [met.id[:met.id.rfind(met.compartment)] for
                                                                        met in metabolites_in_list]
                                    if other_id_without_compartment not in metabolites_without_compartment_in_list:
                                        metabolites_with_same_or_similar_formula[metabolite][
                                            f"similar_formula_{difference_h}H"].append(other_metabolite)


    dict_flatten = {}
    for metabolite in metabolites_with_same_or_similar_formula:

        metabolite_dict = metabolites_with_same_or_similar_formula[metabolite]

        list_of_list_ids = []
        for type_similarity in metabolite_dict:

            list_ids = [similar_metabolite.id for similar_metabolite in metabolite_dict[type_similarity]]
            list_of_list_ids.append(list_ids)

            dict_flatten[metabolite.id] = list_of_list_ids

    dataframe_similar_metabolites = pd.DataFrame.from_dict(dict_flatten, orient='index',
                                                           columns=['same_formula',
                                                                    'similar_formula_+1H',
                                                                    'similar_formula_+2H',
                                                                    'similar_formula_-1H',
                                                                    'similar_formula_-2H',
                                                                    'similar_formula_-3H',
                                                                    'similar_formula_-4H'])
    return dataframe_similar_metabolites

#test
#log_metabolites_test = create_log_for_metabolites(cobra_model, list_metabolites_added)
#log_metabolites_test.to_excel("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/model_modifications/versions/version_test/log_mets_test.xlsx")

def add_metabolites_to_model(
        cobra_model: Model,
        dataframe_metabolites_to_change:pd.DataFrame,
) -> pd.DataFrame:
    """
    adds the metabolites specified in the dataframe to the cobra model
    Changes model in-place
    :param model: cobra model
    :param dataframe_metabolites_to_change:
    excel file (.xlsx) containing the metabolites to add, must include the following columns: metabolite_ID:
        metabolite_name, metabolite_compartment, metabolite_formula and metabolite_charge. See info above for specific
        format of each.
    :return: a cobra model with the new metabolites added and a metabolite log
    """

    model_new = cobra_model.copy()
    start_number_metabolites = len(model_new.metabolites)

    list_metabolites_added = [None for x in range(len(dataframe_metabolites_to_change))]

    for idx in range(len(dataframe_metabolites_to_change)):
        metabolite_to_add = Metabolite(
            dataframe_metabolites_to_change.loc[idx,'after_metabolite_ID'],
            name = dataframe_metabolites_to_change.loc[idx,'after_metabolite_name'],
            compartment = dataframe_metabolites_to_change.loc[idx,'after_metabolite_compartment'],
            formula = dataframe_metabolites_to_change.loc[idx,'after_metabolite_formula'],
            charge = int(dataframe_metabolites_to_change.loc[idx,'after_metabolite_charge']) # Int needed to stop int64 from pd read_excel
            )

        model_new.add_metabolites(metabolite_to_add)

        list_metabolites_added[idx] = metabolite_to_add

    log_metabolites = create_log_for_metabolites(cobra_model, list_metabolites_added)

    end_number_metabolites = len(model_new.metabolites)
    added_number_metabolites = start_number_metabolites - end_number_metabolites
    print(added_number_metabolites, 'metabolites were added to the model')

    return model_new, log_metabolites

def create_dict_stochiometry_by_reaction(
        cobra_model: Model,
        formula: str,
) -> dict:
    """
    creates a dictionary of the stochiometry of each metabolite from a formula string
    Steps:
        create list of metabolites names (without [comp])
        create list of compartments
        replace metabolites names with IDs (exclude [comp])
        add original [comp] to ID
        find if there is a coefficient before the metabolite name
            if yes, stochiometry = coefficient
            if not, stochiometry = 1
        determine if the metabolite is before (substrate) or after (product) the arrow in the formula to determine if
            the stochiometry should be positive or negative

    :param model: Model: a cobra model
    :param formula: str: a string representing the reaction of the formula to add. metabolites names should be spelled
        exactly as the metabolite names in the model or in the metabolites added (aka readable and NOT "MAMXXXXX".
        metabolite names should be following by the compartment between square brackets eg [c], with no space between
        name and the compartment. The arrow in the formula should be "->" for irreversible reactions and "<=>" for
        reversible reactions. Coefficients, "+" and arrows should be separated from names by a space.
        e.g. NADH[m] + 6 H+[m] + lauric acid[c] -> 4 7,12-Dimethylbenz[A]Anthracene-3,4-Diol-1,2-Epoxide[m]
    :return: dict: a dictionary specific to that reaction, with metabolite ids as keys and stochiometry as values
    """

    # split the formula on +, closing brackets and arrow
    # this is necessary because some metabolite names in the model contain brackets, spaces and/or arrows
    pattern_split = r'\]\s\+\s|\]\s\-\-\>\s|\]\s\<\-\-\s|\]\s\<\=\>\s' # matches '] + ' or '] -> ' or '] <- ' or '] <=> '
    split_formula = re.split(pattern_split, formula)
    split_formula = list(filter(bool, split_formula)) # removing empty strings introduced by re.split

    # the very last element of the list will still have [compartment] at the end (instead of [compartment
    # so remove it to make all metabolites the same
    # this is only necessary is there is >1 metabolite in the formula (not exchange reaction)
    if len(split_formula) > 1:
        split_formula[len(split_formula) - 1] = split_formula[len(split_formula) - 1][:-1]

    # add closing bracket to the compartments
    list_names_with_compartements_and_coefficients = [split_formula[x] + ']' for x in range(len(split_formula))]

    # remove coefficients to get names with compartemnt only
    pattern_coefficient = r'(?:^[0-9]\s|^[0-9]*[.][0-9]*\s)' # matches numbers + the space after in looking like 2 or 1.96, only catches numbers at beginning of string
    list_names_with_compartements = [re.sub(pattern_coefficient, '', list_names_with_compartements_and_coefficients[x]) for x in
                  range(len(list_names_with_compartements_and_coefficients))]

    list_location_names = [None for x in range(len(list_names_with_compartements_and_coefficients))]
    for idx in range(len(list_names_with_compartements_and_coefficients)):
        escaped_name = re.escape(list_names_with_compartements_and_coefficients[idx])
        pattern_name = r"\+\s" + escaped_name + r"|\>\s" + escaped_name + r"|^" + escaped_name
        location = re.search(pattern_name, formula).start()
        list_location_names[idx] = location

    # find location of arrow (preceded by ] because some metabolite names contain arrows)
    pattern_arrow = '(?:\]\s\-\-\>\s)|(?:\]\s\<\-\-\s)|(?:\]\s\<\=\>\s)'
    position_arrow = [match.start() for match in re.finditer(pattern_arrow, formula)][0]

    # find coefficients (digit + space), if there is none, coefficient = 1
    list_coefficients = [None for x in range(len(list_names_with_compartements_and_coefficients))]
    for idx in range(len(list_names_with_compartements_and_coefficients)):
        match_object = re.search(pattern_coefficient, list_names_with_compartements_and_coefficients[idx])
        if match_object is None:
            coefficient = float(1)
        elif match_object is not None:
            coefficient = float(
                match_object.group()[:-1]) # the space after the digit is also caught by pattern, so removing it here
        list_coefficients[idx] = coefficient

    # substrates are negative and products are positive
    list_stochiometry = [None for x in range(len(list_coefficients))]
    for idx in range(len(list_coefficients)):
        if list_location_names[idx] < position_arrow:  # if name before the arrow, substrate, so negative coefficient
            list_stochiometry[idx] = -1 * list_coefficients[idx]
        elif list_location_names[idx] > position_arrow:  # if name after the arrow, product, so positive coefficient
            list_stochiometry[idx] = list_coefficients[idx]

    # create dictionnary for all metabolits in model with names[compartement] as keys and ids as values
    list_all_names_with_compartements_in_model = [cobra_model.metabolites[x].name + cobra_model.metabolites[x].id[-3:] for x in range(len(cobra_model.metabolites))]
    list_all_ids_in_model = [cobra_model.metabolites[x].id for x in range(len(cobra_model.metabolites))]
    dict_model_names_and_ids = dict(zip(list_all_names_with_compartements_in_model,list_all_ids_in_model))

    # replace names by ids
    list_ids_with_compartement = [None for x in range(len(list_names_with_compartements))]
    for idx in range(len(list_names_with_compartements)):

        if list_names_with_compartements[idx] not in dict_model_names_and_ids:
            warning_message = f"{list_names_with_compartements[idx]} in formula {formula} was not found in the model"
            warnings.warn(warning_message)
            list_ids_with_compartement[idx] = 'not found'

        elif list_names_with_compartements[idx] in dict_model_names_and_ids:
            list_ids_with_compartement[idx] = dict_model_names_and_ids[list_names_with_compartements[idx]]

    # build final dict
    dict_stochiometry = dict(zip(list_ids_with_compartement, list_stochiometry))

    return dict_stochiometry

# test works
# create_dict_stochiometry_by_reaction(model, "glucose[m] -> 2 3-dehyrdo[m] + 4-hydro[c] + 2 ATP[c] + 2 ADP[c] + H+[c]")

def find_ensemblid_from_gene_symbols(
        GPR: str
) -> str:
    """
    For strings with gene symbols, find the corresponding ENSEMBL id.
    If multiple potential matches are found, take the one with the highest score (first in list so idx = 0)
    If the gene has more than 1 ENSEMBL ID, identify the one that is mapped on the primary assembly = the position of the gene is on a chromosome instead of a scaffold
    :param GPR: string of gene symbols with 'or'/'and' (lower case!) between them, () also ok
    :return: string with the same 'or'/'and'/() but ENSEMBL ids instead of gene symbols
    """
    if GPR.lower() == 'no gene rule':
        return ''

    if GPR.startswith('ENSG'): # if GPR already in GPR id, keep it like that
        return(GPR)

    # raising a warning if the GPR rule contains something else than uppercase letters, digits, hyphens, 'or', 'and' and parentheses
    pattern_allowed_in_GPR = '((?:[A-Z]|[0-9]|\-)+|\sand\s|\sor\s|\(|\))'
    if re.sub(pattern_allowed_in_GPR, '', GPR) != '':  # if something else than the allowed characters
        warning_message = (f"The GPR string seems to contain characters that do not belong to a gene symbol:"
                           f"{re.sub(pattern_allowed_in_GPR, '', GPR)}."
                           f"Make sure that the GPR rule is in the right format")
        warnings.warn(warning_message)

    mg = mygene.MyGeneInfo()

    pattern_gene_symbol = '([A-Z]|[1-9]|\-)+\w'  # matches upper case letters, numbers and hyphens

    list_gene_symbols = [match.group() for match in re.finditer(pattern_gene_symbol, GPR)]
    gene_dict = {gene_symbol: '(no ensembl id found)' for gene_symbol in list_gene_symbols}

    for gene_symbol in list_gene_symbols:

        output_query = mg.query(gene_symbol, scopes='symbol', fields='ensembl.gene,genomic_pos.chr', species='human')

        # if no match is found, the gene symbol is skipped and the value will remain None
        if output_query['max_score'] == None:
            warning_message = f'no ENSEMBL id found for gene symbol {gene_symbol}'
            warnings.warn(warning_message)
            continue

        # we want to know if a gene symbol has more than 1 ENSBL ids
        number_of_ensembl_ids = len(output_query['hits'][0]['ensembl'])

        if number_of_ensembl_ids == 1:  # if only 1 ENSEMBL id is found for that gene symbol, take that ENSEMBL id
            gene_dict[gene_symbol] = output_query['hits'][0]['ensembl']['gene']

        elif number_of_ensembl_ids > 1:  # if more than 1 are found, found the correct one: the one mapped on a chromosome and not a scaffold
            list_multiple_ensembl_ids = [output_query['hits'][0]['ensembl'][x]['gene'] for x in
                                         range(number_of_ensembl_ids)]
            list_multiple_genomic_pos = [output_query['hits'][0]['genomic_pos'][x]['chr'] for x in
                                         range(
                                             number_of_ensembl_ids)]  # this genomic pos is the chromosome number (eg. '1') or an id for the scaffold (eg. 'HSCHR1_1_CTG31')

            # verify is there is a chromosome number in the genomic positions
            list_bolean_genomic_pos_is_numeric = [list_multiple_genomic_pos[x].isnumeric() for x in
                                                  range(len(list_multiple_genomic_pos))]
            if True in list_bolean_genomic_pos_is_numeric:
                for index, value in enumerate(list_bolean_genomic_pos_is_numeric):  # find the index that has a chromosome number
                    if value == True:  # chromosome is only a number scaffold contains letters
                        gene_dict[gene_symbol] = list_multiple_ensembl_ids[index]
                warnings_message = f'the gene symbol {gene_symbol} has returned {number_of_ensembl_ids} ensembl ids: {list_multiple_ensembl_ids}. {gene_dict[gene_symbol]} was chosen'
                warnings.warn(warnings_message)

            elif True not in list_bolean_genomic_pos_is_numeric:  # if there is no numeric in the list of chromosome position
                warning_message = f'none of the ensembl ids found for {gene_symbol} are mapped on a chromosome'
                warnings.warn(warning_message)
                continue

    for key in gene_dict:
        pattern_sub = key + r'\b'  # creating a regex by add the \b to the gene symbol. This allows for example FABP1 to match only FABP1 and not both FABP1 and FABP12
        GPR = re.sub(pattern_sub, gene_dict[key], GPR)

    return GPR

def create_log_for_genes(
        cobra_model: Model,
        list_all_ensembl_ids_added: list
) -> pd.DataFrame:
    """
    using a string of all the GPR rules that are added, find the unique ensembl ids, verify whether they were already in the model, and if so, which reactions they were used in
    :param model: cobra model (before adding the new gpr rules !)
    :param concat_all_GPR: str: a concatenated string of all GPR rules (with 'or'/'and'/() and duplicates) genes should be ensembl ids
    :return: pd.DataFrame: with EMSEMBL IDs as column names and the list of reaction ids below it
    """

    unique_ensembl_ids = set(list_all_ensembl_ids_added)

    list_all_ensembl_ids_in_model = [gene.id for gene in cobra_model.genes]

    dict_ensembl_ids_added = dict()
    for ensembl_id in unique_ensembl_ids:
        if ensembl_id in list_all_ensembl_ids_in_model:
                dict_ensembl_ids_added[ensembl_id] = [reaction.id for reaction in cobra_model.genes.get_by_id(ensembl_id).reactions]

    df_ensembl_ids_added = pd.DataFrame.from_dict(dict_ensembl_ids_added, orient = 'index').transpose()

    return df_ensembl_ids_added

#create_log_for_genes(model, list_all_GPR_added)


def get_lb_and_ub_from_formula_string(
        formula: str,
) -> tuple:

    # determine if the reaction is reversible or irreversebile (arrow <=> or ->)
    directionality_arrow_pattern = '(\-\-\>|\<\=\>|\<\-\-)'
    string_arrow = [match.group() for match in re.finditer(directionality_arrow_pattern, formula)][0]

    if string_arrow == '-->':
        lb = 0
        ub = 1000
    elif string_arrow == '<=>':
        lb = -1000
        ub = 1000
    #elif string_arrow == '<--':
        # raise NotImplementedError('Backwards reactions are not yet supported')
        # lb = -1000
        # ub = 0
    # todo add feature for people that want to define backwards reactions <-

    return (lb, ub)

#(lb, ub) = get_lb_and_ub_from_formula_string('test --> test33 + ff')


def add_reactions_to_model(
        cobra_model: Model,
        dataframe_metabolites_to_change: pd.DataFrame,
        dataframe_reactions_to_change: pd.DataFrame
):
    """
    creates reaction, add stochiometry and gene rules
    then add reaction to model
    :param cobra_model:
        A cobra compatible model
    :param dataframe_metabolites_to_change:
    :param dataframe_reactions_to_change:
    :return:
    """

    # add metabolites to model
    cobra_model, log_metabolites = add_metabolites_to_model(cobra_model,dataframe_metabolites_to_change)

    start_number_reactions = len(cobra_model.reactions)

    list_all_ensembl_ids_added = [] # for making genes log
    # create reactions
    for idx in range(len(dataframe_reactions_to_change)):
        formula = dataframe_reactions_to_change.loc[idx, 'after_reaction_formula_readable']
        (lb, ub) = get_lb_and_ub_from_formula_string(formula)

        # create reaction
        individual_reaction = Reaction(
            dataframe_reactions_to_change.loc[idx, 'after_reaction_ID'],
            name=dataframe_reactions_to_change.loc[idx, 'after_reaction_name'],
            subsystem=dataframe_reactions_to_change.loc[idx, 'after_reaction_subsystem'],
            lower_bound=lb,
            upper_bound=ub
        )

        # add reaction to model
        cobra_model.add_reactions([individual_reaction])

        # add stochiometry to reaction
        dict_stochiometry = create_dict_stochiometry_by_reaction(cobra_model,formula)
        individual_reaction.add_metabolites(dict_stochiometry)

        # add GPRs to reaction
        GPR_gene_symbol = dataframe_reactions_to_change.loc[idx, 'after_reaction_GPR_readable']
        GPR_ensembl = find_ensemblid_from_gene_symbols(GPR_gene_symbol)
        individual_reaction.gene_reaction_rule = GPR_ensembl

        # add the ENSEMB ids in the GPR_ensembl string to a list
        for match in re.finditer('([ENSG]\w+)', GPR_ensembl):
            list_all_ensembl_ids_added.append(match.group())

        # add ec number to reaction
        individual_reaction.EC_number = ""

    # create dataframe with the ENSEMBL ids that were already in the model and the reactions they are involved in
    log_genes = create_log_for_genes(cobra_model, list_all_ensembl_ids_added)

    end_number_reactions = len(cobra_model.reactions)
    added_number_reactions = start_number_reactions - end_number_reactions
    print(added_number_reactions, 'reactions were added to the model')

    return cobra_model, log_genes, log_metabolites

if __name__ == "__main__":

    main_data_folder = r"C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/model_modifications/versions/version_8"

    metabolites_location = os.path.join(main_data_folder, "metabolites_additions_v7.xlsx")
    reactions_location = os.path.join(main_data_folder, "reactions_additions_v7.xlsx")
    model_location = os.path.join(main_data_folder, "model_v17_with_manual_mods.json")

    dataframe_metabolites_to_change = pd.read_excel(metabolites_location)
    dataframe_reactions_to_change = pd.read_excel(reactions_location)
    cobra_model = load_json_model(model_location)
    #change to load_matlab_model if needed

    output_model, output_log_genes, output_log_metabolites = add_reactions_to_model(cobra_model, dataframe_metabolites_to_change, dataframe_reactions_to_change)

    # save everything
    save_json_model(output_model, os.path.join(main_data_folder, "output/output_model.json"))
    output_log_genes.to_excel(os.path.join(main_data_folder, "output/logs/output_log_genes.xlsx"))
    output_log_metabolites.to_excel(os.path.join(main_data_folder, "output/logs/output_log_metabolites.xlsx"))

