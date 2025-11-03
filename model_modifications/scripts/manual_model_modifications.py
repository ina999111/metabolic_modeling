"""
Using the cobra functions modified by Jelle to open and save the model, since the original
cobra functions remove some informations (e.g. subsystems, E.C. numbers).
"""

import sys

sys.path.append(r"C:\Users\inapa\PycharmProjects\SWAMP")

from src.cobrapy_fork._cobra import (
    load_matlab_model,
    save_json_model,
)

'''
This script is to be used to modify reactions in the model, while the model addition pipeline does
not support modifications (only additions)

The changes:
1. changing pyruvate kinase reaction (MAR06627) to be forward only 
2. changing the GPR of the reaction PGP-CL pool --> PG-CL pool (MAR00584):
    was empty, add PTPMT1 (ENSG00000110536)
3. changing the GPR of the reaction phosphatidate-CL pool --> CDP-DAG-CL pool (MAR00581):
    was CDS1 or CDS2, changed to TAMM41 (ENSG00000144559)
4. Adding GPR to MAR11311 (was empty)
    ENSG00000081248 or ENSG00000151067 or ENSG00000157388 or ENSG00000102001
5. Removing MAR07629 
'''

model_v17 = load_matlab_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/model_Human-GEM_COBRA version 17.mat")

# 1. MAR06627 to be forward only
model_v17.reactions.get_by_id('MAR06627').bounds

model_v17.reactions.get_by_id('MAR06627').bounds = (0, 1000)

model_v17.reactions.get_by_id('MAR06627').bounds
model_v17.reactions.get_by_id('MAR06627').reaction

# 2. add GPR to MAR00584

model_v17.reactions.get_by_id('MAR00584').gpr # no GPR
model_v17.genes.get_by_id('ENSG00000110536') # in model

model_v17.reactions.get_by_id('MAR00584').gene_reaction_rule = 'ENSG00000110536'

model_v17.reactions.get_by_id('MAR00584').gpr # worked

# 3. changing GPR of MAR00581 to ENSG00000144559

model_v17.reactions.get_by_id('MAR00581').gpr # cobra.core.gene.GPR('ENSG00000101290 or ENSG00000163624')
model_v17.genes.get_by_id('ENSG00000144559') # not in model

model_v17.reactions.get_by_id('MAR00581').gene_reaction_rule = 'ENSG00000144559'

model_v17.reactions.get_by_id('MAR00581').gpr # worked
model_v17.genes.get_by_id('ENSG00000144559') # worked

# 4. adding GPR to MAR11311

model_v17.reactions.get_by_id('MAR11311').gpr # empty
model_v17.genes.get_by_id('ENSG00000081248')
model_v17.genes.get_by_id('ENSG00000151067')
model_v17.genes.get_by_id('ENSG00000157388')
model_v17.genes.get_by_id('ENSG00000102001')
# none of the genes are in the model

model_v17.reactions.get_by_id('MAR11311').gene_reaction_rule = 'ENSG00000081248 or ENSG00000151067 or ENSG00000157388 or ENSG00000102001'

model_v17.reactions.get_by_id('MAR11311').gpr # worked
model_v17.genes.get_by_id('ENSG00000081248')
model_v17.genes.get_by_id('ENSG00000151067')
model_v17.genes.get_by_id('ENSG00000157388')
model_v17.genes.get_by_id('ENSG00000102001')
# all genes added

# 5. Removing MAR07629

model_v17.reactions.get_by_id('MAR07629').reaction

model_v17.remove_reactions([model_v17.reactions.get_by_id('MAR07629')])

model_v17.reactions.get_by_id('MAR07629').reaction

## save manual in a folder where it is going to be used for the model addition script
save_json_model(model_v17, r"C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/model_modifications/versions/version_8/model_v17_with_manual_mods.json")

