# Model Modifications

These are the scripts and documents necessary for doing model modifications.

Metabolite and reaction additions are made through a script by reading Excel files.  
Eventually, I would like modifications and deletions to also be done by this script, but it is not implemented yet.

---

## Pipeline

The pipeline is as follows:

1. From a basic Human1 model (here **v1.17**), do some modifications and deletions using the `manual_modification` script.  
2. Then run the `model_additions` script to add the metabolites and reactions found in **`metabolite_additions.xlsx`** and **`reaction_addition.xlsx`** respectively.  
   This script inputs the model modified in the previous step.

| script                        | what it does          | inputs                                                                                         | outputs                                                                                                                                                                                                   |
| ----------------------------- | --------------------- | ---------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| manual_model_modifications.py | changes and deletions | basic human gem model (here v1.17)                                                             | model with manual modifications                                                                                                                                                                           |
| model_additions.py            | additions             | model_with_manual_modifications.json<br>metabolites_additions.xlsx<br>reactions_additions.xlsx | output_model<br>logs:<br>\- output_log_genes: the reactions using each gene in the new model. <br>\- output_log_metabolites: for the metabolites added, gives the other metabolites with similar formula. |
