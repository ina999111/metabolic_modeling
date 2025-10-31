# Model Modification Scripts

These are the scripts and documents necessary for doing model modifications.

Metabolite and reaction additions are made through a script by reading Excel files.  
Eventually, I would like modifications and deletions to also be done by this script, but it is not implemented yet.

---

## Pipeline

The pipeline is as follows:

1. From a basic Human1 model (here **v1.17**), do some modifications and deletions using the `manual_modification` script.  
2. Then run the `model_additions` script to add the metabolites and reactions found in **`metabolite_additions.xlsx`** and **`reaction_addition.xlsx`** respectively.  
   This script inputs the model modified in the previous step.

---

## Outputs

The outputs are:

- **Modified model**  
- **`output_log_genes`** — the reactions using each gene in the new model.  
- **`output_log_metabolites`** — for the metabolites added, gives the other metabolites with similar formula.

---
