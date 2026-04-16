These are scripts made to be used with the outputs of the SWAP (or SWAMP... name TBD) pipeline.
Includes:
- creating a summary dataframe from a folder of filtered.xlsx file. For each sample and task combination: objective, time, reactions_used, etc.
- creating a summary of the reactions in common used for each task combinations (eg. task 1 vs task 2, how many common used reactions ?)
- creating a task and sample specific model. uses the reactions_used column of the output summary dataframe. can be used for mapping in escher.
