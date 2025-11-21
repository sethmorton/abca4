can use a parquet file obtained from CADD or any other ranking tool for mis-sense variants. 

this repo uses a from scratch implementation of mis-sense (one unit of alteration) variants using several AI models. refer to the notebooks directory to learn more about this. 

we then add more scores using LLMs / sequence and virtual cell models (in upcoming versions) to add a sense of which variant should be tested in the lab first given an objective like "finding the causal variant for a disease X' 

The updated ranking of the variants using a weighted formulation of all the scores gives the final set of prioritised variants to test. 