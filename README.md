# Carcinoma-in-situ-signature
Code in current repository can be used to regenerate analysis and figures from the study "Gene signatures characterizing driver mutations in lung squamous carcinoma are predictive of the progression of pre-cancer lesions".

In particular, gene-signature-calculation.R contains the steps of deriving gene signatures from raw input
All data can be publicly downloaded based on GEO accession ID.

#### [1] Identify genomic event highly occured in LUSC 
#### [2] Define weighted profiles 
#### [3] Obtain sample-specific enrichment score from BASE algorithm 
#### [4] Downstream analysis: Association with survival 

Figure 2-5.R contains detailed visualization step 
#### Figure 2: Gene signatures for driver genomic aberrations are predictive of patient prognosis in lung squamous cell carcinoma.
#### Figure 3: Gene signatures defined based on LUSC data can be applied to precancer gene expression data to characterize the relevant pathway activity changes across different precancer stages.
#### Figure 4: Driver gene signatures are predictive of the progression of patients with lung carcinoma in situ.
#### Figure 5: Association between gene signature and immune cells show pattern along cancer progression. 

