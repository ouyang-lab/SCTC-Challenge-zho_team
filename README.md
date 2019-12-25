# SCTC-Challenge-zho_team
DREAM Single Cell Transcriptomics Challenge 

The code includes three parts:  
1. Sub challenge
   Input: dge_raw.txt, dge_normalized.txt, binarized_bdtnp.csv, geometry.txt (https://www.synapse.org/#!Synapse:syn16782360).
   Users should register in Synapse (https://www.synapse.org/) using email address and agree to the dataset Conditions for use (https://www.synapse.org/#!Synapse:syn15665609/wiki/583240), then they can download these files
   
   Code: challenge_code_1.R, challenge_code_2.R, challenge_code_3.R, challenge_utils.R. 
   
   Output: 
   - raw.data 
   - normalized.data 
   - geometry
   - insitu.matrix
   - insitu.genes
   - gene60, gene40, gene20
   - dm60, dm40, dm20
   - bin60, bin40, bin20
   - maxc60, maxc40, maxc20
   - sl60, sl40, sl20

2. Post sub challenge - 10 fold cross validation
   Input: folds_train.csv, folds_test.csv
   Code: post_challenge_all.R, newUtils.R
   Output: 
   - binDist60, binDist40, binDist20
   
3. Scores 1, 2, 3 (https://github.com/dream-sctc/Scoring/blob/master/dream_scoring_clean.R)

   Input: output of 1, 2 
   
   Code: estimate_all.R, estimate_util.R 
   
   Output: results in the manuscript 
 
DOI: https://doi.org/10.5281/zenodo.3592532
