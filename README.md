# Phasing

## Basic plans outline
Basic idea is: 
each element of the ordered pair is from either w (parent1)  or x (parent2)
can iterate through each entry s.t. for each y in train_set, 
for sample in train_set:
    for y in sample:
        for w,x in y:
            train model such that it can predict the pair w,x in order from y (unordered pair)
for each SNP, the reference space is the sample space (for simplicity)
so the possibilities are the ordered pairs we see, and the likelihood ~ frequency in the sample set we have
can iterate through the samples to count and --> state space



something like shapeit, where it counts along a "match" with weight ~ frequency of that length 
of 0/1s occurring in that length of ordered pairs (0|1 or 1|0) in the sample, using the homozyg 
0|0 and 1|1 1 to guide the selection of the length.

maybe something like the dmers match or deepu's recombination program? needs to be windowed


### WINDOWS
-STEP 1: WINDOWS
  To Build searchable windows:
 -1) We can first change the datastructure fed in to retain the POS and CHR so we can be aware of the locations by cM of the SNPs.
      a)This means that I'll just move over the changes from the maketestfile into the phasing1.py, so it won't take long
 -2) Then have a cM-POS map fed into the model and figure out how to define cM in the model and keep track of windows by position for each chromosome
 -3) Define the window structure based on the map as 25cM wide, and in terms of where it is centered.
        first window centered at 12.5cM (so it goes from 0-25cM)
        each next window is centered 12.5cM after the previous, so there is an overlap
            2nd window: centered at 25cM (from 12.5cM to 37.5cM)
            3rd window: centered at 37.5cM (from 25cM to 50cM)
            and so on…
 -4) Define a way for the input to be broken into these overlapping windows for searching  

### sorting based on homozygosity
-STEP 2: Sorting subsets of individuals based on homozygosity to reduce search problem
 -1) defining homozygosity within the windows 
 -2) define homozygosity patterns and how to subset
 -3) sorting samples into 'subsets' based on homozygosity matches in windows
    if a sample has the same 001001 homozygosity pattern it goes in subset 1 with all the others that do
    if a sample has an almost same pattern 001000 we allow some mismatch and looseness with recombination/mutation term?

if no matches in the windows, return to step 1, and make window sizes smaller to start, then continue through step 2

### then define search based on windows
-STEP 3: Define search process
 -1) first, catalogue homozygosity for the sample individual and generate the homozygosity pattern for each window
 -2) then break pattern into windows
 -3) define sorting based on shared pattern (similar to Deepu's method, of prioritizing, as explained by Umar)
 -4) search within the subset for best matches at each position (akin to Deepu's recombination method) and 
    phase at heterozygous positions based on best matches overall for both chromosomes (looking at diplotypes)

In terms of the pyro problem of setting up the model: 
treat the subset of the haplotype space (reference) as the prior

Sample from it and see how it fits
then update the output (posterior) based on this

## Backlog/future plans:
- 1)turn windows into _sliding windows_
- 2)error inference - see if we can figure out a way to infer errors in genotype calls and return them as missing data for imputation later
