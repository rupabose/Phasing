import pandas as pd 
import time
pops=['ACB', 'ASW', 'BEB', 'CDX', 'CEU', 'CHB', 'CHS', 'CLM', 'ESN', 'FIN', 'GBR', 'GIH', 'GWD', 'IBS', 'ITU', 'KHV', 'LWK', 'MSL', 'MXL', 'PEL', 'PJL', 'PUR', 'STU', 'TSI', 'YRI']
pop_recomb_maps = {}
recomb_map_filename_pattern = "~/testpy/rupasandbox/Phasing/hg38/{}/{}_recombination_map_hg38_chr_{}.bed" #place chr of choice in last {}

def extract_cM(chr):
    pops=['ACB', 'ASW', 'BEB', 'CDX', 'CEU', 'CHB', 'CHS', 'CLM', 'ESN', 'FIN', 'GBR', 'GIH', 'GWD', 'IBS', 'ITU', 'KHV', 'LWK', 'MSL', 'MXL', 'PEL', 'PJL', 'PUR', 'STU', 'TSI', 'YRI']
    pop_recomb_maps = {}
    recomb_map_filename_pattern = "~/testpy/rupasandbox/Phasing/hg38/{}/{}_recombination_map_hg38_chr_{}.bed" #place chr of choice in last {}
    
    for pop_slug in pops:
        pop_recomb_maps[pop_slug] = pd.read_csv(recomb_map_filename_pattern.format(pop_slug, pop_slug, chr), sep="\t")
    
    start_time=time.time()
    for population in pop_recomb_maps:
        df=pop_recomb_maps[population]
        for i in range(len(df)):
            start_bp=df['Start']
            end_bp=df['End']
            recrate=df['reco_rate_per_base_per_generation']
            cMdist=(end_bp-start_bp)*recrate*100
        del df['#Chromosome']
        del df['End']
        del df['reco_rate_per_base_per_generation']
        df.insert(1,"cM", cMdist)
        pop_recomb_maps[population]=df
        df.to_csv("~/testpy/rupasandbox/recomb_map_chr{chr}_pop{population}", index=False, sep="\t", line_terminator="\n")
    elapsed_time=start_time - time.time()
    print(elapsed_time)
    

pops=['ACB', 'ASW', 'BEB', 'CDX', 'CEU', 'CHB', 'CHS', 'CLM', 'ESN', 'FIN', 'GBR', 'GIH', 'GWD', 'IBS', 'ITU', 'KHV', 'LWK', 'MSL', 'MXL', 'PEL', 'PJL', 'PUR', 'STU', 'TSI', 'YRI']
pop_recomb_maps = {}
recomb_map_filename_pattern = "~/testpy/rupasandbox/Phasing/hg38/{}/{}_recombination_map_hg38_chr_{}.bed" #place chr of choice in last {}
chr=1
for pop_slug in pops:
    pop_recomb_maps[pop_slug] = pd.read_csv(recomb_map_filename_pattern.format(pop_slug, pop_slug, chr), sep="\t")

start_time=time.time()
for population in pop_recomb_maps:
    df=pop_recomb_maps[population]
    for i in range(len(df)):
        start_bp=df['Start']
        end_bp=df['End']
        recrate=df['reco_rate_per_base_per_generation']
        cMdist=(end_bp-start_bp)*recrate*100
    del df['#Chromosome']
    del df['End']
    del df['reco_rate_per_base_per_generation']
    df.insert(1,"cM", cMdist)
    pop_recomb_maps[population]=df
    df.to_csv("~/testpy/rupasandbox/recomb_map_chr{chr}_pop{population}", index=False, sep="\t", line_terminator="\n")