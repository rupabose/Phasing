import numpy as np
import pandas as pd

pops=['ACB', 'ASW', 'BEB', 'CDX', 'CEU', 'CHB', 'CHS', 'CLM', 'ESN', 'FIN', 'GBR', 'GIH', 'GWD', 'IBS', 'ITU', 'KHV', 'LWK', 'MSL', 'MXL', 'PEL', 'PJL', 'PUR', 'STU', 'TSI', 'YRI']
# pops=['ACB']


# for chrom in range(1, 22):
for chrom in [17]:

    pop_recomb_maps = {}
    recomb_map_filename_pattern = "~/testpy/rupasandbox/Phasing/hg38/{}/{}_recombination_map_hg38_chr_{}.bed"

    for pop_slug in pops:
        pop_recomb_maps[pop_slug] = pd.read_csv(recomb_map_filename_pattern.format(pop_slug, pop_slug, chrom), sep="\t")

    pop_recomb_map_max_region_cM     = {pop: 0 for pop in pop_recomb_maps.keys()}
    pop_recomb_map_max_recomb_region = {pop: 0 for pop in pop_recomb_maps.keys()}
    pop_recomb_map_max_region_BP     = {pop: 0 for pop in pop_recomb_maps.keys()}


    for pop_slug in pops:
        df = pop_recomb_maps[pop_slug]
        pop_recomb_map_max_region_cM[pop_slug]     = max(  (df['End']-df['Start'])*df['reco_rate_per_base_per_generation']*100  )
        pop_recomb_map_max_region_BP[pop_slug]     = max(  (df['End']-df['Start'])  )
        pop_recomb_map_max_recomb_region[pop_slug] = max( df['reco_rate_per_base_per_generation'] )


    # print(f"""Chromosome #{chrom}
    # max region in cM: {max(pop_recomb_map_max_region_cM.values())}
    # max recomb rate : {max(pop_recomb_map_max_recomb_region.values())}
    # max BP          : {max(pop_recomb_map_max_region_BP.values())}""")

    print(f"""Chromosome #{chrom}""")
    for pop_slug in pops:
            print(f"""population {pop_slug}
    max region in cM: {pop_recomb_map_max_region_cM[pop_slug]}
    max recomb rate : {pop_recomb_map_max_recomb_region[pop_slug]}
    max BP          : {pop_recomb_map_max_region_BP[pop_slug]}""")
