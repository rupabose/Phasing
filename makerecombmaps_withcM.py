import pandas as pd 
import time
"""
can either use the function to extract individually or run the for loop at the bottom to get all chromosome files
I used the for loop at the bottom in ipython to get all the maps, then zipped them.
"""
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
            cMdist.append(sum(cM_bp[0:i]))
        del df['#Chromosome']
        del df['End']
        df.insert(2,"cM", cMdist)
        pop_recomb_maps[population]=df
        filename_print_pattern="~/testpy/rupasandbox/Chr_maps/Chromosome_{}/recomb_map_chr{}_pop{}"
        df.to_csv(filename_print_pattern.format(chr, chr, pop_slug_new), index=False, sep="\t", line_terminator="\n")
    elapsed_time=start_time - time.time()
    print("chromosome{} took {} minutes".format(chr,elapsed_time))

#can run the below in iPython or comment the above and run the file
for chrom in range(1,23):
    chr=chrom
    start_time=time.time()
    pops=['ACB', 'ASW', 'BEB', 'CDX', 'CEU', 'CHB', 'CHS', 'CLM', 'ESN', 'FIN', 'GBR', 'GIH', 'GWD', 'IBS', 'ITU', 'KHV', 'LWK', 'MSL', 'MXL', 'PEL', 'PJL', 'PUR', 'STU', 'TSI', 'YRI']
    pop_recomb_maps = {}
    recomb_map_filename_pattern = "~/testpy/rupasandbox/Phasing/hg38/{}/{}_recombination_map_hg38_chr_{}.bed" #place chr of choice in last {}
    for pop_slug in pops:
        pop_recomb_maps[pop_slug] = pd.read_csv(recomb_map_filename_pattern.format(pop_slug, pop_slug, chr), sep="\t")
    for population in pop_recomb_maps:
        pop_slug_new=population
        df=pop_recomb_maps[population]
        cM_bp= (df['End']-df['Start'])*df['reco_rate_per_base_per_generation']*100
        cMdist=[]
        for i in range(len(cM_bp)):
            cMdist.append(sum(cM_bp[0:i]))
        del df['#Chromosome']
        del df['End']
        df.insert(2,"cM_cum", cMdist)

        filename_print_pattern="~/testpy/rupasandbox/Chr_maps/Chromosome_{}/recomb_map_chr{}_pop{}"
        df.to_csv(filename_print_pattern.format(chr, chr, pop_slug_new), index=False, sep="\t", line_terminator="\n")
    elapsed_time=(time.time()-start_time)/60
    print("chromosome{} took {} minutes".format(chr,elapsed_time))




