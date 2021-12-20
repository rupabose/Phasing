# type: ignore
import torch
import torch.nn as nn
from torch.distributions import constraints
import functools
import pyro
import pyro.contrib.examples.polyphonic_data_loader as poly
import pyro.distributions as dist
from pyro import poutine
from pyro.infer import SVI, JitTraceEnum_ELBO, TraceEnum_ELBO, TraceTMC_ELBO
from pyro.infer import Trace_ELBO
from pyro.infer.autoguide import AutoDelta
from pyro.ops.indexing import Vindex
import pyro.optim
from pyro.optim import Adam
from pyro.util import ignore_jit_warnings
from pyro.infer.mcmc import MCMC, NUTS
import pyro.contrib.funsor
import pyroapi
from pyroapi import infer, handlers, ops, optim, pyro, pyro_backend

import os 
from functools import partial
import numpy as np
import seaborn as sns 
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 
import argparse
import logging
from subprocess import check_output
from sys import argv
import sys

def clean_inputs(filename): #input file can be a file with the samples as column headers and the columns being the genotype calls at positions indicated by the rows
    #the input file is rearranged into a form that the phasing model and pyro can work with
    df=pd.read_csv(filename, sep="\t")
    df_arrays=np.array(df)
    correct_df=df_arrays.transpose()
    sample_names=df.columns
    clean_input=np.insert(correct_df, 0, sample_names, axis=1)
    return clean_input
    #result of this is the cleaned up samples, with the correct sample names for each array of genotype calls

class reference_panel():
    def make_ref_file(filename, header_position): #for 1kg files, header=95
        df=pd.read_csv(filename, header=header_position, sep='\t')
        del df['#CHROM']
        dropcol=newdf.columns[1:8] #removes the non POS columns before the sample columns with the genotype calls
        df=df.drop(dropcol, axis=1)
        #at this point we have a VCF file in the dataframe with only the samples and calls, and the position
        return df

    def ref_call_positions_bp(df): #store the positions and remove them from the data frame
        positions_bp=df['POS']
        del df['POS']
        ref_sample_df=df
        return ref_sample_df

    def reference_panel_input(ref_sample_df):
        reference_samples=clean_inputs(ref_sample_df)
        reference_positions=

class positions_for_calls():
    def ref_positions():

    def sample_positions():


def windowing(chr, window_size):
    """
    RETURNS the window ranges in terms of cM and bp positions 
    input is the chromosome number and the desired window size
    *can build this into a variable size or sliding window definition
    """
    filename_pattern='~/testpy/rupasandbox/hapmap_recombination_rate_hg38/hapmap_recomb_hg38_chr{}.txt'
    recomb_hapmap_chr20=pd.read_csv(filename_pattern.format(chr), sep=" ")
    x1=0
    x2=25
    window_cM=window_size
    listcM=[x*10 for x in range(x1,x2)]
    windowpoints=[]
    windowpoints_bp=[]
    for element in listcM:
        for i in range(len(recomb_hapmap_chr20)):
            if recomb_hapmap_chr20['Genetic_Map(cM)'][i]>element and recomb_hapmap_chr20['Genetic_Map(cM)'][i-1]<element :
                windowpoints.append(recomb_hapmap_chr20['Genetic_Map(cM)'][i])
                windowpoints_bp.append(recomb_hapmap_chr20['position'][i])
    windowpoints.append(recomb_hapmap_chr20['Genetic_Map(cM)'][len(recomb_hapmap_chr20)-1])
    windowpoints_bp.append(recomb_hapmap_chr20['position'][len(recomb_hapmap_chr20)-1])

    windows=[]
    windows_bp=[]
    for i in range(len(windowpoints)-1):
        newpoint=(windowpoints[i],windowpoints[i+1])
        newbppoint=(windowpoints_bp[i], windowpoints_bp[i+1])
        windows.append(newpoint)
        windows_bp.append(newbppoint)
    return windows, windows_bp

def window_search(windows_bp,reference_samples): # we are searching within the positions that indicate the ()cM window
    """
    INPUTS: 
    the windows_bp defines the BP positions which bound each ()cM window
    the reference_samples = the reference panel 

    what it does:
    searches through the reference panel within a window and returns the homozygosity pattern for each individual
    if the homozygosity pattern has already been found and logged, then it just adds the individual indicator for the signal
        the individual indicator is the first positin in the array

    we don't have the positions stored in the same array as the genotype calls, so we are actually pulling this from a different set
    in the array with the genotype calls, we actually start from index 1 (index 0 is the sample identifier)
    """
    pop_subsets={} 
    for (windowstart,windowend) in windows_bp: #iterating through each window
        for individual in reference_samples: #going individual by individual in the reference panel
            individual_ref=reference_samples[individual] #so we can iterate through the genotype calls of each individual easily
            individual_window_ref=individual_ref[windowstart:windowend]
            homozyg_sig='' #starting a new homozygous signature as a string, which represents the homozyg pattern of this individual
            for g1,g2 in individual_window_ref:
                if g1==g2: #true if it is homozygous (e.g. 0,0 or 1,1)
                    homozyg_sig+=f"{individual_ref['POS']}|{g1}" #we store position and the call, 0 or 1 for ref/alt
            if pop_subsets.get(homozyg_sig) is None: #if this specific homozyg pattern in this window is not already represented, we add it
                pop_subsets[homozyg_sig]=[individual] #we also add the individual under it
            else:
                pop_subsets[homozyg_sig].append(individual) #if this is already represented in the dict, then we just store that the individual has it too
    return pop_subsets

def windowing_samples_on_reference_panel(): 

class phasing_search_model():
    def __init__(reference_panel, chr, unknown_samples):
        
    def windowing_on_panel():

    def window_matching_samples_to_panel():

    def match_find_v2(pop_subsets, unknown_phase_samples):
        for homozygous_sig in pop_subsets:
            pass
