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

def read_in_file(filename): #if using a cleaned up ref file as made in make_ref_file
    df=pd.read_csv(filename, sep="\t")
    df_arrays=np.array(df)
    correct_df=df_arrays.transpose()
    indiv_names=df.columns
    clean_input=np.insert(correct_df, 0, indiv_names, axis=1)
    return correct_df
def individual_names(filename):
    df=pd.read_csv(filename, sep="\t")
    indiv_names=df.columns
    return indiv_names
def cleaned_input(filename):
    correct_df=read_in_file(filename)
    indiv_names=individual_names(filename)
    clean_input=np.insert(correct_df, 0, indiv_names, axis=1)
    return clean_input
#result of this is the cleaned up samples, with the correct sample names for each array of genotype calls
def clean_input(filename):
    reference_file=read_in_file(filename)
    reference_positions=positions_bp
    ref_names=individual_names(filename)
    newarrayshape=reference_file.shape
    newdf=np.zeros(list(newarrayshape) + [2])
    for i, element in enumerate(correct_df):
        for j, w in enumerate(element):
            newdf[i][j]=[int(call) for call in w.split(',')]
    cleaned_input_ref=newdf
    return cleaned_input_ref

class reference_panel():
    def __init__(self, filename, header_position):
        self.filename=filename
        self.header_position=header_position
        
    def make_ref_file(self): #for 1kg files, header=95
        df=pd.read_csv(self.filename, header=self.header_position, sep='\t')
        del df['#CHROM']
        dropcol=df.columns[1:8] #removes the non POS columns before the sample columns with the genotype calls
        df=df.drop(dropcol, axis=1)
        self.df=df#at this point we have a VCF file in the dataframe with only the samples and calls, and the position
        return self.df

    def ref_call_positions_bp(self): #store the positions and remove them from the data frame
        df=self.make_ref_file(self.filename, self.header_position)
        self.positions_bp=df['POS']
        return self.positions_bp

    def cleaned_df(self):
        del self.df['POS']
        return self.df

    def reference_panel_input(self):
        reference_file=self.df
        reference_positions=self.positions_bp
        ref_names=clean_inputs(ref_sample_df)
        newarrayshape=correct_df.shape
        newdf=np.zeros(list(newarrayshape) + [2])
        for i, element in enumerate(correct_df):
            for j, w in enumerate(element):
                newdf[i][j]=[int(call) for call in w.split(',')]
        self.reference_panel=newdf


class samples_from_file(filename):
    def make_sample_file(filename, header_position):
        df=pd.read_csv(filename, header=header_position, sep='\t')
        del df['#CHROM'] # assuming VCF #this may need to change depending on file 
        sample_positions_bp=df['POS']
        dropcol=newdf.columns[1:8] #removes the non POS columns before the sample columns with the genotype calls
        df=df.drop(dropcol, axis=1)
        
        self.get_sample_POS=sample_positions_bp
        #at this point we have a VCF file in the dataframe with only the samples and calls, and the position
        return df, sample_positions_bp
    def sample_indiv(sample_file):
        sample_df=clean_inputs(sample_file)
        return sample_df


class windows_phasing():
    def __init__(self, chrom, window_size, reference_panel, samples_input):
        self.chrom=chrom
        self.window_size=window_size
        self.reference_panel=reference_panel
        self.samples_input=samples_input
    def windowing(self):
        """
        RETURNS the window ranges in terms of cM and bp positions 
        input is the chromosome number and the desired window size
        *can build this into a variable size or sliding window definition
        """
        filename_pattern='~/testpy/rupasandbox/hapmap_recombination_rate_hg38/hapmap_recomb_hg38_chr{}.txt'
        self.recomb_hapmap=pd.read_csv(filename_pattern.format(self.chrom), sep=" ")
        x1=0
        x2=25
        window_cM=window_size
        listcM=[x*10 for x in range(x1,x2)]
        windowpoints=[]
        windowpoints_bp=[]
        for element in listcM:
            for i in range(len(self.recomb_hapmap)):
                if self.recomb_hapmap['Genetic_Map(cM)'][i]>element and self.recomb_hapmap['Genetic_Map(cM)'][i-1]<element :
                    windowpoints.append(self.recomb_hapmap['Genetic_Map(cM)'][i])
                    windowpoints_bp.append(self.recomb_hapmap['position'][i])
        windowpoints.append(self.recomb_hapmap['Genetic_Map(cM)'][len(self.recomb_hapmap)-1])
        windowpoints_bp.append(self.recomb_hapmap['position'][len(self.recomb_hapmap)-1])

        windows=[]
        windows_bp=[]
        for i in range(len(windowpoints)-1):
            newpoint=(windowpoints[i],windowpoints[i+1])
            newbppoint=(windowpoints_bp[i], windowpoints_bp[i+1])
            windows.append(newpoint)
            windows_bp.append(newbppoint)
        self.windows=windows
        self.windows_bp=windows_bp

    def window_search(self): # we are searching within the positions that indicate the ()cM window
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
        for (windowstart,windowend) in self.windows_bp: #iterating through each window
            for individual in self.reference_samples: #going individual by individual in the reference panel
                individual_ref=self.reference_samples[individual] #so we can iterate through the genotype calls of each individual easily
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
