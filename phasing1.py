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



logging.basicConfig(format="%(relativeCreated) 9d %(message)s", level=logging.DEBUG)

# Add another handler for logging debugging events (e.g. for profiling)
# in a separate stream that can be captured.
log = logging.getLogger()
debug_handler = logging.StreamHandler(sys.stdout)
debug_handler.setLevel(logging.DEBUG)
debug_handler.addFilter(filter=lambda record: record.levelno <= logging.DEBUG)
log.addHandler(debug_handler)




#let's build a super tiny example to work with
df=pd.read_csv('~/testpy/rupasandbox/phasing_test_file_pairs.csv', sep="\t")
df_arrays=np.array(df)
correct_df=df_arrays.transpose()
correct_df.shape #to check

"""
need to enter in each pair as a set instead of as a string of o or 1

"""
newarrayshape=correct_df.shape
newdf=np.zeros(list(newarrayshape) + [2])
for i, element in enumerate(correct_df):
    for j, w in enumerate(element):
        newdf[i][j]=[int(call) for call in w.split(',')]


train_set=newdf
data=train_set
#testing df has to somehow encode the unordered pairs rather than the ordered pairs
#can maybe feed in ordered pairs and have it check at each point if the ordering makes sense
#not sure how? 

#building the model
"""

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
"""

#WINDOWS
""" STEP 1: WINDOWS
To Build searchable windows:
 1) We can first change the datastructure fed in to retain the POS and CHR so we can be aware of the locations by cM of the SNPs.
      a)This means that I'll just move over the changes from the maketestfile into the phasing1.py, so it won't take long
 2) Then have a cM-POS map fed into the model and figure out how to define cM in the model and keep track of windows by position for each chromosome
 3) Define the window structure based on the map as 25cM wide, and in terms of where it is centered.
        first window centered at 12.5cM (so it goes from 0-25cM)
        each next window is centered 12.5cM after the previous, so there is an overlap
            2nd window: centered at 25cM (from 12.5cM to 37.5cM)
            3rd window: centered at 37.5cM (from 25cM to 50cM)
            and so on…
 4) Define a way for the input to be broken into these overlapping windows for searching
 
"""
#   input a centimorgan map for chr20
#input maps for each of the populations:
pops=['ACB', 'ASW', 'BEB', 'CDX', 'CEU', 'CHB', 'CHS', 'CLM', 'ESN', 'FIN', 'GBR', 'GIH', 'GWD', 'IBS', 'ITU', 'KHV', 'LWK', 'MSL', 'MXL', 'PEL', 'PJL', 'PUR', 'STU', 'TSI', 'YRI']

pop_recomb_maps = {}
recomb_map_filename_pattern = "~/testpy/rupasandbox/Phasing/hg38/{}/{}_recombination_map_hg38_chr_20.bed"

for pop_slug in pops:
    pop_recomb_maps[pop_slug] = pd.read_csv(recomb_map_filename_pattern.format(pop_slug, pop_slug), sep="\t")
#for the maps, convert recombination rate to centimorgans using simple arithmetic
#((end-start)*recombination rate per bp * 100)=cM 
for population in pop_recomb_maps:
    df=pop_recomb_maps[population]
    for i in range(len(df)):
        start_bp=df['Start']
        end_bp=df['End']
        recrate=df['reco_rate_per_base_per_generation']
        cMdist=(end_bp-start_bp)*recrate*100
    df.insert(4,"cM", cMdist)
    pop_recomb_maps[population]=df

#above part works!    

"""
To build the windows, let's first make an actual map of where every 0.5cM  falls using the dataframe
for each population. 
    To do this, we'll make a sum_cM column which will only have an output every 0.5cM
    Then we'll add the column throughout the dictionary, and take the average positions where 0.5, 1, 1.5 etc occur
    And output that as our window map.
"""
sum=0
x1=0
x2=1000000000
list=[x*0.5 for x in range(2*x1,2*x2 +1)]

for population in pop_recomb_maps:
    df=pop_recomb_maps[population]
    for i in range(len(df)):
        sum=sum+df['cM']
        if sum in list:
            df.insert(5,  sum)

window_bps_map= {}
for population in pop_recomb_maps:
    window_bps_map[population]=
    for i in range(len(df)):
        if df['sum'] in list:
            


#1cM is roughly 100,000 base pairs
#25cM is roughly 2.5million base pairs

#reformat into cM instead of recombination rate per base pair per generation
#
#   define 25cM windows starting at 0, 12.5, 25, 37.5 ...
#   link the windows to chr/bp positions
#   organize SNPs in input into the windows   


#sorting based on homozygosity
""" STEP 2: Sorting subsets of individuals based on homozygosity to reduce search problem
# 1) defining homozygosity within the windows 
# 2) define homozygosity patterns and how to subset
# 3) sorting samples into 'subsets' based on homozygosity matches in windows
    if a sample has the same 001001 homozygosity pattern it goes in subset 1 with all the others that do
    if a sample has an almost same pattern 001000 we allow some mismatch and looseness with recombination/mutation term?

"""

"""
if no matches in the windows, return to step 1, and make window sizes smaller to start, then continue through step 2
"""

#then define search based on windows
""" STEP 3: Define search process
#1) first, catalogue homozygosity for the sample individual and generate the homozygosity pattern for each window
#2) then break pattern into windows
#3) define sorting based on shared pattern (similar to Deepu's method, of prioritizing, as explained by Umar)
#4) search within the subset for best matches at each position (akin to Deepu's recombination method) and 
    phase at heterozygous positions based on best matches overall for both chromosomes (looking at diplotypes)
"""

"""
In terms of the pyro problem of setting up the model: 
treat the subset of the haplotype space (reference) as the prior

Sample from it and see how it fits
then update the output (posterior) based on this

"""

def model_3(sequences, lengths, args, batch_size=None, include_prior=True):
    with ignore_jit_warnings():
        num_samples, sequence_length, data_dim = map(int, np.shape(data))
        assert lengths.shape == (num_samples,)
        assert lengths.max() <= sequence_length
        # lengths=np.zeros(num_samples)
        # for i in range(num_samples):
        #     lengths[i]=sequence_length #change it to np.size(data[1])/2

    hidden_dim = int(2)  # reference or alternative
    #this runs the prior only once to establish the transition probabilities
    """   
    Need to figure out how to write the prior as a function of the possibilities at each position
    s.t. the freq of each pair ~ likelihood of each ordered pair
    so prob(event0,0)=proportion of samples with 0,0 at that position, and so on
    """
    with poutine.mask(mask=include_prior):
        probs_w = pyro.sample(
            "probs_w", dist.Dirichlet(0.5 * torch.eye(hidden_dim)+0.5).expand([2]).to_event(1) #instead of 0.5, 0.5, we assigned to the hidden_dimdim
        )
        # underlying genotype
        probs_x = pyro.sample(
            "probs_x", dist.Dirichlet(0.5 * torch.eye(hidden_dim)+0.5).expand([2]).to_event(1)
        )
        # observed genotype
        probs_y = pyro.sample(
            "probs_y",
            dist.Beta(0.999,0.001).expand([hidden_dim, hidden_dim, data_dim]).to_event(3) #instead of 2,2,2 added hidden_dim, hidden_dim, hidden_dimd
        )
    calls_plate = pyro.plate("calls", data_dim, dim=-1)
    with pyro.plate("sequences", num_samples, batch_size, dim=-2) as batch:
        lengths = lengths[batch]
        w, x = torch.tensor(0).to(torch.long), torch.tensor(0).to(torch.long)
        for t in pyro.markov(range(sequence_length)):
            with poutine.mask(mask=torch.BoolTensor(t < lengths).unsqueeze(-1)) as unsqueezed_mask:
                w = pyro.sample(
                    "w_{}".format(t),
                    dist.Categorical(probs_w[w]), #instead of Bernoulli
                    infer={"enumerate": "parallel"},
                ).to(torch.long)
                # print("w_{}".format(t), w)
                x = pyro.sample(
                    "x_{}".format(t),
                    dist.Categorical(probs_x[x]),
                    infer={"enumerate": "parallel"},
                ).to(torch.long)
                # print("x_{}".format(t), x)
                with calls_plate as calls:
                    y = pyro.sample(
                        "y_{}".format(t),
                        dist.Bernoulli(probs_y[x,w, calls]), #added dep on calls
                        obs=sequences[batch, t],
                        )
                # print("y_{}".format(t), y)
                # if t <= 2:
                #     print("w_{}.shape".format(t), w.shape)
                #     print("x_{}.shape".format(t), x.shape)
                #     print("y_{}.shape".format(t), y.shape)



models = {
    name[len("model_") :]: model
    for name, model in globals().items()
    if name.startswith("model_")
}


def main(args):
    if args.cuda:
        torch.set_default_tensor_type("torch.cuda.FloatTensor")

    logging.info("Loading data")
    train = train_set
    valid=train_set
    test=train_set

    # exit(0)
    data= {
        "train": {
            "sequence_lengths": torch.tensor([train.shape[1] for i in train]).to(torch.long),
            "sequences": torch.tensor(train).to(torch.float32),
        },
        "valid": {
            "sequence_lengths": torch.tensor([valid.shape[1] for i in valid]).to(torch.long),
            "sequences": torch.tensor(valid).to(torch.float32),
        },
        "test": {
            "sequence_lengths": torch.tensor([test.shape[1] for i in test]).to(torch.long),
            "sequences": torch.tensor(test).to(torch.float32),
        },
    }



    assert len(data["train"]["sequence_lengths"]) == len(data["train"]["sequences"])



    logging.info("-" * 40)
    model = models[args.model]
    logging.info(
        "Training {} on {} sequences".format(
            model.__name__, len(data["train"]["sequences"])
        )
    )
    sequences = data["train"]["sequences"]
    lengths = data["train"]["sequence_lengths"]
   

    
    if args.truncate:
        lengths = lengths.clamp(max=args.truncate)
        sequences = sequences[:, : args.truncate]
    num_observations = float(lengths.sum())
    pyro.set_rng_seed(args.seed)
    pyro.clear_param_store()

    # AutoDelta is an automatic guide that learns the point estimates of the conditional probs we named probs_{}. 
    # We'll train using MAP Baum-Welch, which marginalizes on hidden states (w and x).
    # 
    guide = AutoDelta(
         handlers.block(model, expose_fn=lambda msg: msg["name"].startswith("probs_"))
     )
    #guide=AutoDelta(poutine.block((model)))

    # To help debug our tensor shapes, let's print the shape of each site's
    # distribution, value, and log_prob tensor. Note this information is
    # automatically printed on most errors inside SVI.
    if args.print_shapes:
        first_available_dim = 2
        guide_trace = handlers.trace(guide).get_trace(
            sequences, lengths, args=args, batch_size=args.batch_size
        )
        model_trace = handlers.trace(
            handlers.replay(handlers.enum(model, first_available_dim), guide_trace)
        ).get_trace(sequences, lengths, args=args, batch_size=args.batch_size)
        logging.info(model_trace.format_shapes())

    # Bind non-PyTorch parameters to make these functions jittable.
    model = functools.partial(model, args=args)
    guide = functools.partial(guide, args=args)

    # Enumeration requires a TraceEnum elbo and declaring the max_plate_nesting.
    # All of our models have two plates: "data" and "tones".
    optimizer = optim.Adam({"lr": args.learning_rate})
    if args.tmc:
        if args.jit and not args.funsor:
            raise NotImplementedError("jit support not yet added for TraceTMC_ELBO")
        Elbo = infer.JitTraceTMC_ELBO if args.jit else infer.TraceTMC_ELBO
        elbo = Elbo(max_plate_nesting=2)
        tmc_model = handlers.infer_config(
            model,
            lambda msg: {"num_samples": args.tmc_num_samples, "expand": False}
            if msg["infer"].get("enumerate", None) == "parallel"
            else {},
        )  # noqa: E501
        svi = infer.SVI(tmc_model, guide, optimizer, elbo)
    else:
        Elbo = infer.JitTraceEnum_ELBO if args.jit else infer.TraceEnum_ELBO
        # max_plate_nesting = 2
        elbo = Elbo(
            # max_plate_nesting=max_plate_nesting,
            strict_enumeration_warning=True,
            jit_options={"time_compilation": args.time_compilation},
        )
        svi = infer.SVI(model, guide, optimizer, elbo)
    

    # We'll train on small minibatches.
    logging.info("Step\tLoss")
    for step in range(args.num_steps):
        loss = svi.step(sequences, lengths, batch_size=args.batch_size) #inputs here are causing it to fail
        #look into model and guide
        #those are likely suspects causing it to fail
        #check lengths as well

        logging.info("{: >5d}\t{}".format(step, loss / num_observations))

    if args.jit and args.time_compilation:
        logging.debug(
            "time to compile: {} s.".format(elbo._differentiable_loss.compile_time)
        )

    # We evaluate on the entire training dataset,
    # excluding the prior term so our results are comparable across models.
    train_loss = elbo.loss(model, guide, sequences, lengths, include_prior=False)
    logging.info("training loss = {}".format(train_loss / num_observations))

    # Finally we evaluate on the test dataset.
    logging.info("-" * 40)
    logging.info(
        "Evaluating on {} test sequences".format(len(data["test"]["sequences"]))
    )
    sequences = data["test"]["sequences"]
    lengths = data["test"]["sequence_lengths"]
    if args.truncate:
        lengths = lengths.clamp(max=args.truncate)
    num_observations = float(lengths.sum())

    # note that since we removed unseen notes above (to make the problem a bit easier and for
    # numerical stability) this test loss may not be directly comparable to numbers
    # reported on this dataset elsewhere.
    test_loss = elbo.loss(
        model,
        guide,
        sequences,
        lengths,
        batch_size=sequences.shape[0],
        include_prior=False,
    )
    logging.info("test loss = {}".format(test_loss / num_observations))

    # We expect models with higher capacity to perform better,
    # but eventually overfit to the training set.
    capacity = sum(
        value.reshape(-1).size(0) for value in pyro.get_param_store().values()
    )
    logging.info("model_{} capacity = {} parameters".format(args.model, capacity))


if __name__ == "__main__":
    assert pyro.__version__.startswith("1.7.0")
    parser = argparse.ArgumentParser(
        description="MAP Baum-Welch learning Bach Chorales"
    )
    parser.add_argument(
        "-m",
        "--model",
        default="3",
        type=str,
        help="one of: {}".format(", ".join(sorted(models.keys()))),
    )
    parser.add_argument("-n", "--num-steps", default=50, type=int)
    parser.add_argument("-b", "--batch-size", default=20, type=int)
    parser.add_argument("-d", "--hidden-dim", default=16, type=int)
    parser.add_argument("-nn", "--nn-dim", default=48, type=int)
    parser.add_argument("-nc", "--nn-channels", default=2, type=int)
    parser.add_argument("-lr", "--learning-rate", default=0.05, type=float)
    parser.add_argument("-t", "--truncate", type=int)
    parser.add_argument("-p", "--print-shapes", action="store_true")
    parser.add_argument("--seed", default=0, type=int)
    parser.add_argument("--cuda", action="store_true")
    parser.add_argument("--jit", action="store_true")
    parser.add_argument("--time-compilation", action="store_true")
    parser.add_argument("-rp", "--raftery-parameterization", action="store_true")
    parser.add_argument(
        "--tmc",
        action="store_true",
        help="Use Tensor Monte Carlo instead of exact enumeration "
        "to estimate the marginal likelihood. You probably don't want to do this, "
        "except to see that TMC makes Monte Carlo gradient estimation feasible "
        "even with very large numbers of non-reparametrized variables.",
    )
    parser.add_argument("--tmc-num-samples", default=10, type=int)
    parser.add_argument("--funsor", action="store_true")
    args = parser.parse_args()

    PYRO_BACKEND = "pyro"

    with pyro_backend(PYRO_BACKEND):
        main(args)



"""
At the end/updates for later: 
error inference - see if we can figure out a way to infer errors in genotype calls 
and return them as missing data for imputation later

"""