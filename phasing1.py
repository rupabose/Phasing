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
# 1)define 25cM wide windows (or maybe based on #SNPs), first one is centered at 12.5cM, and then each successive one is centered 10cM 
    after the previous one's center. window 1 is centered at 12.5, W2 at 22.5, W3 at 32.5 and so on ...
 
"""
#  

#sorting based on homozygosity
""" STEP 2: Sorting subsets of individuals based on homozygosity to reduce search problem
# 1) defining homozygosity within the windows 
# 2) sorting samples into 'subsets' based on homozygosity matches in windows

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

    hidden_dim = 2  # reference or alternative
    #this runs the prior only once to establish the transition probabilities
    """   
    Need to figure out how to write the prior as a function of the possibilities at each position
    s.t. the freq of each pair ~ likelihood of each ordered pair
    so prob(event0,0)=proportion of samples with 0,0 at that position, and so on
    """
    with poutine.mask(mask=include_prior):
        probs_w = pyro.sample(
            "probs_w", dist.Beta(0.5, 0.5).expand([2]).to_event(1)
        )
        # underlying genotype
        probs_x = pyro.sample(
            "probs_x", dist.Beta(0.5,0.5).expand([2]).to_event(1)
        )
        # observed genotype
        probs_y = pyro.sample(
            "probs_y",
            dist.Beta(0.99,0.1).expand([2,2]).to_event(2)
        )
    with pyro.plate("sequences", num_samples, batch_size, dim=-2) as batch:
        lengths = lengths[batch]
        w, x = torch.tensor(0).to(torch.long), torch.tensor(0).to(torch.long)
        for t in pyro.markov(range(sequence_length)):
            with poutine.mask(mask=torch.BoolTensor(t < lengths).unsqueeze(-1)) as unsqueezed_mask:
                w = pyro.sample(
                    "w_{}".format(t),
                    dist.Bernoulli(probs_w[w]),
                    infer={"enumerate": "parallel"},
                ).to(torch.long)
                # print("w_{}".format(t), w)
                x = pyro.sample(
                    "x_{}".format(t),
                    dist.Bernoulli(probs_x[x]),
                    infer={"enumerate": "parallel"},
                ).to(torch.long)
                # print("x_{}".format(t), x)
                y = pyro.sample(
                    "y_{}".format(t),
                    dist.Bernoulli(probs_y[x,w]),
                    obs=sequences[batch, t],
                    )
                # print("y_{}".format(t), y)
                if t <= 2:
                    print("w_{}.shape".format(t), w.shape)
                    print("x_{}.shape".format(t), x.shape)
                    print("y_{}.shape".format(t), y.shape)



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

    # We'll train using MAP Baum-Welch, i.e. MAP estimation while marginalizing
    # out the hidden state x. This is accomplished via an automatic guide that
    # learns point estimates of all of our conditional probability tables,
    # named probs_*.
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

    logging.info("training loss = {}".format(train_loss / num_observations))

    # Finally we evaluate on the test dataset.
    logging.info("-" * 40)
    logging.info(
        "Evaluating on {} test sequences".format(len(data["test"]["sequences"]))
    )
    sequences = data
    lengths = data
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



