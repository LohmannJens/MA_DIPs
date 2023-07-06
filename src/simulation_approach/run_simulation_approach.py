'''

'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


sys.path.insert(0, "..")
from utils import RESULTSPATH, SEGMENTS, COLORS, NUCLEOTIDES, QUANT, N_SAMPLES


def initialise_population(n_gain: int, n_loss: int)-> pd.DataFrame:
    '''
    
    '''
    # create n candidates and assign a random size between to 1 and 100 to it
    gain = np.array([["gain", np.random.randint(1, 100)] for _ in range(n_gain)])
    loss = np.array([["loss", np.random.randint(1, 100)] for _ in range(n_loss)])
    data = np.concatenate((gain, loss))
    df = pd.DataFrame(data=data, columns=["class", "count"])
    df["count"] = df["count"].astype(np.float16)

    return df


def calculate_change(row, prob_gain, prob_loss, rate_increase, rate_decrease):
    '''
    
    '''
    # calculate random number
    thresh = np.random.random()

    if row["class"] == "gain":
        prob = prob_gain
    elif row["class"] == "loss":
        prob = prob_loss

    if prob > thresh:
        rate = rate_increase
    else:
        rate = rate_decrease

    change = round(row["count"] * rate)

    if change < 0:
        change = min(change, -3)

    return change


def run_single_step(pop, prob_gain, prob_loss, rate_increase, rate_decrease)-> pd.DataFrame:
    '''
    
    '''
    
    # calculate change of each candidate
    pop["change"] = pop.apply(calculate_change, args=[prob_gain, prob_loss, rate_increase, rate_decrease], axis=1)
    pop["count"] = pop["count"] + pop["change"]

    # kick out candidates with count <= 0
    pop.drop(pop[pop["count"] <= 0].index, inplace=True)

    return pop


def run_simulation(pop: pd.DataFrame, n: int, prob_gain, prob_loss, rate_increase, rate_decrease)-> pd.DataFrame:
    '''
    
    '''

    for i in range(n):
        print(f"### iteration {i} ###")
        pop = run_single_step(pop, prob_gain, prob_loss, rate_increase, rate_decrease)

        #print(pop[pop["class"] == "loss"])

        # calculate overall number of counts
        # calculate overall number of candidates (number of rows of df)
        print(pop.shape[0])
        print(pop["count"].sum())


if __name__ == "__main__":
    # Define rates and initialise starting population
    steps = 100

    n_gain = 261
    n_loss = 373

    prob_gain = 0.5
    prob_loss = 0.5

    rate_increase = 0.1
    rate_decrease = -0.1

    population = initialise_population(n_gain, n_loss)

    # Run simulation
    run_simulation(population, steps, prob_gain, prob_loss, rate_increase, rate_decrease)
