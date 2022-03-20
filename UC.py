""" UNIT COMMITMENT

CODE BY: Ricardo Chu Zheng, 2022
GITHUB: https://github.com/kypexfly
REFERENCE: Power Generation, Operation, and Control, Allen J. Wood, Bruce F. Wollenberg, Gerald B. Sheblé (2013) """

import numpy as np
import itertools
import pandas as pd
from ELD import iterative_lambda


def unit_commitment(data: np.ndarray, PD: float):
    Pmin = data[:, 3]
    Pmax = data[:, 4]
    Ngen = data.shape[0]

    Pmatrix = np.zeros((2**Ngen, Ngen))
    Cmatrix = np.zeros(2**Ngen)
    unitState = np.array(list(itertools.product([False, True], repeat=Ngen)))*1
    isFeasible = np.zeros(2**Ngen, dtype=int)

    Plimits_state = np.empty((0, 2))
    for i in range(2**Ngen):
        Pmax_state = np.sum(Pmax*unitState[i, :])
        Pmin_state = np.sum(Pmin*unitState[i, :])
        Plimits_state = np.append(
            Plimits_state, [[Pmin_state, Pmax_state]], axis=0)
        if (PD < Pmax_state and PD > Pmin_state):
            isFeasible[i] = 1

    feasible = np.nonzero(isFeasible)[0]

    for i in feasible:
        m = np.nonzero(unitState[i, :])  # ON generator
        k = np.nonzero(1-unitState[i, :])  # OFF generator (deleted)
        modData = np.copy(data)
        modData = np.delete(modData, k, 0)
        P0 = np.zeros(Ngen)
        P0[m] = iterative_lambda(modData, PD)
        cost = np.zeros(Ngen)
        cost[m] = data[m, 0]*P0[m]**2 + data[m, 1]*P0[m] + data[m, 2]
        Pmatrix[i, :] = P0
        Cmatrix[i] = np.sum(cost)

    unit = [f"U{i+1}" for i in range(Ngen)]
    power = [f"P{i+1}" for i in range(Ngen)]

    table = pd.DataFrame(
        np.concatenate([unitState, Plimits_state, isFeasible.reshape(-1,
                       1), Pmatrix, Cmatrix.reshape(-1, 1)], axis=1),
        columns=[*unit, "Pmin", "Pmax", "Feasible", *power, "Cost"])

    table[[*unit, "Feasible"]] = table[[*unit, "Feasible"]].astype("int")

    return unitState, Pmatrix, Cmatrix, table
