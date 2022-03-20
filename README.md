# Unit Commitment
Python code for unit commitment with thermal units.
Algorithms are based on Power Generation, Operation, and Control, Allen J. Wood, Bruce F. Wollenberg, Gerald B. Sheblé (2013).

The unit commitment problem (UC) in electrical power production is a large family of mathematical optimization problems where the production of a set of electrical generators is coordinated in order to achieve some common target, usually either matching the energy demand at minimum cost or maximizing revenue from electricity production. This is necessary because it is difficult to store electrical energy on a scale comparable with normal consumption; hence, each (substantial) variation in the consumption must be matched by a corresponding variation of the production. 
[Source: Wikipedia](https://en.wikipedia.org/wiki/Unit_commitment_problem_in_electrical_power_production).

## How to use? 
See [Jupyter Notebook example](https://github.com/kypexfly/unit-commitment/blob/master/uc_load_curve.ipynb)
1. Import libraries
2. Set generator cost functions and generation limits in "data"
3. Set load sequence (time & power demand) in "sequence"
4. Run code

## Notes
* ELD code for MATLAB and Python are available [here](https://github.com/kypexfly/economic-load-dispatch)
* If you find any error/problem I would appreciate feedback.
