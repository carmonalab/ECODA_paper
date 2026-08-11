# Gurobi licence expired when running a model on Yggdrasil

- Source: https://hpc-community.unige.ch/t/3854

- Created: 2025-03-11T09:53:53.091Z

- Tags: yggdrasil

- Posts: 11

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Giacomo.Rubino (2025-03-11T09:53:53.132Z)

## Primary informations
Username: rubino
Cluster: yggdrasil
## Description
I am trying to run a data intensive, energy optimization model on Yggdrasil (sbatch abc.sh). This model uses the solver gurobi for optimization. Till December 2024 I run such models without problems. I tried now March 2025 and I get the error : gurobipy.GurobiError: Web license service only available for container environments
## Steps to Reproduce
I can run a simple version of this model, on Yggdrasil or baobab, but I get the same message error in the output.
## Expected Result
I expect the model to solve the optimization problem smoothly like before (In the mean time, I did not change anything)
## Actual Result
Model run crashed with the following error message:
```
Traceback (most recent call last):
  File "<stdin>", line 5, in <module>
  File "/home/users/r/rubino/baobab_python_env/lib/python3.9/site-packages/pyomo/solvers/plugins/solvers/GUROBI_RUN.py", line 66, in gurobi_run
    model = read(model_file)
  File "src/gurobipy/gurobi.pxi", line 3571, in gurobipy.read
  File "src/gurobipy/gurobi.pxi", line 80, in gurobipy.gurobi.read
  File "src/gurobipy/gurobi.pxi", line 32, in gurobipy.gurobi._getdefaultenv
  File "src/gurobipy/env.pxi", line 62, in gurobipy.Env.__init__
gurobipy.GurobiError: Web license service only available for container environments
ERROR: Solver (gurobi) returned non-zero return code (1)
ERROR: See the solver log above for diagnostic information.
Traceback (most recent call last):
  File "/home/users/r/rubino/now_old_SURE_2035_new_model_version/CH_MUN_EXPANSE_V4_2035_new/cluster_original_168h/control/run_Swiss_EXPANSE_2035_MinCost.py", line 147, in <module>
    model,results = solve_model_mincost(model)
  File "/home/users/r/rubino/now_old_SURE_2035_new_model_version/CH_MUN_EXPANSE_V4_2035_new/cluster_original_168h/EXPANSE/solve_model.py", line 52, in solve_model_mincost
    results = SolverFactory(opts['solver_name']).solve(model,
  File "/home/users/r/rubino/baobab_python_env/lib/python3.9/site-packages/pyomo/opt/base/solvers.py", line 627, in solve
    raise ApplicationError("Solver (%s) did not exit normally" % self.name)
pyomo.common.errors.ApplicationError: Solver (gurobi) did not exit normally
```


## Post 2 by @Adrien.Albert (2025-03-11T14:13:45.749Z)

Duplicate thread:
Gurobi license expired since 14 Feb 2025[Gurobi license expired since 14 Feb 2025](https://hpc-community.unige.ch/t/gurobi-license-expired-since-14-feb-2025/3829) HPC issues[HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14)
> Dear HPC team, 
Hope you are doing well! 
Please find below an issue regarding the Gurobi license. 
Thank you very much in advance for your help! 
Best, 
Xin 
Primary informations
Username: $wenx 
Cluster: $Baobab 
Description
It seems like Gurobi license has expired on Baobab since 14 Feb 2025. 
e.g. job submitted: 15144577, 15144581 
Actual Result
Gurobi shell based on Python 3.9.6 can be launched with command gurobi.sh 
Gurobi Python Interface can be loaded in Python 3.9.6 with ‘import gurobi…


## Post 5 by @Adrien.Albert (2025-03-11T14:21:05.787Z)

Sorry it seems the issue is not caused by the license we have.
Could you give more information, how to reproduce the issue ?


## Post 6 by @Giacomo.Rubino (2025-03-17T10:32:32.095Z)

Dear HPC Team,
Thanks for your reply.
I guess that to reproduce the issue one way could be to run an optimization model that uses our version of Gurobi as a solver.
We can select different solver options to use: gurobi-default, gurobi-stability, gurobi-highspeed etc but the problem seems to occur for any of these options selected.
The error message from the model points us at:
Traceback (most recent call last):
File “”, line 5, in
File “/home/users/r/rubino/baobab_python_env/lib/python3.9/site-packages/pyomo/solvers/plugins/solvers/GUROBI_RUN.py”, line 62, in gurobi_run
model = read(model_file)
File “src/gurobipy/gurobi.pxi”, line 3571, in gurobipy.read
File “src/gurobipy/gurobi.pxi”, line 80, in gurobipy.gurobi.read
File “src/gurobipy/gurobi.pxi”, line 32, in gurobipy.gurobi._getdefaultenv
File “src/gurobipy/env.pxi”, line 62, in gurobipy.Env.init
gurobipy.GurobiError: Web license service only available for container environments
Do you know what “container environments” mean here?


## Post 7 by @Adrien.Albert (2025-03-17T14:58:48.877Z)

Hi Giacomo,
I do not know how Gurobi works, so I ask to chatGPT to create a test:
```
(yggdrasil)-[alberta@login1 gurobi]$ cat test.py 
from gurobipy import Model, GRB

# Create the model
model = Model("Test_Gurobi")

# Add variables
x = model.addVar(name="x", vtype=GRB.CONTINUOUS, lb=0)
y = model.addVar(name="y", vtype=GRB.CONTINUOUS, lb=0)

# Set the objective function (maximize)
model.setObjective(3 * x + 2 * y, GRB.MAXIMIZE)

# Add constraints
model.addConstr(x + 2 * y <= 4, "c1")
model.addConstr(4 * x + y <= 5, "c2")

# Optimize the model
model.optimize()

# Display the results
if model.status == GRB.OPTIMAL:
    print(f"Optimal solution found: x = {x.X}, y = {y.X}")
    print(f"Optimal objective value: {model.objVal}")
else:
    print("No optimal solution found.")
```
Load the modules
```
(yggdrasil)-[alberta@login1 gurobi]$ ml GCCcore/12.3.0 Gurobi/11.0.2
Gurobi shell based on Python 3.11.3 can be launched with command `gurobi.sh`
Gurobi Python Interface can be loaded in Python 3.11.3 with 'import gurobipy'
```
Run via slurm
```
(yggdrasil)-[alberta@login1 gurobi]$ srun python test.py
srun: job 39174116 queued and waiting for resources
srun: job 39174116 has been allocated resources
Set parameter WLSAccessID
Set parameter WLSSecret
Set parameter LicenseID to value 2626997
Academic license 2626997 - for non-commercial use only - registered to hp___@unige.ch
Gurobi Optimizer version 11.0.2 build v11.0.2rc0 (linux64 - "Rocky Linux 9.5 (Blue Onyx)")

CPU model: Intel(R) Xeon(R) Gold 6240 CPU @ 2.60GHz, instruction set [SSE2|AVX|AVX2|AVX512]
Thread count: 36 physical cores, 36 logical processors, using up to 32 threads

Academic license 2626997 - for non-commercial use only - registered to hp___@unige.ch
Optimize a model with 2 rows, 2 columns and 4 nonzeros
Model fingerprint: 0x1a79f54a
Coefficient statistics:
  Matrix range     [1e+00, 4e+00]
  Objective range  [2e+00, 3e+00]
  Bounds range     [0e+00, 0e+00]
  RHS range        [4e+00, 5e+00]
Presolve time: 0.01s
Presolved: 2 rows, 2 columns, 4 nonzeros

Iteration    Objective       Primal Inf.    Dual Inf.      Time
       0    3.5000000e+30   2.750000e+30   3.500000e+00      0s
       2    5.7142857e+00   0.000000e+00   0.000000e+00      0s

Solved in 2 iterations and 0.01 seconds (0.00 work units)
Optimal objective  5.714285714e+00
Optimal solution found: x = 0.8571428571428572, y = 1.5714285714285714
Optimal objective value: 5.714285714285714
```
It seems to work as long as chatgpt doesn’t give me anything wrong.
For you specific issue please provide the procedure steps by steps to reproduce it


## Post 8 by @Giacomo.Rubino (2025-03-19T13:23:01.238Z)

Thanks for your support.
We also noticed that some other gurobi-optimization models that we use in our group do work; except this specific model (Swiss EXPANSE). Is it ok if we send this model to you, to be able to reproduce the issue?


## Post 9 by @Adrien.Albert (2025-03-21T12:43:30.546Z)

Hi @Giacomo.Rubino[@Giacomo.Rubino](https://hpc-community.unige.ch/u/giacomo.rubino)
Giacomo.Rubino:
> except this specific model (Swiss EXPANSE). Is it ok if we send this model to you, to be able to reproduce the issue?
I’ve asked multiple times, but I still don’t understand what you’re referring to. I don’t use, develop, or maintain Gurobi, so I can’t follow what you’re talking about.
Please provide a concrete example, some code or documentation about `Swiss EXPANSE` so I can understand and reproduce the issue. Without that, I can’t help.


## Post 10 by @Giacomo.Rubino (2025-03-21T14:03:30.787Z)

Ok, so Swiss EXPANSE is a cost-optimisation model of the future 2035 Swiss electricity system. This model is the one we I have issue with and is currently NOT open source. Please note that until December 2024 I could run the scenarios without issues, so something must have happened after. A reference publication can be found here (see e.g. methods for the model description):
Weather resilience of the future Swiss electricity system with very high shares of variable renewable energy sources
Collin Killenberger, Nik Zielonka, Jan-Phillipp Sasse and Evelina Trutnevyte*
Published 20 January 2025
(Weather resilience of the future Swiss electricity system with very high shares of variable renewable energy sources - IOPscience[Weather resilience of the future Swiss electricity system with very high shares of variable renewable energy sources - IOPscience](https://iopscience.iop.org/article/10.1088/2753-3751/ada77c/meta))
A documentation of a very similar model (pyPSA eur, which open source) can be found here:
https://pypsa.readthedocs.io/en/latest/[https://pypsa.readthedocs.io/en/latest/](https://pypsa.readthedocs.io/en/latest/)
And you can also find the model on git. Not that pyPSA model is referred here just because it is published on git and has relevant documentation but this is NOT the model we are having troubles with.
Now, the model I will send (Swiss EXPANSE) you is organised in 4 folders: control, EXPANSE, network, settings. In control folder, you can run the short .sh (which mentions gurobi) file that then runs a python file in the same folder. This python file then calls the EXPANSE folder. The formulation of the optimization is in EXPANSE/solve_model.py
When you run the .sh file, the model should stop and crash when you reach solve_model.py, and you will receive the error message that I sent in the first message report of this conversation.
Please let me know how do you prefer I send you the model (email or other ways) or, if you prefer, we could do it from my office if you think it’s easier.
Please not that to run this model you should also set up environment with the relevant packages (which I will also send you)
I anyway paste here the initial .sh file, that is sufficiently short to be readable here:
#!/bin/bash
#SBATCH --partition=shared-cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=110000 # in MB
#SBATCH -o run_Swiss_EXPANSE_2035_MinCost-%A_%a.out
# Activate environment
ml GCCcore/11.2.0
ml Gurobi/9.5.0
ml Python/3.9.6
. ~/baobab_python_env/bin/activate
# Send costoptimal run to cluster
echo "Running minimum cost scenario on node " $(hostname)
~/baobab_python_env/bin/python run_Swiss_EXPANSE_2035_MinCost.py


## Post 11 by @Adrien.Albert (2025-03-21T16:19:31.539Z)

Thank @Giacomo.Rubino[@Giacomo.Rubino](https://hpc-community.unige.ch/u/giacomo.rubino)
You can send your source by email  using template[email  using template](mailto:hpc@unige.ch?subject=%5BHPC%5D%5B%24CLUSTER%5D%5BSCOPE%5D%20%3CDescribe%20brively%20your%20issue%3E&body=%3CHello%20formula%20%3E%0A%0AUsername%3A%20%24USER%0ACluster%3A%20%24CLUSTER%0ASubject%3A%20%3Cslurm%7Csoftware%7Csystem%7Cconnexion%7C%20(an%20explanatory%20keyword)%20%3E%0Ajobid%3A%20%3Cseparate%20by%20space%20if%20you%20have%20more%20than%20one%3E%0A%0A%0A%3C%0A%0A%20%20%20%3CHere%20the%20mail%20body%20explaining%20ans%20describing%20%20everething%20with%20the%20MAXIMUM%20details.%20%20the%20MAXIMUM!%3E%0A%0A%20%20%20%20-%20Explain%20us%20your%20experience%20like%20we%20were%20newbies.%20%0A%20%20%20%20%20%20we%20don%27t%20spy%20your%20work%20so%20we%20don%27t%20know%20what%20is%20going%20wrong%20%3A%2F%0A%20%20%20%20-%20Explain%20the%20purpose%20of%20your%20mind%2C%20sometimes%20what%20we%20want%20is%20not%20what%20we%20think%2Fdo%2Fare%20doing.%0A%20%20%20%20-%20Give%20us%20the%20command%20you%20executed%2C%20your%20sbatch%2C%20and%20all%20relevant%20informations%0A%20%20%20%20-%20Give%20us%20the%20path%20to%20the%20relevant%20files%20in%20Baobab%20(sbatch%2C%20error%20logs%2C%20etc.)%20and%2For%20add%20them%20as%20attachments.%0A%0A%3E%0A%0A%0A%3CRespect%20formula%20(it%27s%20always%20appreciated)%3E%0A%0A%3CBye%20Bye%20formula%3E%0A%0A--%0A%24USER) with the link of this issue to keep  history :pray:


## Post 12 by @Adrien.Albert (2025-03-21T16:23:55.931Z)

- TIPS 1: I see you are using a conda environment, did you try to rebuild it after the maintenance. As we update from Rocky8 to Rocky9, it may require to update it too.


## Post 13 by @Giacomo.Rubino (2025-03-24T09:59:40.951Z)

thank you, I sent you a GitHub link by email
