# Gurobi license expired

- Source: https://hpc-community.unige.ch/t/3323

- Created: 2024-02-19T14:17:19.751Z

- Posts: 8

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Zongfei.Wang (2024-02-19T14:17:19.810Z)

## Primary informations
Username: wangzo
Cluster: Yggdrasil and Baobab
## Description
Hi. It seems like Gurobi license has expired both on Yggdrasil and Baobab. Could you please check? Thanks.
(Sorry to cross posting. I previously sent the issue in a wrong category.)
## Steps to Reproduce
To run any optimization problem with Gurobi to see if Gurobi license is still valid.
## Expected Result
Normally, optimization problems can be solved with Gurobi solver.
## Actual Result
E.g., for the job ID 31131275 on Yggdrasil and other jobs I submitted since 16.02.2024, I got the following messages:
> gurobipy.GurobiError: Request denied: license expired
> ERROR: Solver (gurobi) returned non-zero return code (1)
> ERROR: See the solver log above for diagnostic information.


## Post 2 by @Adrien.Albert (2024-02-19T16:13:22.656Z)

Hi @Zongfei.Wang[@Zongfei.Wang](https://hpc-community.unige.ch/u/zongfei.wang)
The license seems to have expired on 15/02/2024.
However, I’ve corrected something on Yggdrasil and communication with the server license (on Baobab) is working again.
I’ve contacted the Gurobi team to extend the license. Keep up to date


## Post 3 by @Adrien.Albert (2024-02-19T16:47:05.451Z)

@Zongfei.Wang[@Zongfei.Wang](https://hpc-community.unige.ch/u/zongfei.wang)
The license have been updated, could you confirm everything is working well?


## Post 4 by @Zongfei.Wang (2024-02-20T08:38:05.116Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
I made test runs again on Yggdrasil and Baobab but they still don’t work, with messages as above saying ‘license expired’.
I checked the license status on Yggdrasil and Baobab and it says license works fine, with messages like below:
```
(yggdrasil)-[wangzo@login1 control]$ module load Gurobi
(yggdrasil)-[wangzo@login1 control]$ gurobi_cl --tokens
Checking status of Gurobi token server 'admin1.baobab'...
Token server functioning normally.
Maximum allowed uses: 4096, current: 0
```
Could you please further check? Thanks! - Zongfei


## Post 5 by @Adrien.Albert (2024-02-20T08:51:03.081Z)

Hi @Zongfei.Wang[@Zongfei.Wang](https://hpc-community.unige.ch/u/zongfei.wang)
Could you give a quick test to get this error message ?


## Post 6 by @Adrien.Albert (2024-02-20T09:05:07.581Z)

Hi @Zongfei.Wang[@Zongfei.Wang](https://hpc-community.unige.ch/u/zongfei.wang)
It’s good to know that when Gurobi staff say extend license, it means generate another one.
I’ve updated the license file, it seems to be working fine.


## Post 7 by @Zongfei.Wang (2024-02-20T09:08:04.599Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
Sure, I ran job 31224141 on Yggdrasil this morning. And the error message is like:
```
Set parameter GURO_PAR_SPECIAL
Set parameter TokenServer to value "admin1.baobab.hpc.unige.ch"
Traceback (most recent call last):
  File "<stdin>", line 5, in <module>
  File "/home/users/w/wangzo/baobab_python_env/lib/python3.9/site-packages/pyomo/solvers/plugins/solvers/GUROBI_RUN.py", line 62, in gurobi_run
    model = read(model_file)
  File "src/gurobipy/gurobi.pxi", line 3571, in gurobipy.read
  File "src/gurobipy/gurobi.pxi", line 80, in gurobipy.gurobi.read
  File "src/gurobipy/gurobi.pxi", line 32, in gurobipy.gurobi._getdefaultenv
  File "src/gurobipy/env.pxi", line 62, in gurobipy.Env.__init__
gurobipy.GurobiError: Request denied: license expired
ERROR: Solver (gurobi) returned non-zero return code (1)
ERROR: See the solver log above for diagnostic information.
Traceback (most recent call last):
  File "/home/users/w/wangzo/V2_Z10_0_N10_Rx_V2_r0/control/run_MinCost.py", line 84, in <module>
    model,results = solve_model_mincost(model)
  File "/home/users/w/wangzo/V2_Z10_0_N10_Rx_V2_r0/EXPANSE/solve_model.py", line 52, in solve_model_mincost
    results = SolverFactory(opts['solver_name']).solve(model,
  File "/home/users/w/wangzo/baobab_python_env/lib/python3.9/site-packages/pyomo/opt/base/solvers.py", line 596, in solve
    raise ApplicationError(
pyomo.common.errors.ApplicationError: Solver (gurobi) did not exit normally
```


## Post 8 by @Zongfei.Wang (2024-02-20T09:22:50.907Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
I just made test runs again both on Yggdrasil and Babab. And Gurobi is now working on both clusters. Thank you so much!
