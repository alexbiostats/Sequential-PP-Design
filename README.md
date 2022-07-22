# Sequential-PP-Design
Accompanying R code to replicate simulation results for the "Bayesian and Frequentist approaches to sequential monitoring for futility in oncology
basket trials: A comparison of Simon’s two-stage design and Bayesian predictive probability monitoring with information sharing across baskets" manuscript to appear in *PLOS ONE*. Briefly, Bayesian predictive probabilities (PP) can be used for sequential monitoring in clinical trials. Specifically, in this manuscript, simulation studies of basket trials with 10 one-arm baskets with binary outcomes are implemented. Comparisons are made between Simon's two-stage minimax design, a Bayesian approach with PP for continual futility monitoring, and a Bayesian approach that further considers information sharing across baskets with multi-source exchangeability models (MEMs).

**function_mem_seq_bt.R**: a set of functions for implementing the simulation study

**snowfall_sim_memseqbt.R**: code to run simulations with parallelization

**seqPP_result_code.R**: code to take simulation results and summarize with tables and figures

**mem_seq_bt_results_part1.txt**: text file of simulation results (part 1 of 2, need to merge together with part 2 for summarizing results using `seqPP_result_code.R`)

**mem_seq_bt_results_part2.txt**: text file of simulation results (part 2 of 2, need to merge together with part 1 for summarizing results using `seqPP_result_code.R`)
