#!/bin/bash
DISORDER_STRENGTHS=(0.5 1 1.5 2 2.5 3 3.5 4 4.5 5)
for r in {0..9..1}
	do
		sbatch --job-name=otoc_simulation_h$r --output=otoc_simulation_h${r}-%j.out slurm/Run_N16.jl ${DISORDER_STRENGTHS[$r]}
		sleep 1
done
