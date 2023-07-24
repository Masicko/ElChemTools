#!/usr/local/pkg/Installs/linux-ubuntu18.04-skylake_avx512/gcc7.4.0/julia/1.3.1/6cvgicvsj2rrblel/bin/julia
#
#SBATCH --job-name=NYQ
#SBATCH -p express3
#SBATCH --time=0-6

include(string(pwd(),"/../src/ElChemTools.jl"))

@time ElChemTools.run_par_study_script_wrap(ARGS[1], ARGS[2], ARGS[3], ARGS[4], ARGS[5], ARGS[6], ARGS[7], ARGS[8], ARGS[9], ARGS[10], ARGS[11])
