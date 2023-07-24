################## Typical basic usage ######################
#
# Function "meta_run_par_study()" in "ElChemTools.jl" file can perform parametrical study.
# Its behavior is still defined in the body of the function. 
# It can 
#   1) call the "run_par_study()" function directly (via parameter direct_bool) 
#   OR 
#   2) by using a custom made external script file (for example used for running a job on a computational cluster).
# 
# You specify the parametrical study by
#   name = "GAS_LoMA_ULTRA_level_2"     # of the parametrical study
#   prms_names = ["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"]
#   prms_lists 
#   scripted_tuple        # which defines, how many cluster jobs (external calls of the script) will be called
#                         # and how many sets of prms will be tried in the called "run_par_study()"
#                         # 
#                         # (1, 1, 0, 0, 0, 0, 0) means al combinations of A0 and R0 will be distributed to 
#                         # differents jobs (calls of the script). And all combinations of the rest of parameters 
#                         # will be tried in each job (e.g. cluster node) by a funtion "run_par_study()"
#                         # ... simply ... if you just want run par_study on one computer (no cluster), use (0,0,0,0,0,0,0) :)
#                         
#   TC, pO2, bias                   # experimental setting lists can be used instead of single value to perform more simulations
#   simulations = ["EIS", "CV"]             # which type of simulation should be performed for each experimental setting
#   save_dir = "../data/par_studies/"       # by default. The par_study will be saved in a folder with its name
#
# All the info about par_study will be stored in "__metafile_par_study.txt"
# The behavior of meta_run_par_study is further specified by variables 
#   bash_command    # which command use to run the outer script
#   mode = "test_one_prms"  --> from the whole prms_lists chooses "the middle" one set or prms and perform par_study
#                               serves as a test of the path before the run of big dataset par_study
#   mode = "only_print"     --> only print the atributes of actually defined par_study 
#                               as well as the count of future (jobs / number of prms per each job)
#   mode = "go"             --> runs the par_study
#   run_file_name   # which script file will be run by the bash_command
# 
#
##############################################################
using ElChemTools

# perform the prepared parametrical study and save it to defined folder (the default is )
ElChemTools.meta_run_par_study()

# you can then import info about par_study from its metafile
# ... and store it in an instance of par_study_struct
# ... but still without any imported simulated data
my_par_study = ElChemTools.par_study_import_info_from_metafile("GAS_LoMA_ULTRA_level_2")

# you can import the whole simulated data to the REPL memory (to "my_par_study")
# but it is not necessary! (and often NOT feasible due to its size)
#ElChemTools.par_study_import_data!(my_par_study)

# you can use the data, compute fitness errors and store (in "my_par_study") only errors (far less memmory demands)
ElChemTools.par_study_err_evaluation!(my_par_study)

# the results of parametrical study can be viewed by error projection to each parameter
# parameter "count" simply says, how many "the best" errors will be shown
ElChemTools.par_study_plot_err_projections(my_par_study, count=8) 

# and you can view the resulting simulation curves (EIS or CV)
# parameters from=4, till=8 says it shows the best {4, 5, 6, 7, 8} experiments
# ... the default value is from=1, till=1, which means "print just the best"
# the error is computed as the sum of fitness function of all simulations from SIM_list
ElChemTools.par_study_plot_the_best(my_par_study, from=1, till=1)

# if you want to consider just some special subset of SIM_list, you can apply filters
#ElChemTools.par_study_plot_the_best(my_par_study, pO2=80)
#ElChemTools.par_study_plot_the_best(my_par_study, simulations=["CV"], pO2=[80])
#ElChemTools.par_study_plot_the_best(my_par_study, simulations=["CV"])




###############################
#
#   ADVANCED
#
###############################
# You can view the SIM_list (which simulations were performed ... EIS or CV for which experimental condition)
#string(my_par_study.SIM_list)

# and then, by tuple, directly choose, which simulations will be viewed
#ElChemTools.par_study_plot_the_best(my_par_study, )
