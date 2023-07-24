################## Typical basic usage ######################
#
# use this funtion to call simple simulation of the experiments
#   ElChemTools.simple_run()
#
# parameters:
#   TC, pO2, (bias=0.0) >>> the experimental setting
#         note: experimental data are for
#                 TC \in    {700, 750, 800, 850}
#                 pO2 \in   {0, 20, 40, 60, 80, 100}
#                 bias \in  {-1.0, -0.5, 0.0, 0.5, 1.0} but for 700 is missing {-1.0, 1.0}
#         note: You can write an array instead of 1 value for this parameters (as for prms_values described below)
#         WARNING: pO2 < 1.0e-5 => pO2 = 1.0e-5 ... regularization treashold
#   simulations     >>> =["EIS"] or =["CV"] or =["EIS", "CV"] for perform both methods for all experimental settings
#                       or now also ["CAP"] means differential capacitance computed with R0 = 0
#                       =["CV(f)"]   
#                               \-------TODO !!!
#                       
#   pyplot          >>> 0 = plot nothing, 1 = finall plot of the experiment, 2 = plot details
#   prms_names      >>> in array or tuple of strings specify, which parametr of the model should be changed
#   prms_values     >>> spedify, to which value it will be changed
#         note: You can change any of the parameters in YSZParameters class by typing its name, e.g. x_frac = 0.50
#         note: You can do small parametrical study by typing an array of values instead of one value
#              For example:  prms_names = ["DGA", "DGR", "DGO"]
#                            prms_values = [0.4, 0.5, 0.3]                       # 1 set of parameters
#                            prms_values = [0.4, [0.4, 0.5, 0.6], 0.3]           # 3 sets of parameters
#                            prms_values = [0.4, collect(0.4 : 0.1 : 0.6), 0.3]  # 3 sets of parameters
#   show_experiment >>> if the experimental setting is in the database, it shows fitting error and (if pyplot > 1) plot experimental data
#   fig_size        >>> standard PyPlot tuple (x,y)
#
# side-notes:
#   - this is using LoMA model
#   - expA \in {0, 1} determines, if the kinetics of the reaction will be switched to pure exponential form (from LoMA)
#   - "L" inductance
#   - TC = (700, 750, 800, 850)°C  => use => DD = ( 1.277, 2.92, 5.35, 9.05)e-13
#       - I know it should influence also the other parameters like \nu, but this is the first approximation to match overall resistance of the system
#
#
# #### update: CAP stuff ####
# 
# ["CAP"] simulation is also implemented
#   - overpotential = 0 is the equilibrium voltage !!!
#   - parametr R0 is set to 0 (but it can be overwriten to non-zero by prms_values)
#   - there are no experimental data, you must use "use_experiment=false" switch! 
#   >>> ["CAP"] Analytical version  - using direct capacitance analytical solution for defined voltage checknodes in simulation_struct
#
#   >>> ["CAP-CV"]     - using voltammetry simulation with default voltrate = 0.000001
#                      - calculate the capacitance from a rescaled current ... C = I / (voltrate * direction_of_sweep)
#                           
#                         
#                       
##############################################################

using ElChemTools

ElChemTools.simple_run(EIS_simulation(850, [60], 0.0, physical_model_name="ysz_model_GAS_LoMA_Temperature", plot_legend=false), pyplot=1,  
                      prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                          "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                          "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                          "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                          "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
                        
                      prms_values=(1, 0.0, 0.0, 0.0, 0.0, 
                          0.0, 24.4626,            0.0, 0.33,
                          0.0, 22.63,              0.0, 0.0088,
                          0.0, 21.94,              0.0, 0.17,
                          0.0, 0.85,       0.0, 3.95,      0.0, 6.818*0.15)
                      ,use_experiment=false);
