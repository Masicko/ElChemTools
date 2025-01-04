################## Typical basic usage ######################
#
# three modes:  "EEC" - equivalent electric circuit of RC & RCPE in series
#                       R1, C1        = paralel RC circuit
#                       R2, C2, alpha = paralel R + CPE circuit
#               "exp" - experimetnal data
#               "sim" - simulation
#  Algoritms:
#  - alg = "Tikhonov":
#      lambda - regularization parameter ... minimization of F = || A*x - b || + || lambda*I*x ||
#  - alg = "Gold" ... deconvolution https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9400974/pdf/CPHC-23-0.pdf
#      k_Gold - iteration number (in the role of regularization parameter)
# 
#  One can use EIS_preprocessing (viz EIS_preprocessing_usage.jl)
#
#
#      - There is also a possibility to subtract indunctance using the keyword
#             subtract_inductance_fac=1.0
#             where the factor 1.0 means it subtracts the inductance which comes out from gDRT analysis. 
#             0.9 takes a little bit lower inductance.
#      - export_in_f ... means the DRT spectrum h(tau) becames h(f)
#
#  DRT Peak analysis:      
#      - divide_R_peaks = "<DIV> <METHOD>" ... the DRT spectrum is divaded by <DIV> algoritm to segments 
#                    each containing (hopefully :) ) one peak and then the total resistance is divided 
#                    into "peaks" by specified <METHOD>
#             <DIV>  = "valley" -> separators are placed in the valley (local minima of h(t))
#                    = "curv"   -> separators are found by curvature (2nd derivative) of h(t)
#             <METHOD> = "sum"  -> h(t) is summed withing each "peak segment" nad frequency f is 
#                                  at the highest point of a segment h(t)
#                      = "gauss"-> gaussian peaks are fitted to h(t) and then the h(t) is 
#                                  redistributed for each "t" to appropriate "peaks" according to the ratio
#                                  of present peaks in the point "t",
#                                  frequency f is the mean of the gaussian peak
#            and than the f and R of each peak is evaluated and exported to a file with prefix "Rx..."
#      - R_peaks_plot_bool = true ... the division and gaussian peak fitting is shown
##############################################################
using ElChemTools

#
ElChemTools.test_DRT(mode="EEC", R1=0, C1=0.001, R2=1.0, C2=0.01, alpha=1.0, lambda=0.0)

#
ElChemTools.test_DRT(mode="EEC", R1=0, C1=0.001, R2=1.0, C2=0.01, alpha=0.9, lambda=0.01)

#
ElChemTools.test_DRT(mode="exp", TC=800, pO2=40, bias=0.0, lambda=0.001)

# new data_set called "POLY" for polycrystaline. The monocrystaline is "MONO_110"
ElChemTools.test_DRT(mode="exp", TC=800, pO2=20, bias=collect(0.0 : 0.05 : 0.2), lambda=0.0, 
       backward_check=false, plot_option="DRT Nyq Bode Rtau", data_set="examples/default_EIS_example.z", 
       tau_min_fac=1.1, tau_max_fac=1.1, tau_range_fac=20.0, f_range=(0.6, 9000, 1.2))
       
# #
# # CAP_comparison feature ---------------------------> CAP_comparison switch
# # This means
# #  - the parameter "rR" is set to "1"
# #  - and CAP_simulation(analytical=false) is called
# # 
# ElChemTools.test_DRT(mode="sim",       TC=850, pO2=100, bias=collect(-0.5 : 0.02 : 0.5),  
                                   
#                                    prms_names=["kappaA", "kappaR", "kappaO", 
#                                                "rA", "rR", "rO",         "rB", "rC",     
#                                                "DGA", "DGR", "DGO",     
#                                                "DD", "nu", "separate_vacancy",       "sites_Om0", "sites_Om1"  ], 
#                                    prms_values=(0.0, 0.0, 0.0,
#                                                 27.5, 17.5, 27.5,        1, 1,#collect(21 : -0.5 : 21),       
#                                                 0.4, -0.5, collect(-0.0 : 0.01 : 0.0),     
#                                                 [90]*1.0e-14, collect(0.85 : 0.05 : 0.85), true,       1/4, 1/2    ),
#                                    CAP_comparison=true, plot_bool=false);
#  
# export file option:
#       - export_file = "" is defaul -> means export Nothing
#                     = "foo.txt" -> will export data from DRT with information as TC, pO2, bias, data_set, lambda atc.
#       - export_append = "true" is default -> if file with the specified name exists, it add new rows to it
#                       = "false" -> overrides file if exists
