#
#   TODO: - Some instructions for this
#         - finish the Snehurka control
#
#

# # # # # # function get_snehurka_control()
# # # # # #   #######################################################  
# # # # # #   # preparing bash output ###############################
# # # # # #   
# # # # # #   ### if true, no script is called! Just direclty run_par_study_script_wrap()
# # # # # #   direct_bool = true
# # # # # #   
# # # # # #             SIM_fitting_mode = true    #!#!#!#!#!#!#!#!#!#!
# # # # # #   
# # # # # #   bash_command = "sbatch"
# # # # # #   #bash_command = "echo"
# # # # # #   #bash_command = "julia"
# # # # # #   
# # # # # #   #mode = "test_one_prms"
# # # # # #   #mode = "only_print"
# # # # # #   mode = "go"
# # # # # #   
# # # # # #   express3_bool = true
# # # # # # 
# # # # # # 
# # # # # # 
# # # # # #   ##### TODO!!!! tohle bych mohl udelat taky obecne, at skrz jeden skriptovaci soubor muze projit vse
# # # # # #   
# # # # # #   if SIM_fitting_mode
# # # # # #     if express3_bool
# # # # # #       run_file_name = "../snehurka/run_EX3_ElChemTools_SIM_fitting.jl"
# # # # # #     else
# # # # # #       run_file_name = "../snehurka/run_ElChemTools_SIM_fitting.jl"
# # # # # #     end
# # # # # #   else
# # # # # #     if express3_bool
# # # # # #       run_file_name = "../snehurka/run_EX3_ElChemTools_par_study-prms-.jl"
# # # # # #     else
# # # # # #       run_file_name = "../snehurka/run_ElChemTools_par_study-prms-.jl"
# # # # # #     end  
# # # # # #   end
# # # # # #   
# # # # # #   return snehurka_control
# # # # # # end





function assemble_meta_SIM_fitting_single_TC(;only_return_SIM_fitting=false, mode="go", bash_command="echo", direct_bool=true)  
  #######################################################
  #######################################################
  #######################################################
  ########### SIM_fitting definition ####################


  TC = 750
  pO2 = [60]
  bias = 0.0

  data_set = "OLD_MONO_100"
#   simulations = ["EIS", "CV(f)"]
#   fitness_factors = [1.0, 5.0]
  
  simulations = ["EIS", "CV(f)"]
  fitness_factors = [1.0, 1.0]

  physical_model_name = "ysz_model_GAS_LoMA_generic"

  #####

  prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
              "A.r",      "A.DG",   "A.S",
              "R.r",      "R.DG",   "R.S",
              "O.r",      "O.DG",   "O.S",
              "nu",     "CO",     "COmm"]

#   prms_lists=      (0, 0.8, 0, 0, 0,
#                     21.9,       [-0.13, -0.05, 0.05, 0.13],    # -0.33,
#                     21.9,       [-0.13, -0.05, 0.05, 0.13],    # 0.33,
#                     21.9,       [-0.13, -0.05, 0.05, 0.13],    # 0.0,
#                     0.7,         [20.0],        [1.0])

#   prms_lists=(0, 0.0, 0, 0, 0, 22.986833011976863, 0.09356569189369754, 21.955873372516326, 0.07492013587080583, 22.67, -0.6154237534608018, 0.3414399051362852, 1, 17.26576136395294)
  prms_lists=(0, 0.0, 0, 0, 0, 22.686833011976863, 0.09356569189369754, 21.755873372516326, 0.07492013587080583, 22.47, -0.6154237534608018, 0.3414399051362852, 1, 17.26576136395294)
  prms_lists=(1, 0.17, 1, 1, 1, 21.9, 0.161, -0.416, 21.23, -0.126, 0.58, 20.977, 0.071, -0.005, 0.42, 14.052, 1.203)
  
  mask          =(0, 0, 0, 0, 0,
                  1,      1,    0,
                  1,      1,    0,
                  1,      1,    0,
                  1,    1,    1)
  lower_bounds=(0.0, -0.8,  0.0, 0.0, 0.0,
                15.5,       -0.8,       -1.0,
                15.9,       -0.8,       -1.0,
                15.7,       -0.8,       -1.0,
                0.3,    0.05,     0.05)

  upper_bounds=(1.0,  2.0,  1.0, 1.0, 1.0,
                27.5,       0.8,        1.0,
                27.9,       0.8,        1.0,
                27.7,       0.8,        1.0,
                0.95,   150,     150)


  scripted_tuple =(1, 1, 1, 1, 1,
                  1,      1,    1,
                  1,      1,    1,
                  1,      1,    1,
                  1,    1,    1)


  #######################################################

  name = "cosi_kdesi"

  #######################################################
  #############  SIM_fitting construction ###############
  #######################################################
  #######################################################

  SIM_fitting = ElChemTools.build_SIM_fitting(
                                    TC=TC,
                                    pO2=pO2,
                                    bias=bias,
                                    data_set=data_set,
                                    simulations=simulations,
                                    fitness_factors=fitness_factors,
                                    physical_model_name=physical_model_name,
                                    #
                                    prms_names=prms_names,
                                    #x0=output_prms_lists,
                                    mask=mask,
                                    lower_bounds=lower_bounds,
                                    upper_bounds=upper_bounds,
                                    #
                                    print_to_file=false,
                                    save_dir="../data/SIM_fitting/temp/"*name*"/",
                                    file_name="SIM_fitting_default.txt",
                                    #
                                    bboptimize_bool=false,
                                    iteration_count=5000,
                                    )

  pyplot = true
  plot_each_x_th = 50
  print_only_result = true

  #
  if only_return_SIM_fitting
    return SIM_fitting
  else
    return ElChemTools.meta_run_par_study(only_return_SIM_fitting=false,
                            prms_lists=prms_lists,
                            pyplot=pyplot,
                            plot_each_x_th=plot_each_x_th,
                            print_only_result=print_only_result,
                            SIM_fitting=SIM_fitting,
                            scripted_tuple=scripted_tuple,

                              ### if true, no script is called! Just direclty run_par_study_script_wrap()
                              direct_bool = direct_bool,

                            SIM_fitting_mode = true,    #!#!#!#!#!#!#!#!#!#!

                            #bash_command = "sbatch",
                            #bash_command = "echo",
                            #bash_command = "julia",
                            bash_command = bash_command,

                            #mode = "test_one_prms",
                            #mode = "only_print",
                            #mode = "go",
                            mode = mode,

                            express3_bool = true
                            )


    return SIM_fitting
  end
end





function assemble_meta_SIM_fitting_TEMPERATURE(;only_return_SIM_fitting=false, mode="go", bash_command = "sbatch")  
  #######################################################
  #######################################################
  #######################################################
  ########### SIM_fitting definition ######################


  TC = [700, 750, 850]
  pO2 = [20, 60]
  bias = 0.0

  data_set = "OLD_MONO_100"
  simulations = ["EIS", "CV(f)"]
  fitness_factors = [1.0, 5.0]

  physical_model_name = "ysz_model_GAS_LoMA_Temperature"

  #####

  prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"]

  prms_lists=(1, 0.4, 0, 0, 0,
                      0.332, 21.717,              [-0.05, 0.05], [-0.1, 0.1],
                      0.397, 21.771,              [-0.05, 0.05], [-0.1, 0.1],
                      0.337, 21.758,              [-0.05, 0.05], [-0.1, 0.1],
                      -0.05, 0.7,      9.996, 30.591,       0.713, 0.878)


### best fit BOTH ... correct zeta
#  prms_lists=      (1, 0.38, 0, 0, 0,
#                           0.2014, 21.5222,       [-0.1, -0.0618, 0.03],    [0.0558328, 0.1],
#                           0.2818, 21.2101,       [ -0.05,  0.133, 0.18],    [      0.0480154, 0.13],
#                           0.2607, 21.809,        [-0.00615, 0.1, 0.15],    [     -0.0438488, -0.1],
#                           -0.1,  0.765282,        13.1233,     37.9658,         0.15*3.86936,  0.15*5.55807)

### best fit BOTH ... wrong zeta stuff #####
#  prms_lists=(1, 0.336, 0, 0, 0,
#                      0.332, 21.817,              -0.017, 0.024,
#                      0.397, 21.171,              -0.027, 0.042,
#                      0.937, 21.858,               0.195, -0.001,
#                      -0.15, 0.65,      9.996, 30.591,       0.713, 0.878)

#  prms_lists=      (1, 0.32, 0, 0, 0,
#                           0.2414, 21.9222,       [-0.1, -0.0618,  0.1],    [0.0,  0.0558328, 0.1],
#                           0.2818, 21.2101,       [ 0.0,  0.133,   0.2],    [      0.0480154, 0.13],
#                           0.2607, 21.809,        [-0.1, -0.00615, 0.1],    [     -0.0438488, -0.1],
#                           -0.133,  0.765282,     13.1233,     37.9658,         0.15*3.86936,  0.15*5.55807)

#  prms_lists=(1, 0.27, 0.0, 0.0, 0.0,
#               0.22236513700278862, 21.92,               0.00465798980401979, 0.055,
#               0.33958803427894996, 21.21,              -0.06332811104733665, 0.048,
#               0.33062408382716696, 21.809,              0.15698240545620798, -0.04,
#               -0.253637378507218, 0.76,        9.999491268858524, 37.97,         0.21670880659050948, 0.834)

  mask          =(0, 0, 0, 0, 0,
                  1, 1,           1, 1,
                  1, 1,           1, 1,
                  1, 1,           1, 1,
                  1, 1,     1, 1,     1, 1)

  lower_bounds=(0.0, -0.4, 0.0, 0.0, 0.0,
                -3, 20.5,         -0.5, -0.8,
                -3, 20.5,         -0.5, -0.8,
                -3, 20.5,         -0.5, -0.8,
                -0.15, 0.65,     0, 0.05,    0, 0.05)


  upper_bounds=(1.0, 1.0, 1.0, 1.0, 1.0,
                3, 26.5,         0.5, 0.8,
                3, 26.5,         0.5, 0.8,
                3, 26.5,         0.5, 0.8,
                0.00, 0.95,     30, 100.0,    30, 100.0)


  scripted_tuple =(1, 1, 1, 1, 1,
                  1, 1,           1, 1,
                  1, 1,           1, 1,
                  1, 1,           1, 1,
                  1, 1,       1, 1,     1, 1)


  #######################################################
  name = "cosi_kdesi"

  #######################################################
  #############  SIM_fitting construction ###############
  #######################################################
  #######################################################

  SIM_fitting = ElChemTools.build_SIM_fitting(
                                    TC=TC,
                                    pO2=pO2,
                                    bias=bias,
                                    data_set=data_set,
                                    simulations=simulations,
                                    fitness_factors=fitness_factors,
                                    physical_model_name=physical_model_name,
                                    #
                                    prms_names=prms_names,
                                    #x0=output_prms_lists,
                                    mask=mask,
                                    lower_bounds=lower_bounds,
                                    upper_bounds=upper_bounds,
                                    #
                                    print_to_file=false,
                                    save_dir="../data/SIM_fitting/temp/"*name*"/",
                                    file_name="SIM_fitting_default.txt",
                                    #
                                    bboptimize_bool=false,
                                    iteration_count=15000,
                                    )

  pyplot = false
  plot_each_x_th = 50
  print_only_result = true

  #
  if only_return_SIM_fitting
    return SIM_fitting
  else
    return ElChemTools.meta_run_par_study(only_return_SIM_fitting=false,
                            prms_lists=prms_lists,
                            pyplot=pyplot,
                            plot_each_x_th=plot_each_x_th,
                            print_only_result=print_only_result,
                            SIM_fitting=SIM_fitting,
                            scripted_tuple=scripted_tuple,

                              ### if true, no script is called! Just direclty run_par_study_script_wrap()
                              direct_bool = false,

                            SIM_fitting_mode = true,    #!#!#!#!#!#!#!#!#!#!

                            bash_command = bash_command,
                            #bash_command = "echo",
                            #bash_command = "julia",

                            #mode = "test_one_prms",
                            #mode = "only_print",
                            mode = mode,

                            express3_bool = true
                            )


    return SIM_fitting
  end
end


function assemble_meta_SIM_fitting_TEMPERATURE_with_S(;only_return_SIM_fitting=false, mode="go", bash_command = "sbatch")  
  #######################################################
  #######################################################
  #######################################################
  ########### SIM_fitting definition ######################


  TC = [700, 750, 850]
  pO2 = [20, 60]
  bias = 0.0

  data_set = "OLD_MONO_100"
  simulations = ["EIS", "CV(f)"]
  fitness_factors = [1.0, 5.0]

  physical_model_name = "ysz_model_GAS_LoMA_Temperature"

  #####

  prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",     "A.S",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",     "R.S",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",     "O.S",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"]

  prms_lists=(1, 0.4, 0, 0, 0,
                      0.332, 21.717,              [-0.05, 0.05], [-0.1, 0.1],     0,
                      0.397, 21.771,              [-0.05, 0.05], [-0.1, 0.1],     0, 
                      0.337, 21.758,              [-0.05, 0.05], [-0.1, 0.1],     0,
                      -0.05, 0.7,      9.996, 30.591,       0.713, 0.878)


### best fit BOTH ... correct zeta
#  prms_lists=      (1, 0.38, 0, 0, 0,
#                           0.2014, 21.5222,       [-0.1, -0.0618, 0.03],    [0.0558328, 0.1],
#                           0.2818, 21.2101,       [ -0.05,  0.133, 0.18],    [      0.0480154, 0.13],
#                           0.2607, 21.809,        [-0.00615, 0.1, 0.15],    [     -0.0438488, -0.1],
#                           -0.1,  0.765282,        13.1233,     37.9658,         0.15*3.86936,  0.15*5.55807)

### best fit BOTH ... wrong zeta stuff #####
#  prms_lists=(1, 0.336, 0, 0, 0,
#                      0.332, 21.817,              -0.017, 0.024,
#                      0.397, 21.171,              -0.027, 0.042,
#                      0.937, 21.858,               0.195, -0.001,
#                      -0.15, 0.65,      9.996, 30.591,       0.713, 0.878)

#  prms_lists=      (1, 0.32, 0, 0, 0,
#                           0.2414, 21.9222,       [-0.1, -0.0618,  0.1],    [0.0,  0.0558328, 0.1],
#                           0.2818, 21.2101,       [ 0.0,  0.133,   0.2],    [      0.0480154, 0.13],
#                           0.2607, 21.809,        [-0.1, -0.00615, 0.1],    [     -0.0438488, -0.1],
#                           -0.133,  0.765282,     13.1233,     37.9658,         0.15*3.86936,  0.15*5.55807)

#  prms_lists=(1, 0.27, 0.0, 0.0, 0.0,
#               0.22236513700278862, 21.92,               0.00465798980401979, 0.055,
#               0.33958803427894996, 21.21,              -0.06332811104733665, 0.048,
#               0.33062408382716696, 21.809,              0.15698240545620798, -0.04,
#               -0.253637378507218, 0.76,        9.999491268858524, 37.97,         0.21670880659050948, 0.834)

  mask          =(0, 0, 0, 0, 0,
                  1, 1,           1, 1,     1,
                  1, 1,           1, 1,     1,
                  1, 1,           1, 1,     1,
                  1, 1,     1, 1,     1, 1)

  lower_bounds=(0.0, -0.4, 0.0, 0.0, 0.0,
                -3, 20.5,         -0.5, -0.8,   -1.0,
                -3, 20.5,         -0.5, -0.8,   -1.0,
                -3, 20.5,         -0.5, -0.8,   -1.0,
                -0.15, 0.65,     0, 0.05,    0, 0.05)


  upper_bounds=(1.0, 1.0, 1.0, 1.0, 1.0,
                3, 26.5,         0.5, 0.8,       1.0,
                3, 26.5,         0.5, 0.8,       1.0,
                3, 26.5,         0.5, 0.8,       1.0,
                0.00, 0.95,     30, 100.0,    30, 100.0)


  scripted_tuple =(1, 1, 1, 1, 1,
                  1, 1,           1, 1,       1,
                  1, 1,           1, 1,       1,
                  1, 1,           1, 1,       1,
                  1, 1,       1, 1,     1, 1)


  #######################################################
  name = "cosi_kdesi"

  #######################################################
  #############  SIM_fitting construction ###############
  #######################################################
  #######################################################

  SIM_fitting = ElChemTools.build_SIM_fitting(
                                    TC=TC,
                                    pO2=pO2,
                                    bias=bias,
                                    data_set=data_set,
                                    simulations=simulations,
                                    fitness_factors=fitness_factors,
                                    physical_model_name=physical_model_name,
                                    #
                                    prms_names=prms_names,
                                    #x0=output_prms_lists,
                                    mask=mask,
                                    lower_bounds=lower_bounds,
                                    upper_bounds=upper_bounds,
                                    #
                                    print_to_file=false,
                                    save_dir="../data/SIM_fitting/temp/"*name*"/",
                                    file_name="SIM_fitting_default.txt",
                                    #
                                    bboptimize_bool=false,
                                    iteration_count=15000,
                                    )

  pyplot = false
  plot_each_x_th = 50
  print_only_result = true

  #
  if only_return_SIM_fitting
    return SIM_fitting
  else
    return ElChemTools.meta_run_par_study(only_return_SIM_fitting=false,
                            prms_lists=prms_lists,
                            pyplot=pyplot,
                            plot_each_x_th=plot_each_x_th,
                            print_only_result=print_only_result,
                            SIM_fitting=SIM_fitting,
                            scripted_tuple=scripted_tuple,

                              ### if true, no script is called! Just direclty run_par_study_script_wrap()
                              direct_bool = false,

                            SIM_fitting_mode = true,    #!#!#!#!#!#!#!#!#!#!

                            bash_command = bash_command,
                            #bash_command = "echo",
                            #bash_command = "julia",

                            #mode = "test_one_prms",
                            #mode = "only_print",
                            mode = mode,

                            express3_bool = true
                            )


    return SIM_fitting
  end
end






function assemble_meta_SIM_fitting_TEMPERATURE_free_nu(;only_return_SIM_fitting=false, mode="go", bash_command = "sbatch")  
  #######################################################
  #######################################################
  #######################################################
  ########### SIM_fitting definition ######################


  TC = [700, 750, 800, 850]
  pO2 = [20, 60]
  bias = 0.0

  data_set = "OLD_MONO_100"
  simulations = ["EIS", "CV(f)"]
  fitness_factors = [1.0, 3.]

  physical_model_name = "ysz_model_GAS_LoMA_Temperature"

  #####

  prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_700", "nu_750", "nu_800", "nu_850",     "CO_B", "CO_C",     "COmm_B", "COmm_C"]

  prms_lists=      (1, 0.38, 0, 0, 0,
                           0.2014, 21.5222,       [-0.1, -0.0618, 0.03],    [0.0558328, 0.1],
                           0.2818, 21.2101,       [ -0.05,  0.133, 0.18],    [      0.0480154, 0.13],
                           0.2607, 21.809,        [-0.00615, 0.1, 0.15],    [     -0.0438488, -0.1],
                           0.8,  0.7, 0.6, 0.5,        13.1233,     37.9658,         0.15*3.86936,  0.15*5.55807)

  mask          =(0, 0, 0, 0, 0,
                  1, 1,           1, 1,
                  1, 1,           1, 1,
                  1, 1,           1, 1,
                  1, 1, 1, 1,     1, 1,     1, 1)

  lower_bounds=(0.0, 0.0, 0.0, 0.0, 0.0,
                -3, 20.5,         -0.5, -0.8,
                -3, 20.5,         -0.5, -0.8,
                -3, 20.5,         -0.5, -0.8,
                0.3, 0.3, 0.3, 0.3,     0, 0.05,    0, 0.05)


  upper_bounds=(1.0, 1.0, 1.0, 1.0, 1.0,
                3, 26.5,         0.5, 0.8,
                3, 26.5,         0.5, 0.8,
                3, 26.5,         0.5, 0.8,
                0.95, 0.95, 0.95, 0.95,     30, 100.0,    30, 100.0)


  scripted_tuple =(1, 1, 1, 1, 1,
                  1, 1,           1, 1,
                  1, 1,           1, 1,
                  1, 1,           1, 1,
                  1, 1, 1, 1,      1, 1,     1, 1)


  #######################################################

  name = "cosi_kdesi"

  #######################################################
  #############  SIM_fitting construction ###############
  #######################################################
  #######################################################

  SIM_fitting = ElChemTools.build_SIM_fitting(
                                    TC=TC,
                                    pO2=pO2,
                                    bias=bias,
                                    data_set=data_set,
                                    simulations=simulations,
                                    fitness_factors=fitness_factors,
                                    physical_model_name=physical_model_name,
                                    #
                                    prms_names=prms_names,
                                    #x0=output_prms_lists,
                                    mask=mask,
                                    lower_bounds=lower_bounds,
                                    upper_bounds=upper_bounds,
                                    #
                                    print_to_file=false,
                                    save_dir="../data/SIM_fitting/temp/"*name*"/",
                                    file_name="SIM_fitting_default.txt",
                                    #
                                    bboptimize_bool=false,
                                    iteration_count=9000,
                                    )

  pyplot = false
  plot_each_x_th = 50
  print_only_result = true

  #
  if only_return_SIM_fitting
    return SIM_fitting
  else
    return ElChemTools.meta_run_par_study(only_return_SIM_fitting=false,
                            prms_lists=prms_lists,
                            pyplot=pyplot,
                            plot_each_x_th=plot_each_x_th,
                            print_only_result=print_only_result,
                            SIM_fitting=SIM_fitting,
                            scripted_tuple=scripted_tuple,

                              ### if true, no script is called! Just direclty run_par_study_script_wrap()
                              direct_bool = false,

                            SIM_fitting_mode = true,    #!#!#!#!#!#!#!#!#!#!

                            bash_command = bash_command,
                            #bash_command = "echo",
                            #bash_command = "julia",

                            #mode = "test_one_prms",
                            #mode = "only_print",
                            #mode = "go",
                            mode = mode,

                            express3_bool = true
                            )


    return SIM_fitting
  end
end  

