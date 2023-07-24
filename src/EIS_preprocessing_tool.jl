using PyPlot

# TESTING DATA
# - (105 ks): 800 - 660 // 1,2,3,4,5            // 0.3 // data_set="HebbWagner_111", data_set="HebbWagner_110",  data_set="HebbWagner_100",
#
# - (?? ks): 800    // 1,2,3,4,5                // 0.3 - 0.5      // "HebbWagner_110_more_biases"  -> worse, but could be somewhat preprocessed
# - (?? ks): 700    // 1,2,3,4,5                // 0.35 - 0.55    // "HebbWagner_110_more_biases"  -> quite bad quality -> need a lot of preprocessing
#
# - (?? ks): 700    // 1,2,3,4,5                // 0.0, 0.1 - 0.55     // "HebbWagner_111_more_biases"
# - (?? ks): 600    // 1,2,3,4,5                // 0.0, 0.1 - 0.55     // "HebbWagner_111_more_biases"
#
# --> (5*3*2 = 24 ks) 700 //  1,2,3,4,5     // 0.35 : 0.05 : 0.55 //  """110, 111, more_biases"""
  
# Structure
# 
# - Saving tool
#   - together with the structure of preprocessed tree
#   - the structure should be custom to user ... using wildcard-style
#     - sth. like "/nova_data/!$(erldfkjj)
#   
# - Main Loop (for menu waiting for instructions from the user)
#   -  update the graphics (clear the window and add new graphs - Old + New )
#   -  the new is done with default settins. But the default settings are written down -> easy to be changed (commenting lines, changing values of procedures)
#   -  the changes are written in the terminal? Or do I want to construct new window also with some text field????
#       - NOT YET, but maybe later
#   - enter to apply the new filtering

# TODO ~~~~~~~
# - coding of each procedure with parameters


function processing_procedure(EIS_exp, processing_token)
  if processing_token == "cosi(3)"
    
  else
    
  end
end


function run_EIS_preprocessing_tool(;TC, pO2, bias, data_set,
                           fig_num=333, plot_option="Bode Nyq DRT RC", 
                           EIS_preprocessing_control = ElChemTools.EIS_preprocessing_control()
#                            EIS_preprocessing_control = ElChemTools.EIS_preprocessing_control(
#                                   f_interval=Nothing, 
#                                   add_inductance=0,
#                                   trim_inductance=false, 
#                                   outlayers_threshold=5.5,                                    
#                                   use_DRT=false, DRT_control=ElChemTools.DRT_control_struct()
#                            )
          )
        
  println("Processing_procedure. (Write \"abort\" for abort) ... :")
  upper_level_token = Nothing 
  abort_flag = false
  
  actual_figure = PyPlot.figure(fig_num)
  
  data_set = make_array_from_string(data_set)
  
  for   (TC_idx, TC_item) in enumerate(TC), 
        (pO2_idx, pO2_item) in enumerate(pO2),
        (bias_idx, bias_item) in enumerate(bias), 
        (data_set_idx, data_set_item) in enumerate(data_set)    
    
    SIM = EIS_simulation(TC_item, pO2_item, bias_item, data_set=data_set_item,                  
                  plot_option=plot_option, fig_num=fig_num)[1]    
    
    try
      EIS_exp = ElChemTools.import_data_to_DataFrame(SIM)
    catch 
      warning_buffer *= " ->->-> ERROR: file for (TC, pO2, bias, data_set_item) = ($(TC_item), $(pO2_item), $bias_item, $data_set_item) ... NOT FOUND! \n"                  
      continue
    end
    
    clf()
    ElChemTools.typical_plot_exp(SIM, EIS_exp)
    
    # default processing
    EIS_exp = EIS_preprocessing(EIS_exp, EIS_preprocessing_control)
    
    ElChemTools.typical_plot_exp(SIM, EIS_exp, "!default processing")
    
    command_loop_active = true
    while command_loop_active
      command_loop_active = false
      print("-> (TC, pO2, bias, data_set) = ($(TC_item), $(pO2_item), $bias_item, $data_set_item) -> command: ")
      upper_level_token = readline()
      if     upper_level_token == "abort" || upper_level_token == "a"
        abort_flag = true
        break
      elseif upper_level_token=="n"
        println("next")
        # TODO save file and continue
        # ... can be clever in the sence it will generate appropriate "repetition" number 
      elseif upper_level_token=="g"
        println("garbage")
        continue
      elseif upper_level_token=="p"
        print("processing: ")
        processing_token = readline()
        EIS_exp = processing_procedure(EIS_exp, processing_token)
      else
        println("\"$(upper_level_token)\" is not a valid command. Write \"abort\" for abort")
        command_loop_active = true
      end
    end
    if abort_flag
      break
    end
  end
  return
end











