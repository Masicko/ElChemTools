

abstract type abstract_simulation end

function string(SIM_list::Array{abstract_simulation})
  string_to_print = ""
  string_to_print = "$(string_to_print)["
  for SIM in SIM_list
    string_to_print = "$(string_to_print)$(string(SIM)), "
  end
  if size(SIM_list,1) > 0 
    string_to_print = string_to_print[1:end-2]
  end
  string_to_print = "$(string_to_print)]"
end

function TCtoT(TC)
  return TC+273.15
end

function TtoTC(T)
  return T-273.15
end

function pO2tosim(pO2)
  # here we can add some regularization
  treashold = 1.0e-5
  if pO2 < treashold
    treashold/100.0
  else
    return (pO2/100.0)
  end
end


function get_experiment(SIM::abstract_simulation)
  apply_checknodes(SIM, import_data_to_DataFrame(SIM), SIM.checknodes)
end

function get_sim_list(SIM_list)
  SIM_list = Array{abstract_simulation}(undef,0)
  #TODO!!
end

function general_simulation_get_fitness_factors(SIM_name::String, simulations, fitness_factors)
  if length(simulations) == length(fitness_factors)
    for (i, name) in enumerate(simulations)
      if name==SIM_name
        return fitness_factors[i]
      end
    end
    println("ERROR: SIM_name \"$(SIM_name)\" not found!")
    return throw(Exception)
  else
    println("ERROR: shape mismatch: length(simulations) == length(fitness_factors) ... $(length(simulations)) != $(length(fitness_factors))")
    return throw(Exception)
  end
end

function get_SIM_list_rectangle(TC,pO2, bias, data_set, simulations::Array{String}, fitness_factors, physical_model_name)    
    SIM_list = Array{abstract_simulation}(undef,0)
    for SIM_name in simulations
      if "CV"==SIM_name
        append!(SIM_list,[
          CV_simulation(TC, pO2, data_set=data_set, 
            fitness_factor=general_simulation_get_fitness_factors("CV", simulations, fitness_factors),
            physical_model_name=physical_model_name)...
        ]) 
        continue
      end
      if "CV(f)"==SIM_name
        append!(SIM_list,[
          CV_simulation(TC, pO2, data_set=data_set, fast_mode=true, 
            fitness_factor=general_simulation_get_fitness_factors("CV(f)", simulations, fitness_factors),
            physical_model_name=physical_model_name)...
        ])
        continue
      end
      if "EIS"==SIM_name
        append!(SIM_list,[
          EIS_simulation(TC, pO2, bias, data_set=data_set, 
            fitness_factor=general_simulation_get_fitness_factors("EIS", simulations, fitness_factors),
            physical_model_name=physical_model_name)...
        ])
        continue
      end
      if "CAP"==SIM_name
        append!(SIM_list,[
          CAP_simulation(TC, pO2)...
        ])
        continue
      end
      if "CAP-CV"==SIM_name
        append!(SIM_list,[
          CAP_simulation(TC, pO2, analytical=false)...
        ])
        continue
      end
      println("ERROR: SIM_name \"$(SIM_name)\" not found!")
      return throw(Exception)
    end
    return SIM_list
end


function filename_format_prms(; save_dir="./nouze/", prefix="", prms=Nothing, prms_names=Nothing, scripted_tuple)

  function consistency_check()
    if (size(prms_names,1) != size(scripted_tuple,1) || 
        size(prms_names,1) != size(prms,1))
      return false
    end
    return true
  end
  
  if !consistency_check()
    println("ERROR: file_format_prms(): shape mismatch (!consistency_check())")
    return throw(Exception)
  end

  scripted_dir = ""
  out_name = prefix
  for i in 1:size(scripted_tuple, 1)
    if scripted_tuple[i] == 1
      scripted_dir = string(scripted_dir,"_",prms_names[i],@sprintf("%.2f",prms[i]))
    end
    out_name = string(out_name, "_", prms_names[i], @sprintf("%.2f",prms[i]))
  end
  
  out_path = string(save_dir,scripted_dir,"/")
  out_name = string(out_name, ".csv")

  (out_path, out_name)
end

function is_between(x, a, b)
    eps = 0.000000001
    if (b - a) > 0
        if (a - eps <= x) & (x <= b + eps)
            return true
        end
    else
        if (b - eps <= x) & (x <= a + eps)
            return true
        end
    end
    return false
end
