using CSV
using DataFrames
using PyPlot



const  CAP_standard_figure_num = 7



######################################


mutable struct CAP_simulation <: abstract_simulation
  TC::Float32
  pO2::Float32  
  data_set::String
  #
  physical_model_name::String
  #
  upp_bound::Float32
  low_bound::Float32
  dx_exp::Float64
  sample::Int64
  voltrate::Float32
  fig_size::Tuple
  fig_num::Int32
  #
  analytical::Bool
  checknodes::Any
  fitness_factor::Float64
  #
  name::String
  ID::Int16
  
  CAP_simulation(n::Nothing) = new()
end

function string(SIM::CAP_simulation)
  return "CAP_sim_TC_$(SIM.TC)_pO2_$(SIM.pO2)"
end

function CAP_simulation(TC, pO2, bias=0.0; data_set="MONO_110", analytical=true, physical_model_name="ysz_model_GAS_LoMA_Temperature", upp_bound=1.0, low_bound=-1.0, dx_exp=-9, sample=40, voltrate=0.00001, fig_size=(9, 6), fig_num=CAP_standard_figure_num, fitness_factor=10.0, checknodes=get_shared_checknodes(CAP_simulation(nothing)))
  output = Array{abstract_simulation}(undef,0)
  for TC_item in TC
    for pO2_item in pO2
      this = CAP_simulation(nothing)
      
      this.TC = TC_item
      this.pO2 = pO2_item  
      this.data_set = data_set
      #
      this.physical_model_name = physical_model_name
      #
      this.upp_bound = upp_bound
      this.low_bound = low_bound
      this.dx_exp = dx_exp
      this.sample = sample
      this.voltrate = voltrate
      this.fig_size = fig_size
      this.fig_num = fig_num
      #
      this.analytical = analytical
      this.checknodes = checknodes
      this.fitness_factor = fitness_factor
      #
      this.name = analytical ? "CAP(a)" : "CAP(CV)"
      this.ID = 3
      
      push!(output, this)
    end
  end
  return output
end 


#######################################
#### run simulation stuff ####
function typical_run_simulation(SIM::CAP_simulation, prms_names_in, prms_values_in, pyplot::Int=0) 
  if SIM.analytical
        
    ysz_experiments.run_new(
        out_df_bool=true, voltammetry=false, pyplot=(pyplot == 2 ? true : false),
        dlcap=true, dlcap_analytical=true, checknodes=SIM.checknodes[:,1], physical_model_name=SIM.physical_model_name,
        dx_exp=SIM.dx_exp, sample=SIM.sample, upp_bound_eta=SIM.upp_bound, low_bound_eta=SIM.low_bound, 
        T=TCtoT(SIM.TC), pO2=pO2tosim(SIM.pO2), data_set=SIM.data_set,
        prms_names_in=prms_names_in,
        prms_values_in=prms_values_in,
    )
  else
    ysz_experiments.run_new(
        out_df_bool=true, voltammetry=true, pyplot=(pyplot == 2 ? true : false),
        dlcap=true,
        dx_exp=SIM.dx_exp, sample=SIM.sample, upp_bound=SIM.upp_bound, low_bound=SIM.low_bound, voltrate=SIM.voltrate,
        T=TCtoT(SIM.TC), pO2=pO2tosim(SIM.pO2),
        prms_names_in=prms_names_in,
        prms_values_in=prms_values_in,
    )
  end
end


#### checknodes stuff ####
function CAP_get_checknodes(start_n,upper_n,lower_n,end_n,step_n)
    # nodes of voltage U with direction of sweep where I will be fitted
    # [U direction; ... ]
    # the list must be sorted !!! w.r.t. U and direction 1 -> -1 -> 1
    # like in the real CAP measurement
    ll = start_n : step_n : upper_n
    res = [ll [1 for i in 1:size(ll,1)]]
    ll = upper_n : -step_n : lower_n
    res = vcat(res, [ll [-1 for i in 1:size(ll,1)]])
    ll = lower_n : step_n : end_n
    res = vcat(res, [ll [1 for i in 1:size(ll,1)]])
    return res
end

function get_shared_checknodes(SIM::CAP_simulation)
  return CAP_get_checknodes(0.0,0.99,-0.99, 0.0,0.01)
end

function get_nice_checknodes(sim::CAP_simulation)
  # experimental range is 0.00104 | 0.99990 | -0.9979 | -0.003.21
  return CAP_get_checknodes(0.0,1.0,-1.0, 0.0,0.01)
end

function CAP_get_C_values(CAP_raw, checknodes)
    Cx_values = []
    #println(CAP_raw.U)
    #println(CAP_raw.C)
    
    start_i = 1
    for row = 1:size(checknodes,1)
        direction = 1
        Ux = checknodes[row,1]
        xk = -1
        i = start_i
        while i < length(CAP_raw.U)
            #println(Ux, " // ", i, " // [", CAP_raw.U[i], ", ", CAP_raw.U[i+1], "]")
            if direction*(CAP_raw.U[i+1] - CAP_raw.U[i]) < 0
                direction *= -1
            end
            if is_between(Ux, CAP_raw.U[i], CAP_raw.U[i+1]) & (direction == checknodes[row,2]) 
                xk = i
                start_i = i
                break
            end
            i += 1
        end
        
        if xk == -1
            println("Error: CAP_get_C_values: Ux ",Ux," is not in CAP_raw.U: ",CAP_raw.U)
            return Exception()
        end
        Ck = CAP_raw.C[xk]
        Cl = CAP_raw.C[xk+1]
        Cx = Ck + ((Cl - Ck)/(CAP_raw.U[xk+1] - CAP_raw.U[xk])) * (Ux - CAP_raw.U[xk])
        
        Cx_values = [Cx_values; Cx]
        end
    
    return Cx_values
end

function apply_checknodes(sim::CAP_simulation, CAP_in, checknodes)
    if checknodes == Nothing
      return CAP_in
    else
      return DataFrame( U = checknodes[:,1], C = CAP_get_C_values(CAP_in, checknodes))
    end
end


#### plotting stuff ####
function setting_legend(SIM::CAP_simulation; latex=true)
  if latex
    return "\$\\theta=$(SIM.TC)\$°C \$\\mathrm{O}_2=$(SIM.pO2)\\%\$"
  else
    return "TC=$(SIM.TC)°C pO2=$(SIM.pO2)%"
  end
end

function typical_plot_sim(SIM::CAP_simulation, CAP_df, additional_string="", to_standard_figure=true)
  if to_standard_figure
    figure(SIM.fig_num, figsize=SIM.fig_size)
  end
  if SIM.analytical
    my_label = "ANL $(setting_legend(SIM))$(additional_string)"
  else
    my_label = "sim $(setting_legend(SIM))$(additional_string)"
  end
  
  PyPlot.title("Differential Capacitance")
  PyPlot.xlabel(L"\eta \ [V]")
  PyPlot.ylabel(L"C \ [F]")
  
  PyPlot.plot(CAP_df.U, CAP_df.C, label=my_label)
  PyPlot.grid(true)
  if !(my_label == "")
      legend(loc="best")
  end
end

function typical_plot_exp(SIM::CAP_simulation, CAP_df, additional_string="", to_standard_figure=true) 
  if to_standard_figure
    figure(SIM.fig_num, figsize=SIM.fig_size)
  end
  my_label = "exp $(setting_legend(SIM))$(additional_string)"
  
  PyPlot.title("Differential Capacitance")
  PyPlot.xlabel(L"\eta \ [V]")
  PyPlot.ylabel(L"C \ [F]")
  
  PyPlot.plot(CAP_df.U, CAP_df.C, "--", label=my_label)
  PyPlot.grid(true)

  if !(my_label == "")
      legend(loc="best")
  end
end

