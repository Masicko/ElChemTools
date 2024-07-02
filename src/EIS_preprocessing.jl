using PyPlot

##
##  Notes:
#### - Potrebuju nejake spolecne veci pro procesing, ktere by se mely zavolat hned ze startu .. a pak predavat pripadne hloubeji?
#### - Ukladani obrazku budu delat pak mimo
#### - nejdriv zde sepisu veskere veci pro automatiku. Az pak budu resit interaktivni veci.
#### - at existuje i "do nothing" varianta ... konstruktor control tridy
##
## EIS_preprocessing(EIS_df, control... )
##    - f_interval_trim(EIS_df, control)         d=
##    - dropping specified points( 50 hz)
##    - trim inductance( ) -> {zpetny ocas, vse, na zaklade variace?}
##    - trim low_frequency_points -> {zpetny ocas, vse, na zaklade variace?}
##    - dropping really bad points(EIS_df, control) (vysoka variace?)
##    - convolution(EIS_df, control)
##    - 



mutable struct EIS_preprocessing_control
  f_interval
  scale_factor
  add_inductance
  trim_inductance
  subtract_inductance_fac
  leave_only_inductance_points
  outlayers_threshold
  use_DRT
  DRT_control
  DAN_prepro
  #
  debug_plot
  
  EIS_preprocessing_control(d) = new()
end

function EIS_preprocessing_control(;
                                    f_interval=Nothing, 
                                    add_inductance=0, 
                                    scale_factor=Nothing,
                                    trim_inductance=false,
                                    subtract_inductance_fac=0.0,
                                    leave_only_inductance_points=Nothing,
                                    outlayers_threshold=1.5,                                    
                                    use_DRT=false, DRT_control=DRT_control_struct(),
                                    DAN_prepro=false,
                                    #
                                    debug_plot=false
                                    )
  this = EIS_preprocessing_control(1)
  this.f_interval = f_interval
  this.scale_factor = scale_factor
  this.add_inductance = add_inductance
  this.trim_inductance = trim_inductance
  this.subtract_inductance_fac = subtract_inductance_fac
  this.leave_only_inductance_points = leave_only_inductance_points
  this.outlayers_threshold = outlayers_threshold
  this.use_DRT = use_DRT
  this.DRT_control = DRT_control
  this.DAN_prepro = DAN_prepro
  #
  this.debug_plot = debug_plot
  return this
end

# procedure CODING
# - f_interval

##

function EIS_crop_to_f_interval(EIS_exp, f_interval)
  return filter(row -> (f_interval[1] < row.f) && (row.f < f_interval[2]), EIS_exp)
end

function EIS_filter_to_idx_blacklist(EIS_exp, idx_list)
  out_df = DataFrame(f=[], Z=[])
  for i in 1:(size(EIS_exp)[1])
    if !(i in idx_list)
      push!(out_df, (EIS_exp.f[i], EIS_exp.Z[i]))
    end
  end
  return out_df
end

##

function get_lowest_freq_idx(EIS_df; find_at_least_negative=10)
  lowest_freq_idx = -1
  negative_counter = 0
  for (i, Z) in enumerate(EIS_df.Z)
    if imag(Z) < 0
      negative_counter += 1
    else
      negative_counter = 0
    end
    if negative_counter == find_at_least_negative
      lowest_freq_idx = i - find_at_least_negative + 1
      break
    end
  end
  return lowest_freq_idx
end

function get_high_freq_intersection_idx(EIS_df, find_at_least_negative=5)
  HF_intersection_idx = -1
  negative_counter = 0
  for (i, Z) in enumerate(reverse(EIS_df.Z))
    if imag(Z) < 0
      negative_counter += 1
    else
      negative_counter = 0
    end
    if negative_counter == find_at_least_negative
      HF_intersection_idx = length(EIS_df.Z) - (i - find_at_least_negative)
      break
    end
  end
  return HF_intersection_idx
end

#####################
#####################
# old version
function f_interval_auto_processing(EIS_df, trim_inductance=false, leave_only_inductance_points=false)
  
  function get_lowest_freq_idx(EIS_df; find_at_least_negative=10)
    lowest_freq_idx = -1
    negative_counter = 0
    for (i, Z) in enumerate(EIS_df.Z)
      if imag(Z) < 0
        negative_counter += 1
      else
        negative_counter = 0
      end
      if negative_counter == find_at_least_negative
        lowest_freq_idx = i - find_at_least_negative + 1
        break
      end
    end
    return lowest_freq_idx
  end
  
  
  #ElChemTools.typical_plot_exp(EIS_simulation(800, 100, 0.0, use_DRT=false)[1], EIS_df, "!f_interval_auto_processing")

  
  #lowest frequency cut off
  lowest_freq_idx = get_lowest_freq_idx(EIS_df, find_at_least_negative=10)
  if lowest_freq_idx == -1
    lowest_freq_idx = get_lowest_freq_idx(EIS_df, find_at_least_negative=4)
  end
  if lowest_freq_idx == -1
    println("ERROR: lowest_freq_idx not found!")
    
    lowest_freq_idx=1
    #return throw(Exception)
  end
  
  # intersection with x axis
  x_intersection_freq_idx = -1
  positive_counter = 0
  for i in (lowest_freq_idx + 5):length(EIS_df.f)
    if imag(EIS_df.Z[i]) > 0
      positive_counter += 1
    else  
      positive_counter = 0
    end
    if positive_counter == 8
      x_intersection_freq_idx = i - 7
      break
    end
  end

  #@show x_intersection_freq_idx
  #@show lowest_freq_idx
  
  if x_intersection_freq_idx == -1
    return DataFrame(f = EIS_df.f[lowest_freq_idx:end], Z = EIS_df.Z[lowest_freq_idx:end])
  else
    # (non)-inductance cut off
    if trim_inductance
      highest_freq_idx = x_intersection_freq_idx
    elseif leave_only_inductance_points!=Nothing   
      lowest_freq_idx = x_intersection_freq_idx
      highest_freq_idx = min(lowest_freq_idx + leave_only_inductance_points - 1, length(EIS_df.f))
      return DataFrame(f = EIS_df.f[lowest_freq_idx:highest_freq_idx], Z = EIS_df.Z[lowest_freq_idx:highest_freq_idx])
    else
      accepted_inductance_real_axis_threshold = 0.00*real(EIS_df.Z[lowest_freq_idx]) + 1.00*real(EIS_df.Z[x_intersection_freq_idx])
      highest_freq_idx = -1
      for i in (x_intersection_freq_idx + 1 ):length(EIS_df.f)
        if real(EIS_df.Z[i]) > accepted_inductance_real_axis_threshold
          highest_freq_idx = i 
          break
        end
      end
      if highest_freq_idx == -1
        return DataFrame(f = EIS_df.f[lowest_freq_idx:end], Z = EIS_df.Z[lowest_freq_idx:end])
      end
    end
    return DataFrame(f = EIS_df.f[lowest_freq_idx:highest_freq_idx], Z = EIS_df.Z[lowest_freq_idx:highest_freq_idx])
  end
end
##################
##################


function f_interval_trim(EIS_exp, EIS_preprocessing_control)
  f_interval = EIS_preprocessing_control.f_interval  
  if f_interval!=Nothing
    if f_interval == "auto"      
      #typical_plot_exp(SIM, EIS_exp, "! before")
      EIS_exp = f_interval_auto_processing( EIS_exp, 
                                            EIS_preprocessing_control.trim_inductance,
                                            EIS_preprocessing_control.leave_only_inductance_points
                                          )               
      #typical_plot_exp(SIM, EIS_exp, "! after")
    else
      EIS_exp = EIS_crop_to_f_interval(EIS_exp, f_interval)      
    end
  end  
  return EIS_exp
end


function EIS_add_inductance!(new_EIS_df::DataFrame, EIS_preprocessing_control::EIS_preprocessing_control)    
  L = EIS_preprocessing_control.add_inductance
  for (i, f) in enumerate(new_EIS_df.f)
      new_EIS_df.Z[i] += im*2*pi*f*L
  end
end

function get_average_Z_difference(Z_list; pyplot=true)
  average_Z_difference = [0.,0.]
  
  for (i, Z) in enumerate(Z_list[1:end-1])
    average_Z_difference[1] += abs(real(Z_list[i] - Z_list[i+1]))
    average_Z_difference[2] += abs(imag(Z_list[i] - Z_list[i+1]))
  end  
  average_Z_difference ./= length(Z_list) - 1  
  return average_Z_difference
end

function dropping_really_bad_points(EIS_df, EIS_preprocessing_control)
  average_Z_difference = get_average_Z_difference(EIS_df.Z)
  threshold = EIS_preprocessing_control.outlayers_threshold
  frequency_blacklist = []
  
  pyplot = EIS_preprocessing_control.debug_plot
  fig_num=135
  
  diff_storage = [[],[]]
  for (i, Z) in enumerate(EIS_df.Z[1:end-1])    
    push!(diff_storage[1], real(EIS_df.Z[i] - EIS_df.Z[i+1]))
    push!(diff_storage[2], imag(EIS_df.Z[i] - EIS_df.Z[i+1]))      
  end
  
  function add_to_unique_list(list, idx)    
    if length(list) < 1 || list[end] != idx
      push!(list, idx)
    end
  end
  
  function diff_storage_outlayer_test(idx, i)  
#     @show diff_storage[idx][i]
#     @show threshold*average_Z_difference[idx]
    if diff_storage[idx][i] > threshold*average_Z_difference[idx]
      if diff_storage[idx][i+1] < -threshold*average_Z_difference[idx]
        # indivitual point        
        add_to_unique_list(frequency_blacklist, i + 1)
      else        
        # shift of the whole segment
      end
    end
    
    if diff_storage[idx][i] < -threshold*average_Z_difference[idx]
      if diff_storage[idx][i+1] > threshold*average_Z_difference[idx]
        # indivitual point
        add_to_unique_list(frequency_blacklist, i + 1)
      else
        # shift of the whole segment
      end
    end    
  end
  
  for i in 1:length(diff_storage[1]) - 1    
    diff_storage_outlayer_test(1, i)
    diff_storage_outlayer_test(2, i)
  end
    
  
  if pyplot
    @show frequency_blacklist
    
    figure(fig_num)   
    plot(diff_storage[1])
    plot(diff_storage[2])
    #
    threshold_constant = deepcopy(diff_storage)
    threshold_constant[1] .= threshold*average_Z_difference[1]
    threshold_constant[2] .= threshold*average_Z_difference[2]
    plot(threshold_constant[1])
    plot(threshold_constant[2])
  end
  
  return EIS_filter_to_idx_blacklist(EIS_df, frequency_blacklist)
end

function convolution(EIS_df, EIS_preprocessing_control)
  if EIS_preprocessing_control.use_DRT
    DRT_actual = get_DRT(EIS_df, EIS_preprocessing_control.DRT_control)
    return deepcopy(DRT_actual.EIS_df)
  else
    return EIS_df
  end
end

function subtract_inductance!(EIS_df, EIS_preprocessing_control)
  if EIS_preprocessing_control.subtract_inductance_fac > 0
    DRT_actual = get_DRT(EIS_df, EIS_preprocessing_control.DRT_control)
    L = DRT_actual.L*EIS_preprocessing_control.subtract_inductance_fac
    
    for (i, f) in enumerate(EIS_df.f)
        EIS_df.Z[i] -= im*2*pi*f*L
    end
  else
    return
  end
end

function EIS_scale_data(EIS_df, EIS_preprocessing_control)
  if EIS_preprocessing_control.scale_factor == Nothing
    return EIS_df
  else
    return DataFrame(f = EIS_df.f, Z = EIS_preprocessing_control.scale_factor .* EIS_df.Z)
  end
end

function DAN_preprocessing(EIS_df, EIS_preprocessing_control)
  
  if EIS_preprocessing_control.DAN_prepro
    # low frequency
   # max_R = -1.0
   # LF_idx = -1
   # for (i, Z) in enumerate(EIS_df.Z)
  #    if real(Z) > max_R
  #      LF_idx = i
 ##       max_R = real(Z)
 #     end
#    end
    LF_intersection_idx = get_lowest_freq_idx(EIS_df, find_at_least_negative = 5)
    if LF_intersection_idx == -1
      LF_intersection_idx = 1
    end


    # high frequency
    HF_intersection_idx = get_high_freq_intersection_idx(EIS_df, 5)
    if HF_intersection_idx == -1
      HF_intersection_idx = length(EIS_df.Z)
    end

    Z_list = [EIS_df.Z[LF_intersection_idx], EIS_df.Z..., EIS_df.Z[HF_intersection_idx]]
    f_list = [EIS_df.f[1]*0.9, EIS_df.f..., EIS_df.f[end]*1.1]
    
    EIS_df = DataFrame(f = f_list, Z = Z_list)
    return EIS_df
  else
    return EIS_df
  end
end


function EIS_preprocessing(EIS_df, EIS_preprocessing_control::EIS_preprocessing_control)            
    new_EIS_df = deepcopy(EIS_df)

    subtract_inductance!(new_EIS_df, EIS_preprocessing_control)
    
    new_EIS_df = f_interval_trim(new_EIS_df, EIS_preprocessing_control)    
    #
    new_EIS_df = EIS_scale_data(new_EIS_df, EIS_preprocessing_control)
    #
    #dropping_specified_points!(new_EIS_df, EIS_preprocessing_control)         --->> IMPORTANT
    #
    #assesing_chaos(new_EIS) #    ----> if there is even any useful piece of information?
    new_EIS_df = dropping_really_bad_points(new_EIS_df, EIS_preprocessing_control)      # TODO ->>> implement filtering whole segments out     
    #
    new_EIS_df = convolution(new_EIS_df, EIS_preprocessing_control)
    #
    new_EIS_df = DAN_preprocessing(new_EIS_df, EIS_preprocessing_control)
    #
    EIS_add_inductance!(new_EIS_df, EIS_preprocessing_control)
    
    
    return new_EIS_df
end
























