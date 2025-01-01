# using Printf

const kB = 1.3806488e-23
const N_A = 6.02214129e23
const R = kB*N_A
const e0 = 1.602176565e-19
const eV = e0


mutable struct prms_struct
  names
  values
end

function prepare_prms(mask, x0, x)
    prms = []
    xi = 1
    for i in collect(1 : 1 :length(mask))
        if convert(Bool,mask[i])
            append!(prms, x[xi])
            xi += 1
        else
            append!(prms, x0[i])
        end
    end
    return Tuple(prms)
end

function consistency_check(prms::prms_struct)
  if prms.names == Nothing && prms.values == Nothing
    return true
  end
  
  if size(prms.names,1) != size(prms.values,1)
    println("ERROR: consistency_check: shape mismatch size(prms.names,1) != size(prms.values,1) >>>> $(size(prms.names,1)) != $(size(prms.values,1)) ")
    return throw(Exception)
  end
  
  for i in 1:size(prms.values,1)
    if typeof(prms.values[i]) != Colon
      if size(prms.values[i],1) < 1
        println("ERROR: consistency_check: empty_field prms.values[$i]")
        return throw(Exception)
      end
    end
  end
  
  return true
end

function for_each_prms_in_prms_lists(prms_lists, perform_generic)
  function recursive_call(output_set, active_idx)
    if active_idx > size(prms_lists,1)
      perform_generic(output_set)
    else
      for parameter in prms_lists[active_idx]
        recursive_call(push!(deepcopy(output_set),parameter), active_idx + 1)
      end
    end
  end
  recursive_call([],1)
  return
end

function for_each_indicies_in_prms_lists(prms_lists, perform_generic)
  function recursive_call(output_set, active_idx)
    if active_idx > size(prms_lists,1)
      perform_generic(output_set)
    else
      for i in 1:size(prms_lists[active_idx],1)
        recursive_call(push!(deepcopy(output_set),i), active_idx + 1)
      end
    end
  end
  recursive_call([],1)
  return
end

function get_prms_from_indicies(prms_lists, tuple)
  output_array = []
  for (i,list) in enumerate(prms_lists)
    append!(output_array, list[tuple[i]])
  end
  return output_array
end

function check_equal_size(list_1, list_2)
  if length(list_1) == length(list_2)
    return true
  else
    println("ERROR: check_equal_size: shape mismatch $(list_1) and $(list_2)")
    return throw(Exception)
  end
end

function check_x_in(x, low, upp)
  for (i, item) in enumerate(x)
    if item < low[i] || upp[i] < item
      return false
    end
  end
  return true
end

function take_only_masked(mask, array)
  output = []
  for (i, item) in enumerate(mask)
    if item == 1
      append!(output, [array[i]])
    end
  end
  return output
end


function make_array_from_string(input)
  if typeof(input) == String
    return [input]
  else
    return input
  end
end

export around
function around(A, bias, count)
  return collect(A-bias : 2.0*bias/(count-1) : A+bias)
end
       

function printfields(this)
    for name in fieldnames(typeof(this))
        @printf("%8s = ",name)
        try
          println(getfield(this,name))
        catch
          println(" <<< UNDEFINED !!! >>>")
        end
    end
end

function convert_prms_values_to_LaTeX_format(;
  prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
  prms_values,
  prms_names_output=["kin", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "OC", "ms_par", "e_fac"])
  
  function etk(exp_string)
    exp_id = findall(x -> x==exp_string, prms_names)
    
    return Bool(prms_values[exp_id]...) ? "E" : "L"
  end
  
  found_item = false
  output_string = ""
  for name_out in prms_names_output
    found_item = false
    for (i, name_in) in enumerate(prms_names)
      if name_out == name_in
        output_string *= "$(@sprintf("%.3f",prms_values[i])) & "
        found_item = true
        break
      end
    end
    if name_out == "kin"
      output_string *= "$(etk("A.exp"))$(etk("R.exp"))$(etk("O.exp")) & "
      found_item = true
    end
    
    if !found_item
      println("ERROR: name_out $(name_out) not found!")
      throw(Exception)
    end
  end
  if output_string == ""
    return output_string
  else
    return output_string[1:end-3]
  end
end


# Temperature parametrization
# function Arrhenius_template(T, E, p)
#   return p*exp(-E/(1.380649e-23 * T))
# end

function quadratic_template(T, A, B, C)
  return A*((T-973.15)/50)^2 + B*((T-973.15)/50) + C
end


function EIS_get_checknodes_geometrical(start_n, end_n, n_fac)
    # frequency nodes to compare
    # must be sorted upwards
    w_list = zeros(0)
    w = start_n
    if start_n == end_n
            append!(w_list, start_n)
            return w_list
    end
    if n_fac <= 1 || start_n > end_n
        println("ERROR: get_checknodes: n_fac <= 1 || start_n > end_n")
        return throw(Exception)
    end
    while w < end_n
        append!(w_list, w)
        w *= n_fac
    end    
    return w_list
end


##### fiting parameters
# supposing f = f(x, p) and the question is about how many parameters is in array "p"
function findout_number_parameters(f, max_number_of_parameters=100)
  for i in 1:max_number_of_parameters
      try 
          p = zeros(i)
          f(1,p)
          return i
      catch e
          if typeof(e) != BoundsError
              throw(e)
          end
      end
  end
end

function fit_func_prms(xdata, ydata, func; 
                      boundary_prms=nothing,
                      initial_prms=nothing,
                      weights=nothing,
                      plot_bool=false, max_number_of_parameters=100, penalty_constant = 1000
                      )
  n_p = findout_number_parameters(func, max_number_of_parameters)

  if length(xdata) != length(ydata)
      println("ERROR: length(xdata) != length(ydata), $(length(xdata)) != $(length(ydata))")
      return throw(Exception)
  end
  if typeof(initial_prms) == Nothing
      initial_prms = zeros(n_p)
      initial_prms .= 1.0
  end
  if typeof(boundary_prms) == Nothing
      boundary_prms = Vector{Vector}(undef, n_p)
      [boundary_prms[i] = [-Inf, Inf] for i in 1:n_p]
  end
  if typeof(weights) == Nothing
    weights = Vector{Float64}(undef, length(xdata))
    weights .= 1.0
  else
    if length(weights) != length(xdata)
      println("ERROR: length(weights) != length(xdata), $(length(weights)) != $(length(xdata))")
    end
  end

  function in_bounds(p)
      for i in 1:n_p
          if !((boundary_prms[i][1] <= p[i]) && (p[i] <= boundary_prms[i][2]))
              return false
          end
      end
      return true
  end

  if !in_bounds(initial_prms)
      println("ERROR: initial_prms $(initial_prms) not in bounds $(boundary_prms)")
      return throw(Exception)
  end
  
  function to_optimize(p)
      sum = 0
      if !in_bounds(p)
          return sum += penalty_constant
      end
     
      for i in 1:length(xdata)
          sum += weights[i]*(func(xdata[i], p) - ydata[i])^2
      end
      #println(p, " -> ", sum)
      return sum
  end

  opt = optimize(to_optimize, initial_prms)

  if plot_bool
      plot(xdata, ydata, label="data", "-x")
      plot(xdata, [func(x, initial_prms) for x in xdata], label="init", ":")
      plot(xdata, [func(x, opt.minimizer) for x in xdata], label="fit", "--")
      legend()
  end
  return opt.minimizer
end

function get_f_range_from_triplet(tri)
  f_storage = Float64[]
  f = tri[1]
  while f <= tri[2]
    push!(f_storage, f)
    f = f*tri[3]
  end  
  return f_storage
end