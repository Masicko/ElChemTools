# general function for fitting 
#   fit_func_prms( (x, p) -> f(x,p), initial_prms, plot_results)
#   ---> find parametr count
#   ---> 
#
#
# find_gaussian_peaks(h_tau, tau) -> returns datastructure gaussian_fit => function, prms, number_of_peaks

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
                        plot_bool=false, max_number_of_parameters=100
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
        #println(p)
        penalty_constant = 1000
        if !in_bounds(p)
            return penalty_constant
        end

        sum = 0
        for i in 1:length(xdata)
            sum += (func(xdata[i], p) - ydata[i])^2
        end
        return sum
    end

    opt = optimize(to_optimize, initial_prms)

    if plot_bool
        plot(xdata, ydata, label="data", "-x")
        plot(xdata, [func(x, opt.minimizer) for x in xdata], label="fit", "--")
        legend()
    end
    return opt.minimizer
end










function discrete_integrate(tau_range, h_tau)
    sum = 0
    for i in 1:length(tau_range)
    sum += h_tau[i]
    end
    return sum
end

function discrete_integrate(idx1, idx2, h_tau)
    sum = 0
    for i in idx1:idx2
    sum += h_tau[i]
    end
    return sum
end

function divide_and_evaluate_R_peaks(DRT::DRT_struct)
    
    function valley_check(idx, f)
    valley_th = threshold
    if (f[idx] >= f[idx+1] + valley_th) && (f[idx + 1] + valley_th <= f[idx + 2])
        return true 
    else
        return false
    end
    end

    function hill_start_check(idx, f)
    if (abs(f[idx] - f[idx+1]) <= threshold) && (f[idx + 1] + threshold <= f[idx + 2])
        return true 
    else
        return false
    end
    end
    
    function hill_end_check(idx, f)
    if (f[idx] > f[idx+1] + threshold) && (abs(f[idx + 1] - f[idx + 2]) <= threshold)
        return true 
    else
        return false
    end
    end
    
    threshold = maximum(DRT.h)/10000.0
    # first artefatic h decay
    division_starting_idx = 1
    while division_starting_idx + 1 <= length(DRT.h) && 
        DRT.h[division_starting_idx + 1] - DRT.h[division_starting_idx] < threshold
    division_starting_idx += 1
    end
    if division_starting_idx == length(DRT.h) - 1
    ## TODO
    end



    # ending and check for artefacts
    if length(DRT.h) >= 3
    ending_idx = length(DRT.h) - 2
    while ending_idx >= 1 
        if hill_end_check(ending_idx, DRT.h)
        ending_idx += 1
        break
        end
        ending_idx -=1
    end
    else
    ending_idx = length(DRT.h)
    end


    # main loop (with an assumption that in the low frequenceis there is no artefact)
    active_idx = division_starting_idx
    division_idxs = [division_starting_idx]
    while active_idx < ending_idx - 1
    if valley_check(active_idx, DRT.h) || hill_start_check(active_idx, DRT.h)
        #@show valley_check(active_idx, DRT.h), hill_start_check(active_idx, DRT.h)
        push!(division_idxs, active_idx+1)
    end
    active_idx += 1
    end
    push!(division_idxs, ending_idx)

    if ending_idx < division_starting_idx
    println("WARNING: starting idx is greater than ending idx in DRT!")
    division_idxs = [division_starting_idx, length(DRT.h)]
    return
    end


    @show division_idxs, [log(10, DRT.tau_range[idx]) for idx in division_idxs]
    # summing for R_list
    R_list = []
    for i in 1:length(division_idxs)-1
    R = discrete_integrate(division_idxs[i], division_idxs[i+1] - 1, DRT.h)
    h_max = -1.0
    f_max = -1.0
    for i in division_idxs[i] : division_idxs[i+1]
        if DRT.h[i] > h_max 
        h_max = DRT.h[i]
        f_max = 1/(DRT.tau_range[i] * 2*pi)
        end
    end
    
    #h_max_check = maximum(DRT.h[division_idxs[i] : division_idxs[i+1]])
    #@show h_max, h_max_check
    push!(R_list, (f_max, R))
    end
    @show R_list
    DRT.R_peak_list = R_list
    return
end
  