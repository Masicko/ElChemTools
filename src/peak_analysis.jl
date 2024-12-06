# general function for fitting 
#   fit_func_prms( (x, p) -> f(x,p), initial_prms, plot_results)
#   ---> find parametr count
#   ---> 
#
#
# find_gaussian_peaks(h_tau, tau) -> returns datastructure gaussian_fit => function, prms, number_of_peaks
#
#   TYJO ... fakt nevim, jak to rozumne charakterizovat. Jak urcit, ze jeden inflexni bod je ten, co chci (mezi kopcama) 
#   a jiny inflexni bod nechci
# 
#   ... jo, slo by to tak, ze se budu divat fakt na krivost a peaky najdu jako negativni peaky krivosti. 
#       A udoli muzu najit bud jako pozitivni peaky krivosti.. a nebo nejaky prumer mezi inflexnimi body kladne casti krivosti
#       <> Mozna to jeste nejak normovat, protoze male peaky budou mit mozna mensi krivost nez velke peaky... ale mozna je to blbost
#
#
# TODO TODO ... eis_100 ... je tam jeden peak, ktery oci vidi, ale algoritmus ne.




@kwdef mutable struct gaussian_fit_struct
    xdata :: Vector
    ydata :: Vector
    num_peaks :: Int64
    division_idxs :: Vector
    #
    func = Nothing
    prms = Nothing
end

function gaussian_peak(x, mu, sigma, h, tilt=0.0)
    return h*exp( - (((x - mu)*(1 + tilt*sign(x-mu)))^2 ) / (2*sigma^2)    )
end

function my_deriv(a)
    return [a[i+1] - a[i] for i in 1:length(a)-1]
end

function find_divisions_std(h_tau)
    tol = maximum(h_tau)*1e-8
    dht = my_deriv(h_tau)
    ddht = my_deriv(dht)
    div_idxs = [[1, "0"]]
    i = 1 

    while i <= length(dht)-2
        #@show i
        
        # local minima
        # if (dht[i] + tol < 0.0) && (0.0 < dht[i+1] - tol)
        #     push!(div_idxs, [i+1, "1"])
        #     i += 2
        # end

        # inflex point
        # if (dht[i] - tol > dht[i+1]) && (dht[i+1] < dht[i+2] - tol)
        #     push!(div_idxs, [i+1, "2"])
        #     i += 2
        # end

        # between inflex points
        if (dht[i] - tol > dht[i+1]) && (dht[i+1] < dht[i+2] - tol)
            push!(div_idxs, [i+1, "2"])
            i += 2
        end

        i += 1
    end
    return div_idxs
end

function find_divisions_curv(h_tau)
    tol = maximum(h_tau)*1e-8
    dht = my_deriv(h_tau)
    ddht = my_deriv(dht)
    ddht = [1.0e-6, ddht..., 1.0e-6]
    div_idxs = [1]
    i = 1 

    previous_inflex_id = -1

    while i <= length(ddht)-1
        if ddht[i] <= 0 && ddht[i+1] >= 0
            previous_inflex_id = i
        end

        if ddht[i] >= 0 && ddht[i+1] <= 0 && previous_inflex_id > 0
            #@show " >> ", previous_inflex_id, i+1
            push!(div_idxs, Int64(round((i + 1 + previous_inflex_id)/2)))
            i += 1
            previous_inflex_id = -1
        end

        i += 1
    end

    push!(div_idxs, length(h_tau))
    return div_idxs
end

# the deciding points should be 1) local minima and 2) inflex points
function find_divisions(h_tau)
    find_divisions_curv(h_tau)
end

function plot_divisions(gf::gaussian_fit_struct)
    plot_ht_and_func(gf.ydata, nothing, omit_f=true)
    
    grid(true)
    if length(gf.division_idxs) != 0
        for idx in gf.division_idxs
            plot([idx, idx], [0, maximum(gf.ydata)], "black")
        end
    end
end

function plot_mimina(ht, div_idxs = nothing)
    figure(11) 
    plot(ht)
    #plot(collect(0.5 : 1.0 : length(ht)-1), ElChemTools.my_deriv(ht))
    plot(collect(1.0 : 1.0 : length(ht)-2), ElChemTools.my_deriv(ElChemTools.my_deriv(ht)))
    grid(true)
    if typeof(div_idxs) != Nothing
        for idx in div_idxs
            plot([idx-1, idx-1  ], [0, maximum(ht)], "black")
        end
    end
end

function num_of_prms_per_peak(gf::gaussian_fit_struct)
    return length(gf.prms)/gf.num_peaks
end

function plot_ht_and_func(ht, f; omit_ht=false, omit_f=false, style=nothing)
    xr = collect(1:length(ht))
    if !omit_ht
        plot(xr, ht)
    end
    if !omit_f
        if typeof(style) == Nothing
            plot(xr, [f(x) for x in xr])
        else
            plot(xr, [f(x) for x in xr], style)
        end
    end
    return
end

function find_peak_idxs(ddht, l_idx, r_idx)
    peak_idxs = []
    i = l_idx
    while i < r_idx
        if ddht[i] < 0 && ddht[i+1] > 0
            push!(peak_idxs, i)
        end
        if ddht[i] > 0 && ddht[i+1] < 0
            push!(peak_idxs, i+1)
            i += 1
        end
        i += 1
    end
    #@show peak_idxs
    if length(peak_idxs) != 2
        println("ERROR: peak_idxs $(peak_idxs) doesnt have 2 elements")
        return throw(Exception)
    end
    return peak_idxs
end

function fit_gaussian_peaks(gf::gaussian_fit_struct, func::Function)
    init_prms_tot = Float64[]
    bounds_tot = []

    plot_divisions(gf)

    #for i in 1:3
    for i in 1:gf.num_peaks
        init_prms_loc = Float64[]
        l_idx = gf.division_idxs[i]
        r_idx = gf.division_idxs[i+1]
        
        ddht = my_deriv(my_deriv(gf.ydata))
        ddht = [1.0e-5, ddht..., 1.0e-5]
        plot(collect(1:length(ddht)), ddht)
        peak_idxs = find_peak_idxs(ddht, l_idx, r_idx)

        min = Inf
        idx_min = -1
        for k in peak_idxs[1] : peak_idxs[2]
            if ddht[k] < min
                min = ddht[k]
                idx_min = k
            end
        end
        mean = idx_min

        # mean
        push!(init_prms_loc, (mean))
        push!(bounds_tot, [l_idx, r_idx])

        # var
        push!(init_prms_loc, 0.4*(peak_idxs[2] - peak_idxs[1]))
        push!(bounds_tot, [1.0e-7, 1.0e4])

        # h
        h_at_mean = gf.ydata[mean]
        push!(init_prms_loc, h_at_mean)
        push!(bounds_tot, h_at_mean.*[0.5, 1.5])

        # tilt
        if findout_number_parameters(func) > 3
            push!(init_prms_loc, 0.0)
            push!(bounds_tot, [-1.0, 1.0])
        end


        # local fit
        #plot_ht_and_func(gf.ydata, x -> func(x, init_prms_loc), omit_ht=true, style=":")
      
        #loc_ydata = gf.ydata[l_idx : r_idx]
        loc_ydata = gf.ydata[peak_idxs[1] : peak_idxs[2]]
        
        fitted_prms_loc = fit_func_prms(
            collect(peak_idxs[1] : peak_idxs[2]), 
            loc_ydata, 
            func, 
            initial_prms = init_prms_loc,
            boundary_prms = bounds_tot[end - length(init_prms_loc)+1 : end],
            #plot_bool=true
        )

        append!(init_prms_tot, fitted_prms_loc)

        plot_ht_and_func(gf.ydata, x -> func(x, fitted_prms_loc), omit_ht=true, style=":")
    end

    n_prms_per_peak = Int64(length(init_prms_tot)/gf.num_peaks)
    func_tot(x , p) = sum([func(x, p[i:i+n_prms_per_peak-1])
                            for i in 1: n_prms_per_peak : n_prms_per_peak*gf.num_peaks ]
                        )

    #plot_ht_and_func(gf.ydata, x -> func_tot(x, init_prms_tot), omit_ht=true, style="-b")

    fitted_prms_tot = fit_func_prms(
            collect(1 : length(gf.ydata)), 
            gf.ydata, 
            func_tot, 
            initial_prms = init_prms_tot,
            boundary_prms = bounds_tot,
        )

    plot_ht_and_func(gf.ydata, x -> func_tot(x, fitted_prms_tot), omit_ht=true, style="-.")
    
    #@show fitted_prms_tot

    for i in 1: n_prms_per_peak : n_prms_per_peak*gf.num_peaks
        plot_ht_and_func(gf.ydata, x -> func(x, fitted_prms_tot[i : i+n_prms_per_peak-1 ]), omit_ht=true, style="--")
    end

    return
end

function find_gaussian_peaks(tau, h_tau, control::DRT_control_struct=DRT_control_struct())
    division_idxs = find_divisions(h_tau)
    #plot_mimina(h_tau, division_idxs)
    l_tau = log.(10, tau)

    gaussian_fit = gaussian_fit_struct(
            xdata = l_tau, 
            ydata = h_tau, 
            num_peaks = length(division_idxs) - 1,
            division_idxs = division_idxs
    )

    figure(11)    
    #fit_gaussian_peaks(gaussian_fit, (x, p) -> gaussian_peak(x, p[1], p[2], p[3]))
    fit_gaussian_peaks(gaussian_fit, (x, p) -> gaussian_peak(x, p[1], p[2], p[3], p[4]))

    #part_fits, part_bounds = get_partial_data(tau, h_tau, division_idxs)
    #whole fit
    return gaussian_fit
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

    # @show division_starting_idx, ending_idx
    ht_o = DRT.h[division_starting_idx : ending_idx]
    t_o = DRT.tau_range[division_starting_idx : ending_idx]
    
    serialize("last_lt_ht.dat", [log.(10, t_o), ht_o])
    find_gaussian_peaks(t_o, ht_o, DRT.control)


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


    #@show division_idxs, [log(10, DRT.tau_range[idx]) for idx in division_idxs]
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
    #@show R_list
    DRT.R_peak_list = R_list
    return
end
  