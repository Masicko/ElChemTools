# general function for fitting 
#   fit_func_prms( (x, p) -> f(x,p), initial_prms, plot_results)
#   ---> find parametr count
#   ---> 
#
#
# obtain_gaussian_peaks(h_tau, tau) -> returns datastructure gaussian_fit => function, prms, number_of_peaks
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
# ... rozdelovani na oblasti funguje pekne
# [ ] zlepsit inicialni nastrely ... specialne u "schovanych peaku" v udoli je to hrozne
# [ ] s tim souvisi nejake navazeni, jak je ktery peak dulezity pri fitovani. Sdelit algoritmu, ze ma "zacit" tema dulezityma.. 
#     a pak dofitovat ty ostatni mensi peaky, ktere nejsou tak zretelne. Neco jako odcitani kopecku 
#     (vlastne by bylo zajimave to alespon videt.. co se stane, kdyz se odectou ty dva velke jasne kopce, co zustane?)




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

function find_divisions_maximum(h_tau)
    tol = maximum(h_tau)*1e-8
    ddht = my_2nd_der(h_tau)

    div_idxs = []
    i = 2
    while i <= length(ddht)-1
        if ddht[i-1] + tol < ddht[i] && ddht[i] > ddht[i+1] + tol
            push!(div_idxs, i)
        end
        i += 1
    end

    if length(div_idxs) < 2
        println("ERROR: div_idxs < 2")
        return throw(Exception)
    end

    div_idxs = [1, div_idxs[2:end-1]..., length(h_tau)]
    #@show div_idxs
    return div_idxs  
end

function my_2nd_der(a)
    res = my_deriv(my_deriv(a))
    res = [res[1]/2.0, res..., res[2]/2.0]
    res = map(x -> ( abs(x) < 1.0e-10 ? 0.0 : x), res)
    return res
end

function find_divisions_curv(h_tau)
    ddht = my_2nd_der(h_tau)
    div_idxs = [1]
    i = 1 

    previous_inflex_id = -1

    while i <= length(ddht)-1
        if ddht[i] < 0 && ddht[i+1] > 0
            previous_inflex_id = i
        end

        if ddht[i] > 0 && ddht[i+1] < 0 && previous_inflex_id > 0
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
    find_divisions_maximum(h_tau)
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

function plot_divisions(ht, division_idxs)
    plot_ht_and_func(ht, nothing, omit_f=true)
    
    grid(true)
    if length(division_idxs) != 0
        for idx in division_idxs
            plot([idx, idx], [0, maximum(ht)], "black")
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
    if length(peak_idxs) == 1
        push!(peak_idxs, peak_idxs[1])
    end
    #@show peak_idxs
    if length(peak_idxs) == 0
        println("ERROR: peak_idxs $(peak_idxs) do not exist")
        return throw(Exception)
    end
    return peak_idxs
end

function no_peak_prms(aux_gf)
    return [0.0, 0.0, 1.0]
end

function change_prms_of_peak(gf, prms, peak)
    npp = num_of_prms_per_peak(gf) 
    if peak > gf.num_peaks
        println("ERROR: peak $(peak) > num_peaks $(gf.num_peaks)")
        throw(Exception)
    end
    gf.prms[Int64.((peak -1)*npp + 1 : peak*npp)] = prms
end

function fit_gaussian_peaks(gf::gaussian_fit_struct, func::Function; plot_bool = false)
    init_prms_tot = Float64[]
    bounds_tot = []

    #for i in 1:3
    for i in 1:gf.num_peaks
        init_prms_loc = Float64[]
        l_idx = gf.division_idxs[i]
        r_idx = gf.division_idxs[i+1]
        
        ddht = my_2nd_der(gf.ydata)
        
        plot_bool && plot(collect(1:length(ddht)), ddht)

        the_min = Inf
        idx_min = -1
        for k in l_idx : r_idx
            if ddht[k] < the_min
                the_min = ddht[k]
                idx_min = k
            end
        end
        mean = idx_min
        min_peak_radius = 2
        prop_peak_idxs = find_peak_idxs(ddht, l_idx, r_idx)
        min_prop_peak_idxs = [
            min(mean - min_peak_radius, prop_peak_idxs[1]), 
            max(mean + min_peak_radius, prop_peak_idxs[2])
        ]
        peak_idxs = [max(min_prop_peak_idxs[1], l_idx), min(min_prop_peak_idxs[2], r_idx)]
        

        # meanx
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

        #plot_ht_and_func(gf.ydata, x -> func(x, fitted_prms_loc), omit_ht=true, style=":")
    end

    n_prms_per_peak = Int64(length(init_prms_tot)/gf.num_peaks)
    func_tot(x , p) = sum([func(x, p[i:i+n_prms_per_peak-1])
                            for i in 1: n_prms_per_peak : n_prms_per_peak*gf.num_peaks ]
                        )

    #plot_ht_and_func(gf.ydata, x -> func_tot(x, init_prms_tot), omit_ht=true, style="-b")

    weights = zeros(length(gf.ydata))
    weights .= 1.0
    #weights[gf.division_idxs[2] : gf.division_idxs[3]] .= 0.3

    fitted_prms_tot = fit_func_prms(
            collect(1 : length(gf.ydata)), 
            gf.ydata, 
            func_tot, 
            initial_prms = init_prms_tot,
            boundary_prms = bounds_tot,
            weights = weights
        )
    gf.prms = fitted_prms_tot
    gf.func = func_tot

    plot_bool && plot_ht_and_func(gf.ydata, x -> func_tot(x, fitted_prms_tot), omit_ht=true, style="-.")
    
    #@show fitted_prms_tot

    plot_bool && for i in 1: n_prms_per_peak : n_prms_per_peak*gf.num_peaks
        plot_ht_and_func(gf.ydata, x -> func(x, fitted_prms_tot[i : i+n_prms_per_peak-1 ]), omit_ht=true, style="--")
    end

    return
end


function subtract_peaks(gf::gaussian_fit_struct, sub_peaks)
    aux_gf = deepcopy(gf)
    for i in 1:aux_gf.num_peaks
        !(i in sub_peaks) && change_prms_of_peak(aux_gf, no_peak_prms(aux_gf), i)
    end

    
    res = [aux_gf.ydata[i] - aux_gf.func(i, aux_gf.prms) for i in 1:length(aux_gf.xdata)]
    plot(1:length(res), res, "black") 
    return 
end

function get_Rt_info(gf::gaussian_fit_struct)
    peaks_funcs = []
    for peak_num in 1:gf.num_peaks
        aux_gf = deepcopy(gf)
        for i in 1:gf.num_peaks
            i != peak_num && change_prms_of_peak(aux_gf, no_peak_prms(aux_gf), i)
        end
        res = [aux_gf.func(i, aux_gf.prms) for i in 1:length(aux_gf.xdata)]
        #plot(1:length(res), res, "black")
        push!(peaks_funcs, res)
    end
    
    peak_Rs = zeros(gf.num_peaks)
    for x in 1:length(gf.xdata)
        sum_on_x = 0.0
        for peak_num in 1:gf.num_peaks
            sum_on_x += peaks_funcs[peak_num][x]
        end
        for peak_num in 1:gf.num_peaks
            peak_Rs[peak_num] += gf.ydata[x] * (peaks_funcs[peak_num][x]/sum_on_x)
        end
    end

    peak_ts = []
    npp = num_of_prms_per_peak(gf)
    for peak_num in 1:gf.num_peaks
        mean = gf.prms[Int((peak_num -1)*npp + 1)]
        mean_tau = exponential_interpolation(mean, 1, length(gf.xdata), 10.0^gf.xdata[1], 10.0^gf.xdata[end])        
        push!(peak_ts, mean_tau)
    end
    return [z for z in zip(1.0 ./(2 .*pi.* peak_ts), peak_Rs)]
end

function obtain_gaussian_peaks(tau, h_tau, division_idxs, control::DRT_control_struct=DRT_control_struct())
    l_tau = log.(10, tau)
    gaussian_fit = gaussian_fit_struct(
            xdata = l_tau, 
            ydata = h_tau, 
            num_peaks = length(division_idxs) - 1,
            division_idxs = division_idxs
    )

    # without tilt
    fit_gaussian_peaks(gaussian_fit, (x, p) -> gaussian_peak(x, p[1], p[2], p[3]), plot_bool=control.R_peaks_plot_bool)
    # with tilt
    #fit_gaussian_peaks(gaussian_fit, (x, p) -> gaussian_peak(x, p[1], p[2], p[3], p[4]))

    #subtract_peaks(gaussian_fit, [1,3])
    Rt_res = get_Rt_info(gaussian_fit)

    return Rt_res
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


function valley_check(idx, f, threshold)
    valley_th = threshold
    if (f[idx] >= f[idx+1] + valley_th) && (f[idx + 1] + valley_th <= f[idx + 2])
        return true 
    else
        return false
    end
end

function hill_start_check(idx, f, threshold)
    if (abs(f[idx] - f[idx+1]) <= threshold) && (f[idx + 1] + threshold <= f[idx + 2])
        return true 
    else
        return false
    end
end

function hill_end_check(idx, f, threshold)
    if (f[idx] > f[idx+1] + threshold) && (abs(f[idx + 1] - f[idx + 2]) <= threshold)
        return true 
    else
        return false
    end
end

function crop_to_true_DRT(DRT, threshold)
    # first artefatic h decay
    division_starting_idx = 1
    while division_starting_idx + 1 <= length(DRT.h) && 
        DRT.h[division_starting_idx + 1] - DRT.h[division_starting_idx] < threshold
        division_starting_idx += 1
    end

    # ending and check for artefacts
    if length(DRT.h) >= 3
        ending_idx = length(DRT.h) - 2
        while ending_idx >= 1 
            if hill_end_check(ending_idx, DRT.h, threshold)
                ending_idx += 1
                break
            end
            ending_idx -=1
        end
    else
        ending_idx = length(DRT.h)
    end

    #@show division_starting_idx, ending_idx

    if ending_idx < division_starting_idx
        println("WARNING: starting idx is greater than ending idx in DRT!")
        return [], []
    end

    ht_cr = DRT.h[division_starting_idx : ending_idx]
    t_cr = DRT.tau_range[division_starting_idx : ending_idx]
    return ht_cr, t_cr
end

function find_divisions_valley(h_tau, threshold)
    # main loop (with an assumption that in the low frequenceis there is no artefact)
    active_idx = 1
    ending_idx = length(h_tau)
    division_idxs = [1]
    while active_idx < ending_idx - 1
        if valley_check(active_idx, h_tau, threshold) || hill_start_check(active_idx, h_tau, threshold)
            #@show valley_check(active_idx, h_tau, threshold), hill_start_check(active_idx, h_tau, threshold)
            push!(division_idxs, active_idx+1)
        end
    active_idx += 1
    end
    push!(division_idxs, ending_idx)
    return division_idxs
end


function obtain_summed_peaks(tau, h_tau, division_idxs, control::DRT_control_struct=DRT_control_struct())
    peak_list = []
    for i in 1:length(division_idxs)-1
        R = discrete_integrate(division_idxs[i], division_idxs[i+1] - 1, h_tau)
        h_max = -1.0
        f_max = -1.0
        for i in division_idxs[i] : division_idxs[i+1]
            if h_tau[i] > h_max 
                h_max = h_tau[i]
                f_max = 1/(2*pi*tau[i])
            end
        end
        push!(peak_list, (f_max, R))
    end
    return peak_list
end

function divide_and_evaluate_R_peaks(DRT::DRT_struct)
    overall_threshold = maximum(DRT.h)*1e-5
    ht_cr, t_cr = crop_to_true_DRT(DRT, overall_threshold)
    threshold = maximum(ht_cr)/10000.0
    
    if length(ht_cr) < 1
        DRT.R_peak_list = []
        return 
    end

    plot_bool = DRT.control.R_peaks_plot_bool
    if plot_bool
        act_fig = PyPlot.gcf()
        figure(11)
    end
    
    #serialize("last_lt_ht.dat", [log.(10, t_cr), ht_cr])
    
    if occursin("valley",DRT.control.divide_R_peaks)
        division_idxs = find_divisions_valley(ht_cr, threshold)
    elseif occursin("curv", DRT.control.divide_R_peaks)
        division_idxs = find_divisions_curv(ht_cr)
        #division_idxs = find_divisions_maximum(ht_cr)
    else
        println("ERROR: divide_R_peaks-division = \"$(DRT.control.divide_R_peaks)\", which is not a valid option!")
        throw(Exception)
    end

    plot_bool && plot_divisions(ht_cr, division_idxs)

    if occursin("gauss", DRT.control.divide_R_peaks)
        DRT.R_peak_list = obtain_gaussian_peaks(t_cr, ht_cr, division_idxs, DRT.control)
    elseif occursin("sum", DRT.control.divide_R_peaks)
        DRT.R_peak_list = obtain_summed_peaks(t_cr, ht_cr, division_idxs, DRT.control)
    else
        println("ERROR: divide_R_peaks-method = \"$(DRT.control.divide_R_peaks)\", which is not a valid option!")
        throw(Exception)
    end

    plot_bool && figure(act_fig)
    return
end
  