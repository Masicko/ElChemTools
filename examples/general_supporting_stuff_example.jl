# if you try to find parameters of any function that would fit the data, use the following
ElChemTools.fit_func_prms([1,2,3,4], [1,3,5,8], (x, p) -> p[1]*x + p[2]*x*x + p[3] + 748, plot_bool=true)

# if the fit is out of the bounds, the algorith gives it a constant penalty value 