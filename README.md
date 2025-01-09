# ElChemTools

A toolkit for processing electrochemistry data or running FVM simulations using custom implemented models. In particular, it offers
- preprocessing of EIS (Electrochemical Impedance Spectroscopy) data
- viewing of EIS response to EEC (Electrical Equivalent Circuit) consisting of serial combination of
  - R resistance
  - C capacitance
  - L inductance
  - RC
  - RCPE resistor with constant phase element in paralel
- fitting the EIS data to a predefined EEC structure using Nelder-Mead simplex optimisation method
- DRT (Distribution of Relaxation Times) analysis of EIS data
- implementation of a custom model using definition of 
  - storage term
  - flux term
  - reaction term
  - source term
- numerical simulation (using the defined model) of 
  - EIS (using linearization and Fourier transformation), 
  - CV (Cyclic Voltammetry)
  - differential capacitance measurement of a blocking electrode regime
  - ... more custom made experiments
- running large numerical parametric studies using personal notebook or a computional cluster
  
## Instalation

If you do not have a python installed on you system and you want julia to take care of it (and install its own python), write
```julialang
ENV["PYTHON"] = ""
```
prior to the adding the package. Then the package can be then installed via
```julialang
]
add https://github.com/Masicko/ElChemTools.jl
```

## Usage

A user guide for every feature is documented in `./examples/` folder, e.g. `./examples/EEC_usage.jl`.

## Acknowledgement

This work was supported by the German Research Foundation, DFG project no. FU 316/14-1, and by the Czech Science Foundation, GAČR project no. 19-14244J.





