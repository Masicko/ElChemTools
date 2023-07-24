using ElChemTools
using Test
using DataFrames

@testset "general testset" begin
    @test ElChemTools.test_DRT()

    @test typeof(ElChemTools.simple_run(
            EIS_simulation( 850, [60], 0.0, 
                            physical_model_name="ysz_model_GAS_LoMA_Temperature", 
                            plot_legend=false), 
            pyplot=1,  
            prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                        "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                        "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                        "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                        "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"],
        
            prms_values=(1, 0.0, 0.0, 0.0, 0.0, 
                        0.0, 24.4626,            0.0, 0.33,
                        0.0, 22.63,              0.0, 0.0088,
                        0.0, 21.94,              0.0, 0.17,
                        0.0, 0.85,       0.0, 3.95,      0.0, 6.818*0.15)
            ,use_experiment=false)
        ) == DataFrames.DataFrame
end
