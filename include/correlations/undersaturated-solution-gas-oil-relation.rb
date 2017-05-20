 # coding: utf-8

declare_correlation_subtype("UndersaturatedSolutionGasOilRelation",
                            "OilCorrelation", "R_s");

begin_correlation("RsAbovePb", "SCF_STB")
add_title("Faked constant for internal computations above bubble point")
add_parameter("pb", "psia", "Bubble point pressure")
add_parameter("p", "psia", "Pressure")
add_parameter("rsb", "SCF_STB", "Solution GOR at Pb")
add_precondition("p", "pb")
set_hidden()
end_correlation()

################################################################


