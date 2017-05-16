
declare_correlation_subtype("GasPseudoreducedTemperature", "GasFunction",
                            "T_{pr}");

begin_correlation("Tpr", "PseudoReducedTemperature") 
add_title("CALCULATION OF THE PSEUDOREDUCED TEMPERATURE")
add_note("The temperature is expressed as a ratio of its critical value.")
add_ref("ahmed:1989")
# add_ref("banzer:1996") Secondary reference
add_parameter("t", "Rankine", "Temperature")
add_parameter("tpc", "Rankine", "Gas pseudocritical temperature")
add_author("Standard Equation")
end_correlation()

