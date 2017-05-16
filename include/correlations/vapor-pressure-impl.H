# ifndef VAPOR_PRESSURE_IMPL_H
# define VAPOR_PRESSURE_IMPL_H

# include <correlations/vapor-pressure.H>

inline double
PvSpiveyMN::impl(const double & t) noexcept // t comes in Kelvin
{
  const double & tk = t;

  constexpr double tc = 647.096; // °K
  constexpr double pc = 22.064; // MPa

  const double v = 1 - (tk/tc);
  const double v2 = v*v;
  const double v3 = v2*v;
  const double v4 = v3*v;
        
  const double lnpv = ((tc/tk) *
		       (-7.85951783*v + 1.84408259*pow(v, 1.5) +
			-11.7866497*v3 + 22.6807411*pow(v, 3.5) +
			-15.9618719*v4 + 1.80122502 * pow(v, 7.5))) + log(pc);

  const double pv = exp(lnpv);

  return pv;
}


# endif
