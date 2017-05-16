# ifndef PURE_WATER_DENSITY_IMPL_H
# define PURE_WATER_DENSITY_IMPL_H

# include <correlations/pure-water-density.H>

inline double
PpwSpiveyMN::impl(const double & t, const double & p) noexcept
{
  const double t100 = t/100;
  const double t100_2 = t100*t100;

  // Pure water density at the reference pressure (70 Mpa) [g/cm³]
  const double ppwr = (-0.127213*t100_2 + 0.645486*t100 + 1.03265) /
    (-0.070291*t100_2 + 0.639589*t100 + 1);
        
  // The coefficients of pure water density, expressed as functions of temperature
  const double epw = (4.221*t100_2 + -3.478*t100 + 6.221) /
    (0.5182*t100_2 + -0.4405*t100 + 1);
     
  const double fpw = (-11.403*t100_2 + 29.932*t100 + 27.952) /
    (0.20684*t100_2 + 0.3768*t100 + 1);
        
  const double ipwr = (1/epw)*log(fabs(epw + fpw));
  const double ipw = (1/epw)*log(fabs((epw * (p/70))+ fpw));
        
  // Pure water density [g/cm³]
  const double ppw = ppwr * exp(ipw - ipwr);
        
  return ppw;
}

# endif
