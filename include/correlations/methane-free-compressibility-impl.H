# ifndef METHANE_FREE_COMPRESSIBILITY_IMPL_H
# define METHANE_FREE_COMPRESSIBILITY_IMPL_H

# include <correlations/methane-free-compressibility.H>

inline double
CgfwSpiveyMN::impl(const double & t, // comes in Celsius
		  const double & p,
		  const double & nacl) noexcept
{
  const double & m = nacl; 

  // Coefficients of compressibility of pure water
  const double t_100 = t/100;
  const double t_100_2 = t_100*t_100;

  const double epw = (4.221*t_100_2 + -3.478*t_100 + 6.221) / 
    (0.5182*t_100_2 + -0.4405*t_100 + 1);

  const double fpw = (-11.403*t_100_2 + 29.932*t_100 + 27.952) /
    (0.20684*t_100_2 + 0.3768*t_100 + 1);

  // Coefficients of gas-free brine compressibility
  constexpr double e = 0.1353;

  const double f1 = (-1.409*t_100_2 + -0.361*t_100 + -0.2532)/(9.216*t_100 + 1);
  const double f2 = (5.614*t_100 + 4.6782) / (-0.307*t_100_2 + 2.6069*t_100 + 1);
  const double f3 = -0.1127*t_100_2 + 0.2047*t_100 + -0.0452;
        
  const double ew = epw + e*m;
  const double fw = fpw + f1*pow(m, 1.5) + f2*m + f3*sqrt(m);

  // Compressibility of gas-free brine
  constexpr double inv_pref = 1.0/70.0; //Pressure of reference: 70 [MPa]
  const double p_70 = p/70;
  const double cgfw1 = 1.0/(ew*p_70 + fw);
  const double cgfw = inv_pref * cgfw1;

  return cgfw;
}

# endif
