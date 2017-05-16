# ifndef METHANE_SOLUBILITY_IMPL_H
# define METHANE_SOLUBILITY_IMPL_H

# include <correlations/methane-solubility.H>

inline double
MwSpiveyMN::impl(const double & t, // comes in Celsius
		  const double & p,
		  const double & pv,
		  const double & nacl) noexcept
{
  const Quantity<Celsius> tc = t; 

  const double t_100 = t/100;
  const double t_100_2 = t_100*t_100;

  const double a = -0.007751*t_100_2 + 0.013624*t_100 + -0.0781;
  const double b = 0.01193*t_100_2 + 0.0851*t_100 + 1.02766;
  const double c = (1.8316*t_100_2 + -7.8119*t_100 + -3.6231) /
    (-0.10733*t_100_2 + 1.09192*t_100 + 1);

  const double log_p_pv = log(p - pv);
  const double mch4pw = exp(a*log_p_pv*log_p_pv + b*log_p_pv + c);

  const double tk = Quantity<Kelvin>(tc).raw();
  const double p2 = p*p;
  const double c1 = 7.015e-2 + 1.074e-4*tk + 2.260e-1*p/tk + -1.227e-3*p2/tk;
  constexpr double c2 = -6.28e-3;
  const double & m = nacl;
  const double m2 = m*m;

  const double mch4w = mch4pw * exp(-2*c1*m - c2*m2);

  return mch4w;
}




# endif
