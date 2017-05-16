# ifndef PSEUDO_REDUCED_PRESSURE_IMPL_H
# define PSEUDO_REDUCED_PRESSURE_IMPL_H

# include <correlations/gas-pseudo-reduced-pressure.H>

inline double Ppr::impl(const double & p, const double & ppc) noexcept
{
  return p/ppc;
}

# endif
