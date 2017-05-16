# ifndef GAS_PSEUDO_REDUCED_TEMP_IMPL_H
# define GAS_PSEUDO_REDUCED_TEMP_IMPL_H

# include <correlations/gas-pseudo-reduced-temp.H>

inline double Tpr::impl(const double & t, const double & tpc) noexcept
{
  return t/tpc;
}

# endif
