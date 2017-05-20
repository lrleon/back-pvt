# ifndef UNDERSATURATED_SOLUTION_GAS_OIL_RELATION_IMPL_H
# define UNDERSATURATED_SOLUTION_GAS_OIL_RELATION_IMPL_H

# include <correlations/undersaturated-oil-formation-volume-factor.H>

inline double RsAbovePb::impl(const double &,
			      const double &,
			      const double & rsb) noexcept
{
  return rsb;
}


# endif
