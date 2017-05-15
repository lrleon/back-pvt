# ifndef COMPOUND_CORR_H
# define COMPOUND_CORR_H

# include "correlation.H"

template <class CorrelationType1, class CorrelationType2, class ParsSet>
class CompoundCorrelation : protected ParsSet
{
  double min_val = 0;
  double max_val = numeric_limits<double>::max();

  const Unit & unit; // unit of pivot parameter

  const CorrelationType1 * corr_ptr1 = nullptr;
  const CorrelationType2 * corr_ptr2 = nullptr;

  double c1 = 0, c2 = 0, m1 = 1, m2 = 1;

  const Quantity<psia> pivot;

public:

  const Unit * result_unit = nullptr; // it will be the result unit of
				      // the first defined correlation

  CompoundCorrelation(const CorrelationType1 * corr_ptr1,
		      const CorrelationType2 * corr_ptr2,
		      const Quantity<psia> & pivot)
    : corr_ptr1(corr_ptr1), corr_ptr2(corr_ptr2), pivot(pivot) {}

  CompoundCorrelation(const CorrelationType1 * corr_ptr1,
		      const CorrelationType2 * corr_ptr2,
		      const Quantity<psia> & pivot,
		      double c1, double m1,
		      double c2, double m2)
    : CompoundCorrelation(corr_ptr1, corr_ptr2, pivot),
      c1(c1), m1(m1), c2(c2), m2(m2) {}

  void set_min(double min_val) { this->min_val = min_val; }

  void set_max(double max_val) { this->max_val = max_val; }

  void set_min(const VtlQuantity & min_val)
  {
    set_min(VtlQuantity(unit, min_val).raw());
  }

  void set_max(const VtlQuantity & max_val)
  {
    set_max(VtlQuantity(unit, max_val).raw());
  }

# define DRY_CALL_CORR()						\
  VtlQuantity val;							\
  if (p_val <= pivot)							\
    val = c1 + m1*corr_ptr1->compute_by_names(pars, check);		\
  else									\
    val = c2 + m2*corr_ptr2->compute_by_names(pars, check);		\
  if (val < min_val)							\
    return VtlQuantity(val.unit, min_val);				\
  if (val > max_val)							\
    return VtlQuantity(val.unit, max_val);				\
  return val

# define DRY_COMPUTE(TYPE)						\
  VtlQuantity								\
  compute_by_names(const TYPE & pars, bool check = true) const		\
  {									\
    auto p_par_ptr = pars.find_ptr([] (auto p) { return p.first == "p"; }); \
    if (p_par_ptr == nullptr)						\
      ZENTHROW(ParameterNameNotFound, "parameter p not found");		\
									\
    VtlQuantity p_val = { p_par_ptr->first, p_par_ptr->second };	\
									\
    DRY_CALL_CORR();							\
  }

  DRY_COMPUTE(DynList<Correlation::ParByName>);
  DRY_COMPUTE(DynList<Correlation::NamedPar>);

  VtlQuantity compute_by_names(const ParList & pars, bool check = true) const
  {
    VtlQuantity p_val = pars.search("p");

    DRY_CALL_CORR();
  }

# undef DRY_CALL_CORR
# undef DRY_COMPUTE
  
  VtlQuantity compute(bool check = true) const
  {
    VtlQuantity p_val = this->p;
    VtlQuantity val;
    if (p_val <= pivot)
      val = c1 + m1*corr_ptr1->compute(check);
    else
      val = c2 + m2*corr_ptr2->compute(check);
    if (val < min_val)
      return VtlQuantity(val.unit, min_val);
    if (val > max_val)
      return VtlQuantity(val.unit, max_val);
    return val;
  }
};


# endif 
