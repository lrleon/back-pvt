
# include <correlations/saturated-oil-formation-volume-factor.H>
# include <correlations/undersaturated-oil-formation-volume-factor.H>
# include "bo-corr.H"

int main()
{
  BoCorr corr(&BobAlShammasi::get_instance(), &BoaDeGhetto::get_instance(), 1000);
 
  // corr_ptr->set_yg(1);
  // corr_ptr->set_api(10);
  // corr_ptr->set_rsb(1000);
  // corr_ptr->set_t(100);
  // cout << "Pb = " << corr_ptr->invoke(false) << endl;
}
