# include <correlations/pvt-correlations.H>

int main()
{
  BubblePointPressure * corr_ptr = &PbAlMarhoun::get_instance();

  corr_ptr->set_yg(1);
  corr_ptr->set_api(10);
  corr_ptr->set_rsb(1000);
  corr_ptr->set_t(100);
  cout << "Pb = " << corr_ptr->invoke(false) << endl;
}
