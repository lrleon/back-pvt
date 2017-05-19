
# include <tclap/CmdLine.h>

# include <correlations/pvt-correlations.H>

# include "bo-corr.H"

using namespace TCLAP;

CmdLine cmd = { "test-compound-call", ' ', "0" };

ValueArg<string> bob_corr_par = { "b", "bob", "bob correlation", false,
				  "BobAlShammasi", "bob-correlation", cmd };
SaturatedOilVolumeFactor * bob_corr = nullptr;

ValueArg<string> boa_corr_par = { "a", "boa", "bob correlation", false,
				  "BobAlShammasi", "bob-correlation", cmd };
UndersaturatedOilVolumeFactor * boa_corr = nullptr;

ValueArg<double> c_bob = { "", "c-bob", "c-bob", false, 0, "c-bob", cmd };

ValueArg<double> m_bob = { "", "m-bob", "m-bob", false, 1, "m-bob", cmd };

ValueArg<double> c_boa = { "", "c-boa", "c-boa", false, 0, "c-boa", cmd };

ValueArg<double> m_boa = { "", "m-boa", "m-boa", false, 1, "m-boa", cmd };

# define Declare_Par(NAME, val)						\
  ValueArg<double> NAME##_par = { "", #NAME, #NAME, false, val, #NAME, cmd }; \
  double NAME = val

Declare_Par(t, 125);
Declare_Par(api, 15);
Declare_Par(pb, 780);
Declare_Par(psep, 100);
Declare_Par(tsep, 100);
Declare_Par(rs, 500);
Declare_Par(rsb, 500);
Declare_Par(yg, 1);
Declare_Par(yo, 0);
Declare_Par(bobp, 0);
Declare_Par(coa, 0);

void init_dft_values()
{
  bob_corr = &BobAlShammasi::get_instance();
  boa_corr = &BoaPetroskyFarshad::get_instance();
  yo = Quantity<Sg_do>(Quantity<Api>(api)).raw();
  bobp = BobAlShammasi::get_instance().call(yg, yo, rsb, t).raw();
  coa = CoaPetroskyFarshad::get_instance().call(yg, api, rsb, t, pb).raw();
}

int main()
{
  init_dft_values();
  BoCorr corr(&BobAlShammasi::get_instance(), &BoaDeGhetto::get_instance(), 1000);

  cout << "Bob pars:" << endl;
  corr.corr1_me()->names.for_each([] (auto &s) { cout << "  " << s << endl; });
  cout << "Boa pars:" << endl;
  corr.corr2_me()->names.for_each([] (auto &s) { cout << "  " << s << endl; });

  corr.set_api(api);
  corr.set_pb(pb);
  corr.set_psep(psep);
  corr.set_tsep(tsep);
  corr.set_rs(rs);
  corr.set_rsb(rsb);
  corr.set_yo(yo);
  corr.set_yg(yg);
  corr.set_bobp(bobp);
  corr.set_t(t);

  for (double p = 300; p <= 3000; p += 40)
    {
      corr.set_p(p);
      cout << "bo = " << corr.compute(p, true) << endl;
    }
}

