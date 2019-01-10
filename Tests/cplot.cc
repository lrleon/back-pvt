
# include <gmock/gmock.h>

# include <iostream>

# include <correlations/cplot.H>

using namespace std;
using namespace Aleph;
using namespace testing;

TEST(Grid, GridClass)
{
  Cplot::Grid grid("t", Fahrenheit::get_instance(),
		   "p", psia::get_instance(),
		   "rs", SCF_STB::get_instance());
  ASSERT_EQ(grid.ncol(), 3);
  ASSERT_EQ(grid.nrow(), 0);
  ASSERT_EQ(grid.col_index("t"), 0);
  ASSERT_EQ(grid.col_index("p"), 1);
  ASSERT_EQ(grid.col_index("rs"), 2);
  ASSERT_EQ(grid.unit("t"), &Fahrenheit::get_instance());
  ASSERT_EQ(grid.unit("p"), &psia::get_instance());
  ASSERT_EQ(grid.unit("rs"), &SCF_STB::get_instance());

  grid.put_col(false, false, 100);
  grid.put_col(false, false, 200, 90);
  grid.put_row(false, false);

  grid.put_row(false, false, 200, 300, 140);
  grid.put_row(false, false, 300, 500, 340);
  grid.put_row(false, false, Quantity<Fahrenheit>(400), Quantity<psia>(700),
	       Quantity<SCF_STB>(500));

  ASSERT_TRUE(eq(grid.col("t"), {100, 200, 300, 400}));
  ASSERT_TRUE(eq(grid.col("p"), {200, 300, 500, 700}));
  ASSERT_TRUE(eq(grid.col("rs"), {90, 140, 340, 500}));
}

struct TestPlot : public Test
{
  BlackoilGrid cplot;

  TestPlot()
    : cplot(Fahrenheit::get_instance(), psia::get_instance(), true)
  {
    set_ttuner_units();
  }
};

struct BlackoilPlot : public TestPlot
{
  BlackoilPlot()
  {
    cplot.set_boa(BoaMillanArcia::get_instance(),
		  -0.260005454479017, 1.23760105495816);
    cplot.set_bob(BobTotalCFP::get_instance(),
		  -0.743320864010538, 1.73939210267341);
    cplot.set_cob(CobMcCainEtAl::get_instance());
    cplot.set_coa(CoaPerezML::get_instance(),
		  8.6296577461858e-06, -0.477971383291644);
    cplot.set_pb(PbSalazar::get_instance(),
		 58.9216674773746, 0.832262895375856);
    cplot.set_rs(RsSalazar::get_instance(),
		 -2.07970430374413, 1.17945155501256);
    cplot.set_uoa(UoaDeGhettoEtAl::get_instance(),
		  32.567065703422, 1.0);
    cplot.set_uob(UobPerezML::get_instance(),
		  -138.470633573862, 1.18174982673783);
    cplot.set_uod(UodPerezML::get_instance(),
		  -4124.34486516297, 1.0);
    cplot.set_zfactor(ZfactorDranchukAK::get_instance());
    cplot.set_api(8.3, Api::get_instance());
    cplot.set_co2(.86, MolePercent::get_instance());
    cplot.set_h2s(0, MolePercent::get_instance());
    cplot.set_n2(.19, MolePercent::get_instance());
    cplot.set_nacl(0, Molality_NaCl::get_instance());
    cplot.set_psep(100, psia::get_instance());
    cplot.set_rsb(79.5, SCF_STB::get_instance());
    cplot.set_tsep(100, Fahrenheit::get_instance());
    cplot.set_yg(.608, Sgg::get_instance());
  }
};

struct SimplePlot : public Test
{
  SimpleGrid cplot;
  
  SimplePlot() : cplot(Fahrenheit::get_instance(), psia::get_instance(), true)
  {
    cplot.set_boa(BoaMillanArcia::get_instance(),
		  -0.260005454479017, 1.23760105495816);
    cplot.set_bob(BobTotalCFP::get_instance(),
		  -0.743320864010538, 1.73939210267341);
    cplot.set_cob(CobMcCainEtAl::get_instance());
    cplot.set_coa(CoaPerezML::get_instance(),
		  8.6296577461858e-06, -0.477971383291644);
    cplot.set_pb(PbSalazar::get_instance(),
		 58.9216674773746, 0.832262895375856);
    cplot.set_rs(RsSalazar::get_instance(),
		 -2.07970430374413, 1.17945155501256);
    cplot.set_uoa(UoaDeGhettoEtAl::get_instance(),
		  32.567065703422, 1.0);
    cplot.set_uob(UobPerezML::get_instance(),
		  -138.470633573862, 1.18174982673783);
    cplot.set_uod(UodPerezML::get_instance(), -4124.34486516297, 1.0);
    cplot.set_zfactor(ZfactorDranchukAK::get_instance());
    cplot.set_api(8.3, Api::get_instance());
    cplot.set_co2(.86, MolePercent::get_instance());
    cplot.set_h2s(0, MolePercent::get_instance());
    cplot.set_n2(.19, MolePercent::get_instance());
    cplot.set_nacl(0, Molality_NaCl::get_instance());
    cplot.set_psep(100, psia::get_instance());
    cplot.set_rsb(79.5, SCF_STB::get_instance());
    cplot.set_tsep(100, Fahrenheit::get_instance());
    cplot.set_yg(.608, Sgg::get_instance());
  }
};

TEST_F(TestPlot, Set_Optional_Correlations)
{
  // This would be the general structure of this test
  ASSERT_THROW(cplot.set_cob(RsAbovePb::get_instance()), InvalidTargetType);
  ASSERT_NO_THROW(cplot.set_cob(CobMcCainEtAl::get_instance()));
  // The remaining tests would be similar but with others correlations
}

TEST_F(TestPlot, Set_Mandatory_Correlations)
{
  // This would be the general structure of this test
  ASSERT_THROW(cplot.set_pb(RsAlMarhoun::get_instance()),
	       InvalidTargetType);
  ASSERT_NO_THROW(cplot.set_pb(PbAlMarhoun::get_instance()));

  ASSERT_THROW(cplot.set_pb_corr(&RsAlMarhoun::get_instance()),
	       InvalidTargetType);
  ASSERT_NO_THROW(cplot.set_pb_corr(&PbAlShammasi::get_instance()));

  ASSERT_NO_THROW(cplot.set_pb(PbAlMarhoun::get_instance(), 2, 2));
  ASSERT_EQ(cplot.c_pb, 2);
  ASSERT_EQ(cplot.m_pb, 2);

  // The remaining tests would be similar but with others correlations
}

TEST_F(TestPlot, Verify_static_method)
{
  ASSERT_NO_THROW(cplot.verify_correlation(&RsAlMarhoun::get_instance(), "rs"));
  ASSERT_THROW(cplot.verify_correlation(&RsAlMarhoun::get_instance(), "pb"),
	       InvalidTargetType);
  ASSERT_NO_THROW(cplot.set_pb(PbAlMarhoun::get_instance()));
  ASSERT_THROW(cplot.verify_correlation(cplot.pb_corr, "rs"), InvalidTargetType);
}

TEST_F(TestPlot, blackoil_correlations_are_set)
{
  ASSERT_THROW(cplot.correlations_are_set(), CorrelationNotFound);
  cplot.set_pb(PbAlMarhoun::get_instance());
  ASSERT_THROW(cplot.correlations_are_set(), CorrelationNotFound);
  cplot.set_rs(RsAlMarhoun::get_instance());
  ASSERT_THROW(cplot.correlations_are_set(), CorrelationNotFound);
  cplot.set_bob(BobAlShammasi::get_instance());
  ASSERT_THROW(cplot.correlations_are_set(), CorrelationNotFound);
  cplot.set_cob(CobMcCainEtAl::get_instance());
  ASSERT_THROW(cplot.correlations_are_set(), CorrelationNotFound);
  cplot.set_coa(CoaDeGhetto::get_instance());
  ASSERT_THROW(cplot.correlations_are_set(), CorrelationNotFound);
  cplot.set_boa(BoaMcCain::get_instance());
  ASSERT_THROW(cplot.correlations_are_set(), CorrelationNotFound);
  cplot.set_uod(UodBeal::get_instance());
  ASSERT_THROW(cplot.correlations_are_set(), CorrelationNotFound);
  cplot.set_uob(UobBeggsRobinson::get_instance());
  ASSERT_THROW(cplot.correlations_are_set(), CorrelationNotFound);
  cplot.set_uoa(UoaAbedini::get_instance());
  ASSERT_THROW(cplot.correlations_are_set(), CorrelationNotFound);
  cplot.set_zfactor(ZfactorBrillBeggs::get_instance());
  ASSERT_NO_THROW(cplot.correlations_are_set());
}

TEST_F(TestPlot, Set_parameters)
{
  ASSERT_THROW(cplot.test_api(), ValueNotFound);
  ASSERT_THROW(cplot.test_rsb(), ValueNotFound);
  ASSERT_THROW(cplot.test_tsep(), ValueNotFound);
  ASSERT_THROW(cplot.test_tsep2(), ValueNotFound);
  ASSERT_THROW(cplot.test_psep(), ValueNotFound);
  ASSERT_THROW(cplot.test_ogr(), ValueNotFound);
  ASSERT_THROW(cplot.test_n2(), ValueNotFound);
  ASSERT_THROW(cplot.test_co2(), ValueNotFound);
  ASSERT_THROW(cplot.test_h2s(), ValueNotFound);
  ASSERT_THROW(cplot.test_nacl(), ValueNotFound);

  ASSERT_THROW(cplot.set_api(12, Fahrenheit::get_instance()), InvalidUnit);
  ASSERT_THROW(cplot.set_rsb(12, Fahrenheit::get_instance()), InvalidUnit);
  ASSERT_THROW(cplot.set_tsep(12, Api::get_instance()), InvalidUnit);
  ASSERT_THROW(cplot.set_tsep2(12, Api::get_instance()), InvalidUnit);
  ASSERT_THROW(cplot.set_psep(12, Fahrenheit::get_instance()), InvalidUnit);
  ASSERT_THROW(cplot.set_ogr(12, Fahrenheit::get_instance()), InvalidUnit);
  ASSERT_THROW(cplot.set_n2(12, Fahrenheit::get_instance()), InvalidUnit);
  ASSERT_THROW(cplot.set_co2(12, Fahrenheit::get_instance()), InvalidUnit);
  ASSERT_THROW(cplot.set_h2s(12, Fahrenheit::get_instance()), InvalidUnit);
  ASSERT_THROW(cplot.set_nacl(12, Fahrenheit::get_instance()), InvalidUnit);

  ASSERT_NO_THROW(cplot.set_api(12, Api::get_instance()));
  ASSERT_NO_THROW(cplot.set_rsb(800, SCF_STB::get_instance()));
  ASSERT_NO_THROW(cplot.set_tsep(120, Fahrenheit::get_instance()));
  ASSERT_NO_THROW(cplot.set_tsep2(120, Fahrenheit::get_instance()));
  ASSERT_NO_THROW(cplot.set_psep(200, psia::get_instance()));
  ASSERT_NO_THROW(cplot.set_ogr(12, STB_MMscf::get_instance()));
  ASSERT_THROW(cplot.set_n2(111, MolePercent::get_instance()), OutOfUnitRange);
  ASSERT_NO_THROW(cplot.set_co2(10, MolePercent::get_instance()));
  ASSERT_NO_THROW(cplot.set_h2s(10, MolePercent::get_instance()));
  ASSERT_THROW(cplot.set_nacl(10, MolePercent::get_instance()), InvalidUnit);
  ASSERT_NO_THROW(cplot.set_nacl(6, Molality_NaCl::get_instance()));  
}

TEST_F(BlackoilPlot, blackoil_ready)
{
  cplot.init();
  ASSERT_NO_THROW(cplot.init());
}

TEST_F(BlackoilPlot, generate)
{
  BlackoilGrid::Grid grid =
    cplot.generate_grid({100, 200, 300}, {100, 300, 700});

  //  cout << grid << endl;
}

TEST_F(SimplePlot, generate)
{
  SimpleGrid::Grid grid =
    cplot.generate_grid({100, 200, 300}, {100, 300, 700});

  cout << grid << endl;
}

