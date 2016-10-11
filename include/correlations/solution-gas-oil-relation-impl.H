
# ifndef SOLUTION_GAS_OIL_RELATION_IMPL_H
# define SOLUTION_GAS_OIL_RELATION_IMPL_H

# include "solution-gas-oil-relation.H"

inline Quantity<SCF_STB>
RsAlShammasi::impl(const Quantity<Sgg> & yg,
		   const Quantity<psia> & pb,
		   const Quantity<psia> & p,
		   const Quantity<Sg_do> & yo,
		   const Quantity<Fahrenheit> & t,
		   const Quantity<SCF_STB> & rsb) const
{
  if (p >= pb)
    return rsb;

  constexpr double c1 = 0.00538088;
  constexpr double c2 = -1.87784;
  constexpr double c3 = 3.1437;
  constexpr double c4 = 1.32657;
  constexpr double c5 = 1 / 0.715082;
        
  const double rs =
    pow(p.get_value() / c1*pow(yg, c2)*pow(yo, c3)*pow(t, c4), c5);

  return Quantity<SCF_STB>(rs);
}


inline Quantity<SCF_STB>
RsAlMarhoun::impl(const Quantity<Sgg> & yg,
		  const Quantity<psia> & pb,
		  const Quantity<psia> & p,
		  const Quantity<Sg_do> & yo,
		  const Quantity<Rankine> & t,
		  const Quantity<SCF_STB> & rsb) const
{
  if (p >= pb)
    return rsb;
    
  constexpr double c1 = 5.527215;
  constexpr double c2 = -1.841408;
  constexpr double c3_inv = 1/0.783716;

  const double raw_yo = yo.get_value();
  const double raw_yg = yg.get_value();
  const double raw_p  = p.get_value();

  const double rs =
    pow(raw_p/(pow(raw_yo, c1)*exp(c2*floor(raw_yo*raw_yg))), c3_inv) /
    ((t.get_value() + 460) * raw_yg); 

  return Quantity<SCF_STB>(rs);
}

inline Quantity<SCF_STB>
RsDeGhetto::impl(const Quantity<Sgg> & yg,
		 const Quantity<psia> & pb,
		 const Quantity<psia> & p,
		 const Quantity<Api> & api,
		 const Quantity<Fahrenheit> & t,
		 const Quantity<Fahrenheit> & tsep,
		 const Quantity<psia> & psep,
		 const Quantity<SCF_STB> & rsb) const
{
  if (p >= pb)
    return rsb;

  const double raw_api = api.get_value();
  const double raw_yg = yg.get_value();
  const double raw_t = t.get_value();
  const double raw_p = p.get_value();

  double rs = 0;
  if (raw_api <= 10)
    rs = raw_yg * pow(raw_p/10.7025*pow(10, 0.0169*raw_api-0.00156*raw_t),
		      1.1128);
  else
    {
      const double raw_tsep = tsep.get_value();
      const double raw_psep = psep.get_value();
      const double ygcorr = raw_yg *
	(1 + 0.5912 * raw_api * raw_tsep * log10(raw_psep/114.7) * 0.0001);
      rs = ((ygcorr * pow(raw_p, 1.2057))/56.434) *
	pow(10, 10.9267*raw_api/(raw_t + 460));
    }

  return Quantity<SCF_STB>(rs);
}


inline Quantity<SCF_STB>
RsDindorukChristman::impl(const Quantity<Sgg> & yg,
			  const Quantity<psia> & pb,
			  const Quantity<psia> & p,
			  const Quantity<Api> & api,
			  const Quantity<Fahrenheit> & t,
			  const Quantity<SCF_STB> & rsb) const
{
  if (p >= pb)
    return rsb;

  constexpr double a1  = 4.86996e-06;
  constexpr double a2  = 5.730982539;
  constexpr double a3  = 0.0099251;
  constexpr double a4  = 1.776179364;
  constexpr double a5  = 44.2500268;
  constexpr double a6  = 2.702889206;
  constexpr double a7  = 0.744335673;
  constexpr double a8  = 3.35975497;
  constexpr double a9  = 28.10133245;
  constexpr double a10 = 1.57905016;
  constexpr double a11 = 0.928131344;

  const double raw_p = p.get_value();
  const double raw_api = api.get_value();
  const double d = a5 + 2*pow(raw_api, a6)/pow(pb, a7);
  const double d2 = d*d;

  const double a = (a1*pow(raw_api, a2) + a3*pow(t, a4)) / d2;

  const double rs = pow((raw_p/a8 + a9) * pow(yg, a10) * pow(10, a), a11);

  return Quantity<SCF_STB>(rs);
}

inline Quantity<SCF_STB>
RsDoklaOsman::impl(const Quantity<Sgg> & yg,
		   const Quantity<psia> & pb,
		   const Quantity<psia> & p,
		   const Quantity<Sg_do> & yo,
		   const Quantity<Rankine> & t,
		   const Quantity<SCF_STB> & rsb) const
{
  if (p >= pb)
    return rsb;

  constexpr double k = 11956.20204068456430404144; // (1/0.836386 * 10 **4)
  constexpr double d = 1/0.724047;

  const double raw_p = p.get_value();

  const double rs = pow(k * raw_p * pow(yg, 1.01049) * pow(yo, -0.107991) *
			pow(t, 0.952584), d);

  return Quantity<SCF_STB>(rs);
}


inline Quantity<SCF_STB>
RsGlaso::impl(const Quantity<Sgg> & yg,
	      const Quantity<psia> & pb,
	      const Quantity<psia> & p,
	      const Quantity<Api> & api,
	      const Quantity<Fahrenheit> & t,
	      const Quantity<SCF_STB> & rsb) const
{
  if (p >= pb)
    return rsb;

  constexpr double a = -0.30218;
  constexpr double b = 1.7447;
  constexpr double b2 = b*b; // b**2
  const double c = 1.7669 - log10(p);

  // TODO: considerar resolver la eq 2do grado y decidir por el
  // valor de p. La solución es log10(p)=(17669*a-2500*b^2)/(10000*a)
  double r = 0;
  if (b2 - 4*a*c < 0)
    r =  -b/(2*a);
  else
    r = (-b + sqrt(b2 - 4*a*c)) / (2*a);

  const double f = pow(10, r);
  const double rs = yg.get_value() *
    pow((f * pow(api, 0.989) ) / ( pow(t, 0.172) ), 1/0.816);

  return Quantity<SCF_STB>(rs);
}


inline Quantity<SCF_STB>
RsHanafy::impl(const Quantity<psia> & pb,
	       const Quantity<psia> & p,
	       const Quantity<SCF_STB> & rsb) const
{
  if (p >= pb)
    return rsb;

  const double rs = -49.069 + 0.312*p.get_value();

  return Quantity<SCF_STB>(rs);
}

inline Quantity<SCF_STB>
RsKartoatmodjoSchmidt::impl(const Quantity<Sgg> & yg,
			    const Quantity<psia> & pb,
			    const Quantity<psia> & p,
			    const Quantity<Api> & api,
			    const Quantity<Fahrenheit> & t,
			    const Quantity<Fahrenheit> & tsep,
			    const Quantity<psia> & psep,
			    const Quantity<SCF_STB> & rsb) const
{
  if (p >= pb)
    return rsb;

  const double raw_t = t.get_value();

  const double c =
    1 + 0.1595 * pow(api, 0.4078) * pow(tsep, -0.2466) * log10(psep/114.7);
  const double ygcorr = c*yg.get_value();

  double rs = 0;
  const double raw_api = api.get_value();
  if (raw_api > 30)
    {
      constexpr double c1 = 0.0315;
      constexpr double c2 = 0.7587;
      constexpr double c3 = 1.0937;
      constexpr double c4 = 11.289;
      rs = c1 * pow(ygcorr, c2) * pow(p, c3) *
	pow(10, c4*raw_api/(raw_t + 460));
    }
  else
    {
      constexpr double c1 = 0.05958;
      constexpr double c2 = 0.7972;
      constexpr double c3 = 1.0014;
      constexpr double c4 = 13.1405;
      rs = c1 * pow(ygcorr, c2) * pow(p, c3) *
	pow(10, c4*raw_api/(raw_t + 460));
    }

  return Quantity<SCF_STB>(rs);
}


inline Quantity<SCF_STB>
RsLasater::impl(const Quantity<Sgg> & yg,
		const Quantity<psia> & pb,
		const Quantity<psia> & p,
		const Quantity<Sg_do> & yo,
		const Quantity<Rankine> & t,
		const Quantity<SCF_STB> & rsb) const
{
  if (p >= pb)
    return rsb;

  const double raw_api = Quantity<Api>(yo).get_value();
  const double mo =
    raw_api <= 40 ? 630 - 10*raw_api : pow(73110*raw_api, -1.562);

  const double raw_p = p.get_value();
  const double raw_yg = yg.get_value();
  const double raw_t = t.get_value();

  // Calculation of Xpb = 
  double sy = 0;
  const double xpb = raw_p*raw_yg/raw_t; // Bubble point pressure factor
  if (xpb < 3.29)
    {
      constexpr double c1 = 0.359;
      const double     c2 = (1.473*raw_p*raw_yg) / raw_t;
      constexpr double c3 = 0.476;
      sy = c1 * log(c2 + c3);
    }
  else
    {        
      const double     c1 = (0.121*raw_p*raw_yg) / raw_t;
      constexpr double c2 = -0.236;
      constexpr double c3 = 0.281;
      sy = pow(c1 + c2, c3);
    }        

  if (sy >= 1)
    sy = 0.99999999;
        
  const double rs = (132755*yo.get_value()*sy) / (mo*(1 - sy));
    
  return Quantity<SCF_STB>(rs >= 0 ? rs : 0);
}

inline Quantity<SCF_STB>
RsPetroskyFarshad::impl(const Quantity<Sgg> & yg,
			const Quantity<psia> & pb,
			const Quantity<psia> & p,
			const Quantity<Api> & api,
			const Quantity<Fahrenheit> & t,
			const Quantity<SCF_STB> & rsb) const
{
  if (p >= pb)
    return rsb;

  const double raw_p = p.get_value();
  const double x = 0.0007916 * pow(api, 1.541) - 0.00004561 * pow(t, 1.3911);
        
  const double rs = pow((raw_p/112.727 + 12.34) * pow(yg, 0.8439) * pow(10, x),
			1.73184);

  return Quantity<SCF_STB>(rs);
}



inline Quantity<SCF_STB> RsStanding::impl(const Quantity<Sgg> & yg,
					  const Quantity<psia> & pb,
					  const Quantity<psia> & p,
					  const Quantity<Api> & api,
					  const Quantity<Fahrenheit> & t,
					  const Quantity<SCF_STB> & rsb) const
{
  const double raw_p = p.get_value();
  if (raw_p >= pb)
    return rsb;

  const double raw_yg = yg.get_value();
  const double raw_api = api.get_value();
  const double raw_t = t.get_value();
   
  const double rs = raw_yg *
    pow((raw_p/18.2 + 1.4)*pow(10, 0.0125*raw_api - 0.00091*raw_t), 1.2048);

  return Quantity<SCF_STB>(rs);
}

inline Quantity<SCF_STB> RsTotalCFP::impl(const Quantity<Sgg> & yg,
					  const Quantity<psia> & pb,
					  const Quantity<psia> & p,
					  const Quantity<Api> & api,
					  const Quantity<Fahrenheit> & t,
					  const Quantity<SCF_STB> & rsb) const
{
  const double raw_p = p.get_value();
  if (raw_p >= pb)
    return rsb;

  double c1 = 0, c2 = 0, c3 = 0, c4 = 0;
  const double raw_api = api.get_value();
  if (raw_api <= 10)
    {
      c1 = 12.2651;
      c2 = 0.030405;
      c4 = 0.9669;
    }
  else if (raw_api <= 35)
    {
      c1 = 15.0057;
      c2 = 0.0152;
      c3 = 0.0004484;
      c4 = 1.095;
    }
  else // API <= 45:
    {
      c1 = 112.925;
      c2 = 0.0248;
      c3 = -0.001469;
      c4 = 1.129;
    }
        
  const double raw_yg = yg.get_value();
  const double raw_t = t.get_value();

  const double rs = raw_yg *
    pow((raw_p/c1) * pow(10, c2*raw_api - c3*raw_t), c4);

  return Quantity<SCF_STB>(rs);
    
}

inline Quantity<SCF_STB>
RsVasquezBeggs::impl(const Quantity<Sgg> & yg,
		     const Quantity<psia> & pb,
		     const Quantity<psia> & p,
		     const Quantity<Api> & api,
		     const Quantity<Fahrenheit> & t,
		     const Quantity<Fahrenheit> & tsep,
		     const Quantity<psia> & psep,
		     const Quantity<SCF_STB> & rsb) const
{
  const double raw_p = p.get_value();
  if (raw_p >= pb)
    return rsb;

  const double raw_api = api.get_value();
  double c1 = 0, c2 = 0, c3 = 0;
  if (raw_api <= 30)
    {
      c1 = 0.0362;
      c2 = 1.0937;
      c3 = 25.724;
    }
  else
    {
      c1 = 0.0178;
      c2 = 1.187;
      c3 = 23.931;
    }

  const double raw_yg = yg.get_value();
  const double ygs = raw_yg *
    (1 + 0.00005912*raw_api*tsep.get_value()*log10(psep/114.7));

  const double rs = c1 * ygs * pow(raw_p, c2) *
    exp((c3*raw_api) / (t.get_value() + 460));

  return Quantity<SCF_STB>(rs > 0 ? rs : 0);
}


inline Quantity<SCF_STB> RsVelarde::impl(const Quantity<Sgg> & yg,
					 const Quantity<psig> & pb,
					 const Quantity<psig> & p,
					 const Quantity<Api> & api,
					 const Quantity<Fahrenheit> & t,
					 const Quantity<SCF_STB> & rsb) const
{
  const double raw_pb = pb.get_value();
  const double raw_p  = p.get_value();
  if (raw_p >= raw_pb)
    return rsb;

  constexpr double a0 = 9.73e-07;
  constexpr double a1 = 1.672608;
  constexpr double a2 = 0.92987;
  constexpr double a3 = 0.247235;
  constexpr double a4 = 1.056052;
        
  constexpr double b0 = 0.022339;
  constexpr double b1 = -1.00475;
  constexpr double b2 = 0.337711;
  constexpr double b3 = 0.132795;
  constexpr double b4 = 0.302065;
        
  constexpr double c0 = 0.725167;
  constexpr double c1 = -1.48548;
  constexpr double c2 = -0.164741;
  constexpr double c3 = -0.09133;
  constexpr double c4 = 0.047094;

  const double raw_yg = yg.get_value();
  const double raw_api = api.get_value();
  const double raw_t = t.get_value();

  const double s1 =
    a0 * pow(raw_yg, a1) * pow(raw_api, a2) * pow(raw_t, a3) * pow(raw_pb, a4);

  const double s2 =
    b0 * pow(raw_yg, b1) * pow(raw_api, b2) * pow(raw_t, b3) * pow(raw_pb, b4);

  const double s3 =
    c0 * pow(raw_yg, c1) * pow(raw_api, c2) * pow(raw_t, c3) * pow(raw_pb, c4);
        
  const double pr = raw_p / raw_pb;
        
  const double rsr = s1*pow(pr, s2) + (1 - s1)*pow(pr, s3);
        
  const double rs = rsr*rsb.get_value();

  return Quantity<SCF_STB>(rs);
}

inline Quantity<SCF_STB> RsCegarra::impl(const Quantity<Sgg> & yg,
					 const Quantity<psia> & pb,
					 const Quantity<psia> & p,
					 const Quantity<Api> & api,
					 const Quantity<Fahrenheit> & t,
					 const Quantity<SCF_STB> & rsb) const
{
  const double raw_pb = pb.get_value();
  const double raw_p  = p.get_value();
  if (raw_p >= raw_pb)
    return rsb;

  double l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0;
  const double raw_api = api.get_value();
  if (raw_api < 29.9)
    {
      l1 = 154.158;
      l2 = 0.4577;
      l3 = 0.0006680;
      l4 = 0.000514;
      l5 = 4.70257;
    }
  else 
    {
      l1 = 809.238;
      l2 = 0.32;
      l3 = 0.00061;
      l4 = 0.011;
      l5 = 1.1142;
    }

  const double raw_yg = yg.get_value();
  const double raw_t = t.get_value();

  const double rs =
    raw_yg * pow( (raw_p/l1 + l5) * pow(10, l4*raw_api - l3*raw_t), 1/l2);

  return Quantity<SCF_STB>(rs); 
}

inline Quantity<SCF_STB> RsPerezML::impl(const Quantity<Sgg> &,
					 const Quantity<psia> & pb,
					 const Quantity<psia> & p,
					 const Quantity<Api> &,
					 const Quantity<Fahrenheit> &,
					 const Quantity<SCF_STB> & rsb) const
{
  const double raw_pb = pb.get_value();
  const double raw_p  = p.get_value();
  if (raw_p >= raw_pb)
    return rsb;

  const double rs = (rsb.get_value()/1.031) * pow((raw_p/raw_pb), 0.83);

  return Quantity<SCF_STB>(rs);
}

inline Quantity<SCF_STB>
RsMillanArcia::impl(const Quantity<Sgg> &,
		    const Quantity<psia> & pb,
		    const Quantity<psia> & p,
		    const Quantity<Api> &,
		    const Quantity<Fahrenheit> &,
		    const Quantity<SCF_STB> & rsb) const
{
  const double raw_pb = pb.get_value();
  const double raw_p  = p.get_value();
  if (raw_p >= raw_pb)
    return rsb;

  const double rs = (rsb.get_value()/1.031) * pow(raw_p/raw_pb, 0.83);

  return Quantity<SCF_STB>(rs);
}

inline Quantity<SCF_STB>
RsManucciRosales::impl(const Quantity<Sgg> & yg,
		       const Quantity<psia> & pb,
		       const Quantity<psia> & p,
		       const Quantity<Api> & api,
		       const Quantity<Fahrenheit> & t,
		       const Quantity<SCF_STB> & rsb) const
{
  const double raw_pb = pb.get_value();
  const double raw_p  = p.get_value();
  if (raw_p >= raw_pb)
    return rsb;

  const double x = 0.000922*t.get_value() - 0.0072*api.get_value();
            
  const double rs = yg.get_value() * pow(raw_p/(84.88 * pow(10, x) ), 1/0.53);

  return Quantity<SCF_STB>(rs);
}

# endif // SOLUTION_GAS_OIL_RELATION_IMPL_H