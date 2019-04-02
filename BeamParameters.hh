//
//  Beam Parameters Evolution for Collision process
//  @authors
//  v.1   (2017) D. Pellegrini, initial port to C++ of github.com/nkarast/beamCal
//  v.1.1 (2018) N. Karastathis, adding emittance evolution
//  v.1.2 (2019) N. Karastathis, adding burn-off of all 4 IPs
//  @lic: GPL3
#pragma once

#include <limits>
#include "gsl_fsolve.hh"
#include "gsl_fint.hh"
#include <tuple>

template <typename T>
constexpr T sqr(const T & t) {return t*t;}

template <typename T>
constexpr T sqr8(const T & t) {return t*t;}

template <typename T>
class CachedVal {
  bool ok_;
  T value_;
  public:
  CachedVal(): ok_(false) {}
  CachedVal(const T & t): ok_(true), value_(t) {}
  CachedVal(const CachedVal &) = default;

  bool ok() const {return ok_;}
  T value() const {if (!ok_) throw; return value_;}
  T value(const T & t) {value_=t; ok_=true; return value_;}
  void clear() {ok_ = false;}
};


struct BeamParameters {
  const double clight = 299792458.;
  const double pi     = 3.1415926535897932385;
  const double pmass  = 0.93827231;
  const double pp_xsect = 110e-31; //m^2
  const double pp_xsect_inel = 81e-31; //m^2
  const double pp_xsect_elast  = 29.7e-31; //m^2
  const double pp_boff = 110e-31;
  std::string beamProfile = "gauss";


  // IP1 & IP5
  double bx = 0.60; // beta x
  double by = 0.15; // beta y
  double dx = 0.;  // separation x
  double dy = 0.;  // separation y
  double alpha = 500e-6;          // crossing angle
  double Nb    = 2736;           // number of collision at IP1 and IP5

  // IP8
  double bx_IP8       = 3.0; // beta x
  double by_IP8       = 3.0; // beta y
  double dx_IP8       = 0.;  // separation x
  double dy_IP8       = 0.;  // separation y
  double polarity_IP8 = -135.0e-6;
  double alpha_IP8    = 2.0*(250e-6 + polarity_IP8); //crossing angle full
  double alpha0_IP8   = 0.0e-6; // for the CC
  double Nb_IP8       = 2374;

  // IP2
  double bx_IP2       = 10.0; // beta x
  double by_IP2       = 10.0; // beta y
  double dx_IP2       = 0.;  // separation x
  double dy_IP2       = 0.;  // separation y
  double polarity_IP2 = -70.0e-6;
  double alpha_IP2    = 2.0*(170e-6 + polarity_IP2); //crossing angle full
  double alpha0_IP2   = 0.0e-6; // for the CC
  double Nb_IP2       = 2258;

  double Npart = 2.2e11;          // bunch charge at begining of coast
  double Nrj   = 7000.;           // collision energy [GeV]
  double gamma = Nrj/pmass;       // relativisitc factor
  double enx   = 2.5e-6;          // r.m.s. horizontal normalised emittance in collision
  double eny   = 2.5e-6;          // r.m.s. horizontal vertical emittance in collision
  double emitx = enx/gamma;       // r.m.s. horizontal physical emittance in collision
  double emity = eny/gamma;       // r.m.s. vertical physical emittance in collision

  double sigz   = 0.076;          // r.m.s. bunch length [m]
  double circum = 26658.8832;     // ring circumference [m]
  double frev   = clight/circum;  // revolution frequency
  double hrf400 = 35640.;         // harmonic number


  double VRFx0 = 6.8;             // actual CC voltage x
  double VRFy0 = 0.0;             // actual CC voltage y

  double omegaCC0 = 2*pi*hrf400/circum; // wave number at 400 Mhz
  double omegaCCx = omegaCC0;
  double omegaCCy = omegaCC0;
  double VRF0   = 6.8; //11.4;    // reference CC voltage
  double alpha0 = 380e-6; //590e-6;  // crossing angle compensated by VRF0

  private:
  CachedVal<double> cachedRloss;
  CachedVal<double> cachedRloss_IP8;
  CachedVal<double> cachedRloss_IP2;
  // Common stuff to compute integrals in other member functions
  const double epsabs = 1e-8;
  const double epsrel = 1e-8;
  const size_t integration_limit;
  gsl_integration_workspace_pp wsp1;
  gsl_integration_workspace_pp wsp2;
  gsl_integration_workspace_pp wsp3;
  gsl_integration_workspace_pp wsp4;
  gsl_integration_workspace_pp wsp5;
  gsl_integration_workspace_pp wsp6;

  public:
  BeamParameters(const size_t limit = 1000000):
    integration_limit(limit),
    wsp1(limit), wsp2(limit),
    wsp3(limit), wsp4(limit),
    wsp5(limit), wsp6(limit)
  {}

  void clear_cache() {
    cachedRloss.clear();
  }

  void clear_cache_IP8() {
    cachedRloss_IP8.clear();
  }

  void clear_cache_IP2() {
    cachedRloss_IP2.clear();
  }

  double maxHcrabbing() const {
    return VRFx0/VRF0*alpha0;
  }

  double maxVcrabbing() const {
    return VRFx0/VRF0*alpha0;
  }

/////////////////

  void setNrj(const double GeV) {
    Nrj = GeV;
    gamma = Nrj/pmass;
    emitx = enx/gamma;
    emity = eny/gamma;
  }

  void setenx(const double um) {
    enx = um;
    emitx = enx / gamma;
  }

  void seteny(const double um) {
    eny = um;
    emity = eny / gamma;
  }

  void setsigz(const double m){
    sigz = m;
  }

/////////////////

  // longitudinal density function
  double rho(double z) const {
    if (beamProfile == "gauss") {
      return exp(-(z*z)/(2*sigz*sigz)) / (sqrt(2*pi)*sigz);
    }
    std::cerr << "Unrecognized profile!\n";
    return std::numeric_limits<double>::quiet_NaN();
  }

  // kernel function
  double kernel(double z, double t) const {
    const double VRFx = std::min(VRFx0, alpha/alpha0*VRF0); //avoid overcrabbing
    const double VRFy = VRFy0;
    return exp(-sqr((dx*sqrt(bx*emitx) + alpha*z - VRFx/VRF0*alpha0/omegaCCx*cos(omegaCCx*t)*sin(omegaCCx*z))/(2*sqrt(bx*emitx))/sqrt(1 + sqr(z/bx))))*
           exp(-sqr((dy*sqrt(by*emity) +           VRFy/VRF0*alpha0/omegaCCy*sin(omegaCCy*t)*cos(omegaCCy*z))/(2*sqrt(by*emity))/sqrt(1 + sqr(z/by))))/
           sqrt(1 + sqr(z/bx))/sqrt(1 + sqr(z/by));
  }

  // Non-normalized 2D pileup density -- Depends on rho function return value and beam profile
  double density(double z, double t) const {
    return 2*kernel(z,t) * rho(z-t) * rho(z+t);
  }

  // Generalized Loss factor - the double integral of density for z, t in the range of (-inf,inf)
  double Rloss() {
    if (cachedRloss.ok()) {
      return cachedRloss.value();
    }
    double result, abserr, inner_result, inner_abserr;
    auto outer = make_gsl_function( [&](double z) {
      auto inner = make_gsl_function( [&](double t) {return density(z,t);} );
      gsl_integration_qagi(inner, epsabs, epsrel, integration_limit,
                           wsp1, &inner_result, &inner_abserr);
      return inner_result;
    } );
    gsl_integration_qagi(outer, epsabs, epsrel, integration_limit, wsp2, &result, &abserr);
    return cachedRloss.value(result);
  }

  // --------------- Same functions for IP8 -----------------
  // kernel function
  double kernel_IP8(double z, double t) const {
    const double VRFx = 0; //std::min(VRFx0, alpha/alpha08*VRF0); //avoid overcrabbing
    const double VRFy = 0; //VRFy0;

    return exp(-sqr8((dx_IP8*sqrt(bx_IP8*emitx) + alpha_IP8*z )/(2*sqrt(bx_IP8*emitx))/sqrt(1 + sqr8(z/bx_IP8))))*
      exp(-sqr8((dy_IP8*sqrt(by_IP8*emity)   )/(2*sqrt(by_IP8*emity))/sqrt(1 + sqr8(z/by_IP8))))/
      sqrt(1 + sqr8(z/bx_IP8))/sqrt(1 + sqr8(z/by_IP8));
  }


  // Non-normalized 2D pileup density -- Depends on rho function return value and beam profile
  double density_IP8(double z, double t) const {
    //std::cout << kernel8(z,t) << '\t' << rho(z-t) << '\t' << rho(z+t) << std::endl;
    return 2*kernel_IP8(z,t) * rho(z-t) * rho(z+t);
  }

  // Loss factor for IP8
  double Rloss_IP8() {
    if (cachedRloss_IP8.ok()) {
      return cachedRloss_IP8.value();
    }
    double result, abserr, inner_result, inner_abserr;
    auto outer = make_gsl_function( [&](double z) {
      auto inner = make_gsl_function( [&](double t) {return density_IP8(z,t);} );
      gsl_integration_qagi(inner, epsabs, epsrel, integration_limit,
                           wsp3, &inner_result, &inner_abserr);
      return inner_result;
    } );
    gsl_integration_qagi(outer, epsabs, epsrel, integration_limit, wsp4, &result, &abserr);
    return cachedRloss_IP8.value(result);
  }

  // --------------- Same functions for IP2 -----------------
  // kernel function
  double kernel_IP2(double z, double t) const {
    const double VRFx = 0; //std::min(VRFx0, alpha/alpha08*VRF0); //avoid overcrabbing
    const double VRFy = 0; //VRFy0;

    return exp(-sqr((dx_IP2*sqrt(bx_IP2*emitx) + alpha_IP2*z )/(2*sqrt(bx_IP2*emitx))/sqrt(1 + sqr(z/bx_IP2))))*
      exp(-sqr((dy_IP2*sqrt(by_IP2*emity)   )/(2*sqrt(by_IP2*emity))/sqrt(1 + sqr(z/by_IP2))))/
      sqrt(1 + sqr(z/bx_IP2))/sqrt(1 + sqr(z/by_IP2));
  }


  // Non-normalized 2D pileup density -- Depends on rho function return value and beam profile
  double density_IP2(double z, double t) const {
    //std::cout << kernel8(z,t) << '\t' << rho(z-t) << '\t' << rho(z+t) << std::endl;
    return 2*kernel_IP2(z,t) * rho(z-t) * rho(z+t);
  }

  // Loss factor for IP2
  double Rloss_IP2() {
    if (cachedRloss_IP2.ok()) {
      return cachedRloss_IP2.value();
    }
    double result, abserr, inner_result, inner_abserr;
    auto outer = make_gsl_function( [&](double z) {
      auto inner = make_gsl_function( [&](double t) {return density_IP2(z,t);} );
      gsl_integration_qagi(inner, epsabs, epsrel, integration_limit,
                           wsp5, &inner_result, &inner_abserr);
      return inner_result;
    } );
    gsl_integration_qagi(outer, epsabs, epsrel, integration_limit, wsp6, &result, &abserr);
    return cachedRloss_IP2.value(result);
  }


  //=========================================================================================
  // ----------- CALCULATIONS FOR IP1/5 ----------
  // Luminosity [Hz/m^2]
  double lumi() {
    return sqr(Npart)*frev*Nb*Rloss()/(4*pi*sqrt(bx*by*emitx*emity));
  }

  // set beta* for the required target luminosity
  double beta4lumi(const double target_lumi) {
    auto lumi_func = [&](const double beta) {
      by = bx = beta;
      clear_cache();
      return lumi()-target_lumi;
    };
    return gsl_fsolve(lumi_func, /*x_lo*/ 0.1, /*x_hi*/ 1.0);
  }

  // set dx for the required target luminosity
  double dx4lumi(const double target_lumi) {
    dx = 0.;
    clear_cache();
    if (lumi() < target_lumi) return 0.;
    auto lumi_func = [&](const double test) {
      dx = test;
      clear_cache();
      return lumi()-target_lumi;
    };
    return gsl_fsolve(lumi_func, /*x_lo*/ 0.0, /*x_hi*/ 10.0);
  }


  double dy4lumi(const double target_lumi) {
    dy = 0.;
    clear_cache();
    if (lumi() < target_lumi) return 0.;
    auto lumi_func = [&](const double test) {
      dy = test;
      clear_cache();
      return lumi()-target_lumi;
    };
    return gsl_fsolve(lumi_func, /*x_lo*/ 0.0, /*x_hi*/ 10.0);
  }

  // r.m.s. Luminous region [m]
  double siglumz() {
    double result, abserr, inner_result, inner_abserr;
    auto outer = make_gsl_function( [&](double z) {
	auto inner = make_gsl_function( [&](double t) {return density(z,t)*z*z;} );
	gsl_integration_qagi(inner, epsabs, epsrel, integration_limit,
			     wsp1, &inner_result, &inner_abserr);
	return inner_result;
      } );
    gsl_integration_qagi(outer, epsabs, epsrel, integration_limit, wsp2, &result, &abserr);
    return sqrt(result/Rloss());
  }

  // total pileup
  double mutot() {
    return lumi()*pp_xsect_inel/(Nb*frev);
  }

  // normalised line PU density [evt/mm] vs. z [m]
  double muz(const double z) {
    double result, abserr;
    auto ddt = make_gsl_function( [&](double t){
	return density(z,t);
      } );
    gsl_integration_qagi(ddt, epsabs, epsrel, integration_limit, wsp2, &result, &abserr);
    return mutot()/Rloss()*result*1e-3;
  }


  // ----------- CALCULATIONS FOR IP8 ----------
  double lumi_IP8() {
    return sqr8(Npart)*frev*Nb_IP8*Rloss_IP8()/(4*pi*sqrt(bx_IP8*by_IP8*emitx*emity));
  }

  // set dx8 for the required target luminosity
  double dx4lumi_IP8(const double target_lumi8) {
    dx_IP8 = 0.;
    clear_cache_IP8();
    if (lumi_IP8() < target_lumi8) return 0.;
    auto lumi_func = [&](const double test) {
      dx_IP8 = test;
      clear_cache_IP8();
      return lumi_IP8()-target_lumi8;
    };
    return gsl_fsolve(lumi_func, /*x_lo*/ 0.0, /*x_hi*/ 3.0);
  }

  // set dx8 for the required target luminosity
  double dy4lumi_IP8(const double target_lumi8) {
    dy_IP8 = 0.;
    clear_cache_IP8();
    if (lumi_IP8() < target_lumi8) return 0.;
    auto lumi_func = [&](const double test) {
      dy_IP8 = test;
      clear_cache_IP8();
      return lumi_IP8()-target_lumi8;
    };
    return gsl_fsolve(lumi_func, /*x_lo*/ 0.0, /*x_hi*/ 3.0);
  }

  // r.m.s. Luminous region IP8 [m]
  double siglumz_IP8() {
    double result, abserr, inner_result, inner_abserr;
    auto outer = make_gsl_function( [&](double z) {
	auto inner = make_gsl_function( [&](double t) {return density_IP8(z,t)*z*z;} );
	gsl_integration_qagi(inner, epsabs, epsrel, integration_limit,
			     wsp3, &inner_result, &inner_abserr);
	return inner_result;
      } );
    gsl_integration_qagi(outer, epsabs, epsrel, integration_limit, wsp4, &result, &abserr);
    return sqrt(result/Rloss_IP8());
  }

  // total pileup IP8
  double mutot_IP8() {
    return lumi_IP8()*pp_xsect_inel/(Nb_IP8*frev);
  }

  // normalised line PU density IP9 [evt/mm] vs. z [m]
  double muz_IP8(const double z) {
    double result, abserr;
    auto ddt = make_gsl_function( [&](double t){
	return density_IP8(z,t);
      } );
    gsl_integration_qagi(ddt, epsabs, epsrel, integration_limit, wsp4, &result, &abserr);
    return mutot_IP8()/Rloss_IP8()*result*1e-3;
  }


  // ----------- CALCULATIONS FOR IP2 ----------
  double lumi_IP2() {
    return sqr(Npart)*frev*Nb_IP2*Rloss_IP2()/(4*pi*sqrt(bx_IP2*by_IP2*emitx*emity));
  }

  // set dx2 for the required target luminosity
  double dx4lumi_IP2(const double target_lumi2) {
    dx_IP2 = 0.;
    clear_cache_IP2();
    if (lumi_IP2() < target_lumi2) return 0.;
    auto lumi_func = [&](const double test) {
      dx_IP2 = test;
      clear_cache_IP2();
      return lumi_IP2()-target_lumi2;
    };
    return gsl_fsolve(lumi_func, /*x_lo*/ 0.0, /*x_hi*/ 10.0);
  }

  // set dy2 for the required target luminosity
  double dy4lumi_IP2(const double target_lumi2) {
    dy_IP2 = 0.;
    clear_cache_IP2();
    if (lumi_IP2() < target_lumi2) return 0.;
    auto lumi_func = [&](const double test) {
      dy_IP2 = test;
      clear_cache_IP2();
      return lumi_IP2()-target_lumi2;
    };
    return gsl_fsolve(lumi_func, /*x_lo*/ 0.0, /*x_hi*/ 10.0);
  }

  // r.m.s. Luminous region IP2 [m]
  double siglumz_IP2() {
    double result, abserr, inner_result, inner_abserr;
    auto outer = make_gsl_function( [&](double z) {
	auto inner = make_gsl_function( [&](double t) {return density_IP2(z,t)*z*z;} );
	gsl_integration_qagi(inner, epsabs, epsrel, integration_limit,
			     wsp5, &inner_result, &inner_abserr);
	return inner_result;
      } );
    gsl_integration_qagi(outer, epsabs, epsrel, integration_limit, wsp6, &result, &abserr);
    return sqrt(result/Rloss_IP2());
  }

  // total pileup IP2
  double mutot_IP2() {
    return lumi_IP2()*pp_xsect_inel/(Nb_IP2*frev);
  }

  // normalised line PU density IP2 [evt/mm] vs. z [m]
  double muz_IP2(const double z) {
    double result, abserr;
    auto ddt = make_gsl_function( [&](double t){
	return density_IP2(z,t);
      } );
    gsl_integration_qagi(ddt, epsabs, epsrel, integration_limit, wsp6, &result, &abserr);
    return mutot_IP2()/Rloss_IP2()*result*1e-3;
  }




  // dedtElastic
  std::tuple<double, double> dedtElastic(){

    double B = 19.9;
    double P = 7000.0;
    double t = 1.0/B;

    double theta   = sqrt(t/sqr(P));
    double theta_x = theta/sqrt(2.0);

    double den_x_dt = 0.5*2.0*gamma*(bx*100.0)*sqr(theta_x)*(lumi()*1.0e-4)*(pp_xsect_elast*1.0e4)/Nb/Npart;
    double den_x_dt_mmmradpersec = den_x_dt*10.0*1000.0;

    double den_dt = 0.5*0.2*gamma*(bx*100.0)*sqr(theta)*(lumi()*1.0e-4)*(pp_xsect_elast*1.0e4)/Nb/Npart;
    double den_dt_mmmradpersec = den_dt*10.*1000.;

    return std::make_tuple(den_x_dt_mmmradpersec*1.0e-6, den_dt_mmmradpersec*1.0e-6);
  }

  // IBS
  std::tuple<double, double, double, double, double> IBSModel(){

    double loc_Npart, loc_enx, loc_eny, loc_bl;
    loc_Npart = Npart/1.0e11; // e11 ppb
    loc_enx   = enx/1.0e-6; // um
    loc_eny   = eny/1.0e-6; // um
    loc_bl    = (sigz*4.0*1.0e9)/clight; // ns

    double V0 = 16.0; // MV
    double t  = 5.0; //minutes

    double a00 = 5.1987e-06*pow(t,1.9415);
    double a01 = 0.018282*pow(t,0.83002)-4.2914;
    double a10 = 0.0032374*pow(t,0.93508);
    double a11 = 0.0062827*pow(t,0.89334)-1.6245;
    double b00 = 1.547e-07*pow(t,1.9416)+1.0489e-05;
    double b01 = 0.064587*pow(t,0.70466)-4.6279;
    double b10 = 1.3308e-07*pow(t,2)-5.3766e-05*t+4.998e-05;
    double b11 = -0.00025467*pow(t,1.7402)+2.6232;
    double b20 = -0.0006784*pow(t,0.98408)+1.0001;
    double b21 = 1.6052e-07*pow(t,2)-5.1973e-06*t+5.2314e-05;

    double a00l=2.5149e-06*pow(t,2)+1.8552e-05*t-2.2517e-05;
    double a01l=-3.8476*pow(t,0.09162)+4.3565;
    double a02l=-9.5113e-08*pow(t,2)+1.6066e-05*t+2.0409e-05;
    double a10l=0.0011588*pow(t,0.92637);
    double a11l=0.0060754*pow(t,1.1078)-1.3471;
    double a12l=-5.7847e-07*pow(t,2.2954);
    double b00l=7.7172e-07*pow(t,2)+3.2533e-05*t-3.1935e-05;
    double b01l=0.00014293*pow(t,2)-0.01935*t-0.70398;
    double b02l=-8.7952e-09*pow(t,2)+1.8045e-06*t-1.6533e-05;
    double b10l=-2.7508e-07*pow(t,2)+0.00075061*t+0.55646;
    double b11l=2.6594e-07*pow(t,2)+7.1229e-06*t-6.3951e-05;
    double b12l=3.0086e-07*pow(t,2)-0.0014477*t+0.44358;

    double a0l = a00l*pow(loc_enx,a01l)+a02l;
    double a1l = a10l*pow(loc_enx,a11l)+a12l;
    double b0l = b00l*pow(loc_enx,b01l)+b02l;
    double b1l = b10l*pow(loc_enx,b11l)+b12l;

    double C0l   = a1l*loc_Npart+a0l;
    double ccSRl = b0l*loc_Npart+b1l;


    double ey=loc_eny;

    double IBSl=C0l*pow(loc_bl,-3.3)+ccSRl;
    double bl=IBSl*loc_bl;

    double a0=a00*pow(loc_bl,a01);
    double a1=a10*pow(loc_bl,a11);
    double b0=b00*pow(loc_bl,b01);
    double b1=b10*pow(loc_bl,b11);
    double b2=b20*pow(loc_bl,b21);

    double C0=a1*loc_Npart+a0;
    double ccSR=b0*pow(loc_Npart,2)+b1*loc_Npart+b2;

    double IBSx=ccSR+C0/pow(loc_enx,2);
    double ex =IBSx*loc_enx;


    double ex_out_norm_m = ex*1.0e-6;
    double ey_out_norm_m = ey*1.0e-6;
    double bl_out_m = (bl*1.0e-9)*clight/4.0;
    //double bl_out_4sigma_s =  bl*1.0e-9;

    return std::make_tuple(ex_out_norm_m, IBSx, bl_out_m, IBSl, ey_out_norm_m);
  }


};
