//
//  Beam Parameters Evolution for Collision process
//  @authors
//  v.1   (2017) D. Pellegrini, initial port to C++ of github.com/nkarast/beamCal
//  v.1.1 (2018) N. Karastathis, adding emittance evolution
//  v.1.2 (2019) N. Karastathis, adding burn-off of all 4 IPs
//  @lic: GPL3

// g++ -O3 -Wall -Wextra -pedantic -march=native beamParameters.cpp $(gsl-config --libs)

#include "BeamParameters.hh"
#include "gsl_interpolate.hh"
#include <iostream>
#include <cmath>

int main() {
  BeamParameters calc;

  // levelling settings
  const double target_lumi = 7.5e38;
  const double target_lumi_IP8 = 2.0e37;
  const double target_lumi_IP2 = 1.0e35;
  const double min_beta = 0.15;

  bool considerIP8 = true;
  bool considerIP2 = true;

  double fill_time = 15; //h
  double ist_L = target_lumi; //Hz/m^2
  double int_L = 0.;
  // for IP8
  double ist_L_IP8 = target_lumi_IP8;
  double int_L_IP8 = 0.0;
 // for IP2
  double ist_L_IP2 = target_lumi_IP2;
  double int_L_IP2 = 0.0;

  const double extra_losses = 1.0;     // multiplicative factor
  const size_t time_step = 5*60;       // integration step [s]
  const size_t time_const_beta = 5*60;   // beta* is changed every x seconds
  const size_t print_step = 5*60;       // data is output every x seconds


  // for emittance
  const double tauSRxy_s   = 51.79833076*3600;
  const double tauSRl_s    = 25.89916538*3600;
  const double tau_empirical_h1 = 0.05e-6;
  const double tau_empirical_h2 = 0.05e-6;
  const double tau_empirical_v1 = 0.10e-6;
  const double tau_empirical_v2 = 0.10e-6;

  auto xing_I = make_gsl_lin_interp (

    // Nominal : adaptive @ 6sigma v13
    //{ 0.1,     0.8,     0.9,     1.0,     1.1,        1.16,          1.9,      2.2,      100 },
    //{ 2.0*209, 2.0*209, 2.0*216, 2.0*223, 2.0*232,  2.0*234 ,    2.0*181.7 , 2.0*154,  2.0*154 }

    // Nominal : adaptive @ 5sigma v13
    //{ 0.1,     0.8,     0.9,     1.0,      1.11,        1.9,     2.2,       100 },
    //{ 2.0*180, 2.0*180, 2.0*185, 2.0*196, 2.0*206 ,    2.0*160.8 , 2.0*132,  2.0*132 }

    // Ultimate : adaptive @ 6sigma v13
    //{ 0.1,     0.8,      0.9,         1.0,          1.1,          1.2,        1.51,        1.9,      2.2,        100 },
    //{ 2.0*209, 2.0*209,  2.0*216,     2.0*223,      2.*232,       2.0*236.6,  2.0*263.5,   2.0*247 , 2.0*206.5,  2.0*206.5 }

    // Ultimate : adaptive @ 5sigma v13
    //{ 0.1,       0.8,        0.9,         1.0,        1.2,      1.4,         1.9,      2.2,        100 },
    //{ 2.0*179.7, 2.0*179.7,  2.0*185,     2.0*196,    2.0*209,  2.0*223.5,   2.0*207 , 2.0*178,    2.0*178 }

   // constant crossing
   { 0.1, 100 },
   { 2.0*250, 2.0*250}
  );


  bool levelling = true;
  bool levelling_IP8 = true;
  for (size_t t = 0; t <= fill_time*3600; t += time_step) {


    // update the value of beta* and checking for end of levelling
    if (levelling and (t % time_const_beta == 0)) {
      calc.alpha = xing_I(calc.Npart*1e-11)*1e-6;
      calc.beta4lumi(target_lumi);
      if(considerIP8){
	calc.clear_cache_IP8();
	calc.dy4lumi_IP8(target_lumi_IP8);
	ist_L_IP8 = calc.lumi_IP8();
      }//consider IP8

      if(considerIP2){
	calc.clear_cache_IP2();
	calc.dy4lumi_IP2(target_lumi_IP2);
	ist_L_IP2 = calc.lumi_IP2();
      }//consider IP8

      if (calc.bx < min_beta) {
        levelling = false; //end of levelling
        calc.bx = calc.by = min_beta;
        calc.clear_cache();
        ist_L = calc.lumi();
	if(considerIP8){
	  calc.clear_cache_IP8();
	  calc.dy4lumi_IP8(target_lumi_IP8);
	  ist_L_IP8 = calc.lumi_IP8();
	}//consider IP8

	if(considerIP2){
	  calc.clear_cache_IP2();
	  calc.dy4lumi_IP2(target_lumi_IP2);
	  ist_L_IP2 = calc.lumi_IP2();
	}//consider IP8

      }
    } else {
      calc.alpha = xing_I(calc.Npart*1e-11)*1e-6;
      // restrict the xing angle to aperture limit
      //if (calc.alpha > 250.0e-6){
	//calc.alpha = 250.0e-6;
      //}//if (calc.alpha > 250.0e-6)

      calc.clear_cache();
      ist_L = calc.lumi();

      if(considerIP8){
	calc.clear_cache_IP8();
	calc.dy4lumi_IP8(target_lumi_IP8);
	ist_L_IP8 = calc.lumi_IP8();
      }//consider IP8

      if(considerIP2){
	calc.clear_cache_IP2();
	calc.dy4lumi_IP2(target_lumi_IP2);
	ist_L_IP2 = calc.lumi_IP2();
      }//consider IP2

    }//else



    // printing
    if (t%print_step == 0) {

      if (considerIP2==true && considerIP8==true){
	if (t == 0) {
	  std::cout << "# 1. t [h]\n"
		    << "# 2. beta [m]\n"
		    << "# 3. Xing [urad]\n"
		    << "# 4. int_L [fm^1]\n"
		    << "# 5. ist_L [1e34 Hz/cm^2]\n"
		    << "# 6. total pileup\n"
		    << "# 7. r.m.s. luminous region [cm]\n"
		    << "# 8. peak evt/mm\n"
		    << "# 9. Intensity [10^{11} p]\n"
		    << "# 10. Enx [um]\n"
		    << "# 11. Eny [um]\n"
		    << "# 12. Sigz [m]\n"
		    << "# 13. beta IP8 [m]\n"
		    << "# 14. Xing IP8 [urad]\n"
		    << "# 15. dy IP8 [sigma]\n"
		    << "# 16. ist_L IP8 [e34 Hz/cm^2]\n"
		    << "# 17. int_L IP8 [fm^-1]\n"
		    << "# 18. rms luminous region IP8 [cm]\n"
		    << "# 19. total pileup IP8\n"
		    << "# 20. peak pileup IP8 [evt/mm]\n"
		    << "# 21. beta IP2 [m]\n"
		    << "# 22.  Xing IP2 [urad]\n"
		    << "# 23. dy IP2 [sigma]\n"
		    << "# 24. ist_L IP2 [Hz/m^2]\n"
		    << "# 25. int_L IP2 [fm^-1]\n"
		    << "# 26. rms luminous region IP2 [cm]\n"
		    << "# 27. total pileup IP2\n"
		    << "# 28. peak pileup IP2 [evt/mm]\n"

	    ;
	}// t==0

	std::cout << t/3600. << ' ' << calc.bx << ' '
		  << xing_I(calc.Npart*1e-11) << ' '
		  << int_L*1e-40 << ' ' << ist_L << ' '
		  << calc.mutot() << ' '
		  << calc.siglumz()*1e2 << ' '
		  << calc.muz(0) << ' '
		  << calc.Npart*1e-11 << ' '
		  << calc.enx*1.0e6 << ' '
		  << calc.eny*1.0e6 << ' '
		  << calc.sigz << ' '
		  << calc.bx_IP8 << ' '
		  << calc.alpha_IP8*1.0e6 << ' '
		  << calc.dy_IP8 << ' '
		  << ist_L_IP8 << ' '
		  << int_L_IP8*1.0e-40 << ' '
		  << calc.siglumz_IP8()*1.0e2 << ' '
		  << calc.mutot_IP8() << ' '
		  << calc.muz_IP8(0) << ' '
		  << calc.bx_IP2 << ' '
		  << calc.alpha_IP2*1.0e6 << ' '
		  << calc.dy_IP2 << ' '
		  << ist_L_IP2 << ' '
		  << int_L_IP2*1.0e-40 << ' '
		  << calc.siglumz_IP2()*1.0e2 << ' '
		  << calc.mutot_IP2() << ' '
		  << calc.muz_IP2(0) << '\n';

	//print progress
	if (t%(time_step*200) == 0)
	  std::cerr << "progress: " << (int) (t/fill_time/36) << " % ...\n";


      }// if (considerIP2 and considerIP8)
      else if (considerIP8==true and considerIP2==false){

	if (t == 0) {
	  std::cout << "# 1. t [h]\n"
		    << "# 2. beta [m]\n"
		    << "# 3. Xing [urad]\n"
		    << "# 4. int_L [pb^1]\n"
		    << "# 5. ist_L [Hz/m^2]\n"
		    << "# 6. total pileup\n"
		    << "# 7. r.m.s. luminous region [cm]\n"
		    << "# 8. peak evt/mm\n"
		    << "# 9. Intensity [10^{11} p]\n"
		    << "# 10. Enx [um]\n"
		    << "# 11. Eny [um]\n"
		    << "# 12. Sigz [m]\n"
		    << "# 13. beta IP8 [m]\n"
		    << "# 14. Xing IP8 [urad]\n"
		    << "# 15. dy IP8 [sigma]\n"
		    << "# 16. ist_L IP8 [Hz/m^2]\n"
		    << "# 17. int_L IP8 [pm^-1]\n"
		    << "# 18. rms luminous region IP8 [cm]\n"
		    << "# 19. total pileup IP8\n"
		    << "# 20. peak pileup IP8 [evt/mm]\n"
	    ;
	}// t==0

	std::cout << t/3600. << ' ' << calc.bx << ' '
		  << xing_I(calc.Npart*1e-11) << ' '
		  << int_L*1e-40 << ' ' << ist_L << ' '
		  << calc.mutot() << ' '
		  << calc.siglumz()*1e2 << ' '
		  << calc.muz(0) << ' '
		  << calc.Npart*1e-11 << ' '
		  << calc.enx*1.0e6 << ' '
		  << calc.eny*1.0e6 << ' '
		  << calc.sigz << ' '
		  << calc.bx_IP8 << ' '
		  << calc.alpha_IP8*1.0e6 << ' '
		  << calc.dy_IP8 << ' '
		  << ist_L_IP8 << ' '
		  << int_L_IP8*1.0e-40 << ' '
		  << calc.siglumz_IP8()*1.0e2 << ' '
		  << calc.mutot_IP8() << ' '
		  << calc.muz_IP8(0) << '\n';

	//print progress
	if (t%(time_step*200) == 0)
	  std::cerr << "progress: " << (int) (t/fill_time/36) << " % ...\n";

      }else{

	if (t == 0) {
	  std::cout << "# 1. t [h]\n"
		    << "# 2. beta [m]\n"
		    << "# 3. Xing [urad]\n"
		    << "# 4. int_L [pm^1]\n"
		    << "# 5. ist_L [Hz/m^2]\n"
		    << "# 6. total pileup\n"
		    << "# 7. r.m.s. luminous region [cm]\n"
		    << "# 8. peak evt/mm\n"
		    << "# 9. Intensity [10^{11} p]\n"
		    << "# 10. Enx [um]\n"
		    << "# 11. Eny [um]\n"
		    << "# 12. Sigz [m]\n"
	    ;
	}// t==0

	std::cout << t/3600. << ' ' << calc.bx << ' '
		  << xing_I(calc.Npart*1e-11) << ' '
		  << int_L*1e-40 << ' ' << ist_L << ' '
		  << calc.mutot() << ' '
		  << calc.siglumz()*1e2 << ' '
		  << calc.muz(0) << ' '
		  << calc.Npart*1e-11 << ' '
		  << calc.enx*1.0e6 << ' '
		  << calc.eny*1.0e6 << ' '
		  << calc.sigz << '\n';
	//print progress
	if (t%(time_step*200) == 0)
	  std::cerr << "progress: " << (int) (t/fill_time/36) << " % ...\n";
      }//else IP8

    }// if(t%print_step == 0)

    //
    // update emittances:
    //

    auto ibs_result = calc.IBSModel();

    double exn_IBS_m, IBSx1, bl_IBS_m, IBSl1, eyn_IBS_m;
    exn_IBS_m = std::get<0>(ibs_result);
    IBSx1     = std::get<1>(ibs_result);
    bl_IBS_m  = std::get<2>(ibs_result);
    IBSl1     = std::get<3>(ibs_result);
    eyn_IBS_m = std::get<4>(ibs_result);

    auto elastic_result = calc.dedtElastic();
    double dexdt1, deconvdt1;
    dexdt1    = std::get<0>(elastic_result);
    deconvdt1 = std::get<1>(elastic_result);


    //// ONLY IBS + SR
    calc.setenx(exn_IBS_m+dexdt1*time_step);
    calc.seteny(calc.eny*exp(-2.0*time_step/tauSRxy_s)+dexdt1*time_step);
    calc.setsigz(bl_IBS_m);


    // IBS+SR + EXTRA GROWTH
    //calc.setenx(exn_IBS_m+dexdt1*time_step+time_step*tau_empirical_h1/3600.);
    //calc.seteny(calc.eny*exp(-2.0*time_step/tauSRxy_s)+dexdt1*time_step+time_step*tau_empirical_v1/3600.);
    //calc.setsigz(bl_IBS_m);


    // update intesity based on lumi value
    if (considerIP8==true and considerIP2==true){
      int_L += ist_L*time_step;
      int_L_IP8 += ist_L_IP8*time_step;
      int_L_IP2 += ist_L_IP2*time_step;
      calc.Npart -= (ist_L*calc.pp_boff*time_step/calc.Nb*extra_losses*2 + ist_L_IP8*calc.pp_boff*time_step/calc.Nb_IP8*extra_losses*2 + ist_L_IP2*calc.pp_boff*time_step/calc.Nb_IP2*extra_losses*2) ; //4 ips
    }
    else if (considerIP8==true and considerIP2==false){
      int_L += ist_L*time_step;
      int_L_IP8 += ist_L_IP8*time_step;
      calc.Npart -= (ist_L*calc.pp_boff*time_step/calc.Nb*extra_losses*2 + ist_L_IP8*calc.pp_boff*time_step/calc.Nb_IP8*extra_losses*2) ; //3 ips
    }
    else{
      int_L += ist_L*time_step;
      calc.Npart -= ist_L*calc.pp_boff*time_step/calc.Nb*extra_losses*2; //2 ips
    }



  }// for loop



  return 0;

}
