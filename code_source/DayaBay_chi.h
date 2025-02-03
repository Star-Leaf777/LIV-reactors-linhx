#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>    /* GLoBES library */


//读取ArM参数的值并将ArM过程重置(trigger_DYB = 0)
/***************************************************************************
 *                      C H I ^ 2   F U N C T I O N S                      *
 ***************************************************************************/


/***************************************************************************
 * Calculate chi^2 for DayaBay, including the following systematical       *
 * errors   (10 errors in total):                                           *
 *   x[0]: Flux normalization of reactor                                   *
 *   x[1]: Fiducial mass error - EH1                                       *
 *   x[2]: Fiducial mass error - EH2                                       *
 *   x[3]: Fiducial mass error - EH3                                       *
 *   x[4]: Energy calibration error - EH1                                  *
 *   x[5]: Energy calibration error - EH2                                  *
 *   x[6]: Energy calibration error - EH3     
 *   x[7]: Scaling error - EH1 detector					   *
 *   x[8]: Scaling error - EH2 detector	                             *
 *   x[9]: Scaling error - EH3 detector					                     *										                                       *
 * The measured data must be included in user_data in order {EH1, EH2,EH3}.*
 ***************************************************************************/

/***************************************************************************
 *                      Created by Hyy                                     *
 ***************************************************************************/

double chiDayaBayNorm(int exp, int rule, int n_params, double *x, double *errors, void *user_data)
{
/* Fit rates are calculated with GLoBES as usual, but true rates are taken from user input */
  int n_bins_EH1 = glbGetNumberOfBins(EXP_DYB_EH1);
  int n_bins_EH2 = glbGetNumberOfBins(EXP_DYB_EH2);
  int n_bins_EH3 = glbGetNumberOfBins(EXP_DYB_EH3);
  
  double *bin_widths_EH1;
  bin_widths_EH1 = glbGetBinSizeListPtr(EXP_DYB_EH1);

  double *bin_widths_EH2;
  bin_widths_EH2 = glbGetBinSizeListPtr(EXP_DYB_EH2);
  
  double *bin_widths_EH3;
  bin_widths_EH3 = glbGetBinSizeListPtr(EXP_DYB_EH3);
  
  double signal_fit_rates_EH1[n_bins_EH1];
  double signal_fit_rates_EH2[n_bins_EH2];
  double signal_fit_rates_EH3[n_bins_EH3];
  double *bg_rates_EH1;
  double *bg_rates_EH2;
  double *bg_rates_EH3;
  double signal_norm_EH1,signal_norm_EH2,signal_norm_EH3;
  int ew_low, ew_high, i, k;
  double emin, emax, fit_rate, true_rate;
  double chi2 = 0.0;
  
  double freedom  = 0;
  /* Initialize bin counter for the user input */
  k = 0;

  /* Apply far detector systematic uncertainties */
  glbGetEminEmax(EXP_DYB_EH1, &emin, &emax);
  glbGetEnergyWindowBins(EXP_DYB_EH1, 0, &ew_low, &ew_high);
  
  glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_DYB_EH1, 0),
                      signal_fit_rates_EH1, n_bins_EH1, emin, emax);
                      
  bg_rates_EH1 = glbGetBGFitRatePtr(EXP_DYB_EH1,0);
  signal_norm_EH1 = 1.0+ x[0] + x[1] + x[7];
  
  for (i=ew_low; i <= ew_high; i++)
  {
    true_rate = ((double*)user_data)[k]*(bin_widths_EH1[i]/1E-3);
    fit_rate  = signal_norm_EH1 * signal_fit_rates_EH1[i]+bg_rates_EH1[i];
    // fit_rate  = fit_rate*(1E-3/bin_widths_EH1[i]);
 
      if(fit_rate>bg_rates_EH1[i]) {chi2 += likelihood(true_rate, fit_rate, true_rate);freedom++;}
      if(fit_rate<=bg_rates_EH1[i]) {chi2 += 0;}
 
    k++;
    
  }

  /*EH2*/
  glbGetEminEmax(EXP_DYB_EH2, &emin, &emax);
  glbGetEnergyWindowBins(EXP_DYB_EH2, 0, &ew_low, &ew_high);
   
  glbShiftEnergyScale(x[5], glbGetSignalFitRatePtr(EXP_DYB_EH2, 0),
                      signal_fit_rates_EH2, n_bins_EH2, emin, emax);


  bg_rates_EH2 = glbGetBGFitRatePtr(EXP_DYB_EH2,0);

  signal_norm_EH2 = 1.0 + x[0] + x[2] + x[8];


  for (i=ew_low; i <= ew_high; i++)
  {
    true_rate = ((double*)user_data)[k]*(bin_widths_EH2[i]/1E-3);
    fit_rate  = signal_norm_EH2 * signal_fit_rates_EH2[i]+bg_rates_EH2[i];
    // fit_rate  = fit_rate*(1E-3/bin_widths_EH2[i]);

      if(fit_rate>bg_rates_EH2[i])  {chi2 += likelihood(true_rate, fit_rate, true_rate);freedom++;}
      if(fit_rate<=bg_rates_EH2[i]) {chi2 += 0;}
   
    k++;
  }

  /*EH3*/
  glbGetEminEmax(EXP_DYB_EH3, &emin, &emax);
  glbGetEnergyWindowBins(EXP_DYB_EH3, 0, &ew_low, &ew_high);
   
  glbShiftEnergyScale(x[6], glbGetSignalFitRatePtr(EXP_DYB_EH3, 0),
                      signal_fit_rates_EH3, n_bins_EH3, emin, emax);


  bg_rates_EH3 = glbGetBGFitRatePtr(EXP_DYB_EH3,0);

  signal_norm_EH3 = 1.0 + x[0] + x[3] + x[9];

  for (i=ew_low; i <= ew_high; i++)
  {
    true_rate = ((double*)user_data)[k]*(bin_widths_EH3[i]/1E-3);
    fit_rate  = signal_norm_EH3 * signal_fit_rates_EH3[i]+bg_rates_EH3[i];
    // fit_rate  = fit_rate*(1E-3/bin_widths_EH3[i]);
      if(k<78) {chi2 += likelihood(true_rate, fit_rate, true_rate);freedom++;}
    k++;
  }



  for (i=0; i < n_params; i++)   chi2 += SQR((x[i]) / errors[i]);

  // /*Obtain the best nuisance parameter*/
  // struct chi2_params ArM_params_t;
  // int n = n_params;
  
  // if(trigger_DYB == 0) {ArM_initial(&ArM_params_DYB,n);trigger_DYB++;}
  // ArM_initial(&ArM_params_t,n);
  // ArM_params_t.chi2 = chi2;
  // for (i=0; i < n_params; i++) {
  //   trigger_DYB++;
  //   ArM_params_t.y[i] = x[i];
  //   }
  // //ArM_load(&ArM_params_t, n);
  // ArM_select(&ArM_params_DYB,&ArM_params_t,n);
  return chi2;

}


/*Data From neutrino 2022 coference DayaBay collaboration*/
// double DayaBay_true_rates[78]={102397.7067,290055.3255,402044.1266,491025.8641,565478.4305,
// 							   626004.0631,660495.3869,660476.1153,647742.4138,610186.9076,
// 							   556282.6714,498140.2918,456345.0364,415761.7496,367307.0901,
// 							   320667.707,269185.8022,215888.6209,170460.1359,136536.2533,
// 							   102006.3864,71108.67837,49294.84726,31717.55347,18983.85202,
// 							   633.55344,
// 							   110598.2656,325460.3781,448988.4663,558058.1322,640842.5377,
// 							   710482.2936,749892.0095,753157.1123,738020.7442,695285.3431,
// 							   632832.9676,569723.4173,521728.887,471106.8116,423112.2813,
// 							   365260.1293,309379.5016,249556.9792,196961.0715,158167.5651,
// 							   117401.3805,83207.52068,57558.08698,39136.42195,24000.63083,
// 							   2110.767704,
// 							   36351.15874,101335.7617,139243.1681,168323.9191,192189.7271,
// 							   213447.1837,221665.5229,223264.0025,217640.6412,206200.5314,
// 							   187939.2564,170279.9627,156231.5015,143186.4013,126730.7185,
// 							   111278.0447,94019.9546,75759.20745,60306.1817,46057.11872,
// 							   35820.09185,24780.65776,17553.00991,11328.01915,7109.046438,
// 							   26.219008};
double DayaBay_true_rates[78]={102397.7067 ,290055.3255 ,402044.1266 ,491025.8641 ,565578.4305,
							   626104.0631 ,660595.3869 ,660476.1153 ,647742.4138 ,610186.9076 ,
							   556382.6714 ,498140.2918 ,456345.0364 ,415761.7496 ,367307.0901 ,
							   320667.707 ,269185.8022 ,215888.6209 ,170460.1359 ,136536.2533 ,
							   102006.3864 ,71108.67837 ,49294.84726 ,31717.55347 ,18983.85202 ,
							   633.55344 ,
							   110598.26560 ,325460.37810 ,448988.46630 ,558058.132200 ,640842.537700+100 ,
							   710582.29360 ,749892.00952 ,753157.11230 ,738020.74420 ,695285.34310 ,
							   632932.96760 ,569723.41730 ,521728.8870 ,471106.81160 ,423112.28130 ,
							   365260.12930 ,309379.50160 ,249556.97920 ,196961.07150 ,158167.56510 ,
							   117401.38050 ,83207.520680 ,57558.086980 ,39136.421950 ,24000.630830 ,
							   2110.7677040 ,
							   36501.15874,101635.7617,139493.1681,168673.9191,192289.7271,
							   213547.1837,221715.5229,223314.0025,217640.6412,206200.5314,
							   187939.2564,170079.9627,156081.5015,143286.4013,126780.7185,
							   111278.0447,94019.9546,75759.20745,60306.1817,46057.11872,
							   35820.09185,24780.65776,17553.00991,11328.01915+50,7109.046438,
							   519.219008};
double DayaBay_NOSC_rates[78]={
                 106327.4,297943.5,411556.3,505239.8,578571.5,
                 640876.4,672227.2,672200.4,657756.7,619144.4,
                 561025.8,503755.2,461325.1,418049.5,370955.6,
                 322165.1,270830.4,219074.3,172403.9,136758.7,
                 101961.1,73101.0,48905.6,31917.2,19169.6,
                 674.5,
                 116777.0,335859.3,466198.9,573888.7,658928.6,
                 730100.8,766141.1,767512.9,749933.2,706468.0,
                 643587.1,577470.6,527995.3,480370.1,425809.5,
                 369400.5,312067.2,252423.3,199712.2,159481.6,
                 119251.1,83643.2,58205.8,39239.0,23969.9,
                 2101.4,
                 38135.6,107203.4,148022.6,181779.7,208192.1,
                 230649.7,242090.4,242655.4,237429.4,223305.1,
                 202966.1,182203.4,166666.7,151129.9,134745.8,
                 116949.2,98163.8,79661.0,62429.4,50000.0,
                 37429.4,26553.7,18361.6,12005.6,7344.6,
                 565.0};// Usual unity

