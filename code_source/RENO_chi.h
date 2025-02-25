#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h> /* GLoBES library */
#define SQR(x) ((x) * (x)) /* macro to calculate squares */
#define MAX_SYS 200
double sys_RNerrors[MAX_SYS];       /* Uncertainties of systematics params */
double sys_RNstartval[MAX_SYS]; 
/***************************************************************************
 *                      C H I ^ 2   F U N C T I O N S                      *
 ***************************************************************************/


/***************************************************************************
 * Calculate chi^2 for RENO, including the following systematical	   *
 * errors   (7 errors in total):                                           *
 *   x[0]: Flux normalization of the reactors                              *
 *   x[1]: Detection error - far detector                                  *
 *   x[2]: Detection error - near detector                                 *
 *   x[3]: Energy calibration error - far detector                         *
 *   x[4]: Energy calibration error - near detector                        *
 *   x[5]: Scaling error - far detector					   *
 *   x[6]: Scaling error - near detector				   *
 *									   *
 * The measured data must be included in user_data in order {FD, ND}.	   *
 ***************************************************************************/
double chiRENONorm(int exp, int rule, int n_params, double *x, double *errors, void *user_data)
{
  /* Fit rates are calculated with GLoBES as usual, but true rates are taken from user input */
  int n_bins_FD = glbGetNumberOfBins(EXP_RENO_FD);
  int n_bins_ND = glbGetNumberOfBins(EXP_RENO_ND);
  double *bin_widths_FD;
  double *bin_widths_ND;
  bin_widths_FD = glbGetBinSizeListPtr(EXP_RENO_FD);
  bin_widths_ND = glbGetBinSizeListPtr(EXP_RENO_ND);

  double signal_fit_rates_FD[n_bins_FD];
  double signal_fit_rates_ND[n_bins_ND];
  double *bg_rates_FD, *bg_rates_ND;
  double signal_norm_FD, signal_norm_ND;
  int ew_low, ew_high, i, k;
  double emin, emax, fit_rate, true_rate;
  double chi2 = 0.0;

  /* Initialize bin counter for the user input */
  k = 0;

  /* Apply far detector systematic uncertainties */
  glbGetEminEmax(EXP_RENO_FD, &emin, &emax);
  glbGetEnergyWindowBins(EXP_RENO_FD, 0, &ew_low, &ew_high);

  glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(EXP_RENO_FD, 0),
                      signal_fit_rates_FD, n_bins_FD, emin, emax);

  bg_rates_FD = glbGetBGFitRatePtr(EXP_RENO_FD, 0);
  signal_norm_FD = 1.0 + x[0] + x[1] + x[5];

  for (i = ew_low; i <= ew_high; i++)
  {
    true_rate = ((double *)user_data)[k]* (bin_widths_FD[i]/0.1E-3);
    fit_rate = signal_norm_FD * signal_fit_rates_FD[i] + bg_rates_FD[i];
    // fit_rate = fit_rate * (0.1E-3 / bin_widths_FD[i]);
    if (fit_rate > bg_rates_FD[i])
      chi2 += likelihood(true_rate, fit_rate, true_rate);
    if (fit_rate <= bg_rates_FD[i])
    {
      chi2 += 0;
    }

    k++;
  }
  /* Apply far detector systematic uncertainties */
  glbGetEminEmax(EXP_RENO_ND, &emin, &emax);
  glbGetEnergyWindowBins(EXP_RENO_ND, 0, &ew_low, &ew_high);

  glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_RENO_ND, 0),
                      signal_fit_rates_ND, n_bins_ND, emin, emax);

  bg_rates_ND = glbGetBGFitRatePtr(EXP_RENO_ND, 0);
  signal_norm_ND = 1.0 + x[0] + x[2] + x[6];

  for (i = ew_low; i <= ew_high; i++)
  {
    true_rate = ((double *)user_data)[k]*(bin_widths_ND[i]/0.1E-3);
    fit_rate = signal_norm_ND * signal_fit_rates_ND[i] + bg_rates_ND[i];
    // fit_rate = fit_rate * (0.1E-3 / bin_widths_ND[i]);
    if (fit_rate > bg_rates_ND[i])
      chi2 += likelihood(true_rate, fit_rate, true_rate);
    if (fit_rate <= bg_rates_FD[i])
    {
      chi2 += 0;
    }

    k++;
  }

  /* Systematical part of chi^2 (= priors) */
  for (i = 0; i < n_params; i++){
    chi2 += SQR(x[i] / errors[i]);
    sys_RNstartval[i] = x[i];
  }
  return chi2;
}

/* Provide measured data for RENO (FD data, ND data) */



/* Provide measured data for RENO (FD_data,ND_data) */
double RENO_true_rates[102] = {
    1219.101086 , 1473.353537 , 1663.175406 , 1908.760431 , 2052.338992 ,
    2288.165047 , 2523.998733 , 2694.912146 , 2865.819836 , 3033.314484 ,
    3127.472003 , 3298.630262 , 3416.457115 , 3406.279041 , 3529.346924 ,
    3430.330582 , 3587.576569 , 3573.9721 , 3525.956128 , 3410.354735  ,
    3434.337595 , 3437.823132 , 3325.137456  , 3226.113483 , 3168.097511  ,
    2983.657849 , 2829.969895 , 2761.695722  , 2666.094329 , 2573.909792  ,
    2509.052476 , 2351.947665 , 2239.26199  , 2153.905442  , 2116.388701  ,
    2027.622929 , 1860.298074 , 1638.337855  , 1354.870404 , 1160.252668  ,
    962.208537  , 822.2509663 , 596.9425716  , 303.4401316 , 37.46904607  ,

    9061.050117, 11071.15066, 12943.51016, 14981.15892, 16853.51842,
    18477.94403, 19826.88755, 21864.53631, 23185.93162, 24589.97156,
    25498.14373, 26626.70158, 27507.32554, 28167.56383, 28662.51286,
    29157.46189, 28798.41643, 28935.2109, 28796.57896, 28520.15028,
    28271.3255, 27746.96295, 26644.11583, 25844.29902, 24879.19296,
    23831.44226, 22811.26762, 21983.87476, 21184.05795, 20632.17503,
    19804.83785, 19197.83067, 18232.7246, 17598.22489, 16577.99456,
    15888.37059, 15061.00557, 14150.99593, 13378.72733, 12220.78379,
    11696.47692, 10510.95734, 9545.851273, 8828.679092, 8083.958701,
    7366.75868, 6649.614339, 6042.634995, 5408.107442, 4856.224517,
    4249.217333, 3752.458667, 2952.168569, 2123.884813, 1571.08315,
    879.6217028, 101.8685175};
double RENO_NOSC_rates[45] = {
    1225.680414, 1516.178081, 1765.677285, 2035.675721, 2274.92054,
    2483.425098, 2664.585265, 2944.83427, 3108.91397, 3293.494809,
    3409.732795, 3549.894499, 3676.38687, 3755.04325, 3830.280865,
    3871.348008, 3844.076113, 3837.32062, 3816.890069, 3803.2856,
    3745.271536, 3707.754795, 3536.976836, 3424.289252, 3308.188627,
    3154.487319, 3014.470607, 2929.123598, 2819.849055, 2734.496324,
    2621.814463, 2556.953332, 2423.766517, 2338.415693, 2222.311253,
    2102.78614, 1949.143973, 1706.673077, 1481.296002, 1235.423511,
    1044.213093, 859.8382955, 620.8624749, 306.8569881, 37.46713827};
