#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h> /* GLoBES library */
#define MAX_SYS 200
double sys_dcerrors[MAX_SYS]; /* Uncertainties of systematics params */
double sys_dcstartval[MAX_SYS];

/***************************************************************************
 *                      C H I ^ 2   F U N C T I O N S                      *
 ***************************************************************************/

/***************************************************************************
 * Calculate chi^2 for Double Chooz, including the following systematical  *
 * errors   (5 errors in total):                                           *
 *   x[0]: Flux normalization of reactor                                   *
 *   x[1]: Fiducial mass error - far detector                              *
 *   x[2]: Fiducial mass error - near detector                             *
 *   x[3]: Energy calibration error - far detector                         *
 *   x[4]: Energy calibration error - near detector                        *
 *									   *
 * The measured data must be included in user_data in order {FD, ND}.	   *
 ***************************************************************************/
double chiDCNorm(int exp, int rule, int n_params, double *x, double *errors, void *user_data)
{
  /* Fit rates are calculated with GLoBES as usual, but true rates are taken from user input */
  int n_bins_ND = glbGetNumberOfBins(EXP_DC_ND);
  int n_bins_FD = glbGetNumberOfBins(EXP_DC_FD);
  double signal_fit_rates_ND[n_bins_ND];
  double signal_fit_rates_FD[n_bins_FD];
  double *bg_rates_ND;
  double *bg_rates_FD;
  double signal_norm_ND, signal_norm_FD;
  int ew_low, ew_high, i, k;
  double emin, emax, fit_rate, true_rate;
  double chi2_ND = 0.0, chi2_FD = 0.0, chi2_error = 0.0;
  double chi2;

  double freedom = 0;
  /* Initialize bin counter for the user input */
  k = 0;

  /* Apply far detector systematic uncertainties */
  glbGetEminEmax(EXP_DC_FD, &emin, &emax);
  glbGetEnergyWindowBins(EXP_DC_FD, 0, &ew_low, &ew_high);

  glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(EXP_DC_FD, 0),
                      signal_fit_rates_FD, n_bins_FD, emin, emax);

  bg_rates_FD = glbGetBGFitRatePtr(EXP_DC_FD, 0);
  signal_norm_FD = 1.0 + x[0] + x[1];
  // signal_norm_FD = 1.0;

  for (i = ew_low; i <= ew_high; i++)
  {
    true_rate = ((double *)user_data)[k];
    fit_rate = signal_norm_FD * signal_fit_rates_FD[i] + bg_rates_FD[i];
    if (k < 38)
    {
      if (fit_rate > bg_rates_FD[i])
      {
        chi2_FD += likelihood(true_rate, fit_rate, true_rate);
        freedom++;
      }
      if (fit_rate <= bg_rates_FD[i])
      {
        chi2_FD += 0;
      }
    }
    k++;
  }
  /* Apply near detector systematic uncertainties */
  glbGetEminEmax(EXP_DC_ND, &emin, &emax);
  glbGetEnergyWindowBins(EXP_DC_ND, 0, &ew_low, &ew_high);
  glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_DC_ND, 0),
                      signal_fit_rates_ND, n_bins_ND, emin, emax);

  bg_rates_ND = glbGetBGFitRatePtr(EXP_DC_ND, 0);

  signal_norm_ND = 1.0 + x[0] + x[2];
  // signal_norm_ND = 1.0;

  for (i = ew_low; i <= ew_high; i++)
  {
    true_rate = ((double *)user_data)[k];
    fit_rate = signal_norm_ND * signal_fit_rates_ND[i] + bg_rates_ND[i];
    if (k < 76)
    {
      if (fit_rate > bg_rates_ND[i])
      {
        chi2_ND += likelihood(true_rate, fit_rate, true_rate);
        freedom++;
      }
      if (fit_rate <= bg_rates_ND[i])
      {
        chi2_ND += 0;
      }
    }
    k++;
  }
  /* Systematical part of chi^2 (= priors) */
  for (i = 0; i < n_params; i++)
  {
    chi2_error += SQR((x[i]) / errors[i]);
    sys_dcstartval[i] = x[i];
  }

  chi2 = chi2_FD + chi2_ND + chi2_error;
  return chi2;
}

/* Provide measured data for Double Chooz (FD data, ND data),Data from DooubleChooze collbration 2020 */
double DoubleChooz_true_rates[76] = {4718.384683, 5967.473439, 6873.508091, 7798.976045, 8549.674769,
                                     9203.282113, 9578.584028, 9267.793554, 8879.345001, 8452.045661,
                                     7694.649071, 6911.378316, 6361.10157, 6005.015492, 5500.050165,
                                     4839.741001, 3978.800305, 3447.941044, 2878.262627, 2256.788435,
                                     1803.599115, 1395.733076, 1071.993658, 819.4556402, 579.8705204,
                                     476.1998902, 307.8161705, 223.5630255, 158.6759655, 139.0449709,
                                     112.9414812, 80.35363467, 60.62379314, 53.70846391, 40.12689965,
                                     32.78850569, 57.81258726, 37.53711083,
                                     11903.96203, 16353.48067, 20589.1539, 24118.90103, 27263.66516,
                                     28889.69814, 29809.8196, 29746.013, 28762.40572, 27265.44405,
                                     25233.72301, 22731.38458, 20742.45887, 19117.11121, 17128.1855,
                                     15203.41634, 12893.59127, 10626.57612, 8744.558541, 7183.425741,
                                     5772.040142, 4360.610799, 3077.552865, 2265.097748, 1580.941133,
                                     1239.059669, 747.4164233, 491.0964572, 277.7467931, 192.8851853,
                                     193.5850749, 130.1575845, 88.41208913, 68.47981803, 49.18911234,
                                     73.437369, 76.29525132, 57.76275931};
/* Provide measured data for NOSC Double Chooz (FD data, ND data),Data from DooubleChooze collbration 2020 */
double DoubleChooz_NOSC_rates[76] = {
    4951.394507, 6394.654162, 7365.421673, 8400.910182, 9261.641323,
    9876.421604, 10174.05358, 9889.157039, 9487.751633, 8950.427785,
    8212.456588, 7416.221073, 6671.781335, 5959.704072, 5305.871358,
    4626.152617, 3998.213836, 3357.326112, 2787.643742, 2347.395458,
    1965.415446, 1564.010041, 1188.494616, 890.6609942, 670.4814981,
    482.6723853, 346.6353256, 230.0276129, 171.6170019, 145.5135121,
    112.9414812, 86.82217591, 60.62379314, 60.18095902, 40.12689965,
    39.2610008, 57.81654114, 37.53315696,
    11732.85362, 16674.33629, 20845.82381, 24846.18834, 27798.40994,
    29852.23582, 30879.33831, 30644.37956, 29532.4592, 27992.71677,
    25704.31125, 23137.84543, 20785.23962, 18646.55213, 16465.0839,
    14476.12903, 12358.86108, 10369.90621, 8659.011625, 7290.406774,
    5964.568092, 4831.228197, 3612.326805, 2692.890653, 1923.201703,
    1410.197245, 1004.144658, 662.277776, 448.9135308, 321.2274313,
    321.9273208, 279.8464617, 195.3639608, 196.822064, 113.3602353,
    116.2181177, 161.8567486, 121.9338823};
/***************************************************************************
 *     U S E R - D E F I N E D   P R O B A B I L I T Y   E N G I N E       *
 ***************************************************************************/
