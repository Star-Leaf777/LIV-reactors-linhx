/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* ##################################################################### *
 *									 *
 * 	This program is to check the antineutrino flux to Daya Bay, Double Chooz and RENO.	 *
 * 			 *
 *									 *
 * ##################################################################### *
 */

#include "myio.h" /* my input-output routines */

#include <globes/globes.h> /* GLoBES library */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Lorentz_violation_lat.h"

#define SQR(x) ((x) * (x)) /* macro to calculate squares */

/* If filename given, write to file; for empty filename write to screen */
char MYFILE1[] = "../Data/DoubleChooz/Double_Chooz_ND_mc.dat";
char MYFILE2[] = "../Data/DoubleChooz/Double_Chooz_FD_mc.dat";
char MYFILE3[] = "../Data/RENO/RENO_ND_mc.dat";
char MYFILE4[] = "../Data/RENO/RENO_FD_mc.dat";
char MYFILE5[] = "../Data/DayaBay/DayaBay_EH1_mc.dat";
char MYFILE6[] = "../Data/DayaBay/DayaBay_EH2_mc.dat";
char MYFILE7[] = "../Data/DayaBay/DayaBay_EH3_mc.dat";
char th13file[] = "../Data/th13_DC_RN_DYB.dat";
char TwoDfile[] = "../Data/TwoD_RN_DYB.dat";

/***************************************************************************
 *                      C H I ^ 2   F U N C T I O N S                      *
 ***************************************************************************/

/* Pre-defined indexing for the reactor experiments */
#define EXP_DC_ND 0   /* Double Chooz near detector */
#define EXP_DC_FD 1   /* Double Chooz far detector */
#define EXP_RENO_ND 2 /* RENO near detector */
#define EXP_RENO_FD 3 /* RENO far detector */
#define EXP_DYB_EH1 4 /* Daya Bay EH1 detectors */
#define EXP_DYB_EH2 5 /* Daya Bay EH2 detectors */
#define EXP_DYB_EH3 6 /* Daya Bay EH3 detectors */

/* Gaussian likelihood function */
inline double likelihood(double true_rate, double fit_rate, double sqr_sigma)
{
  if (sqr_sigma > 0)
    return SQR(true_rate - fit_rate) / sqr_sigma;
  else
    return 0.0;
}
#include "DayaBay_chi.h"
#include "DoubleChooz_chi.h"
#include "RENO_chi.h"

/***************************************************************************
 *                      M A I N   P R O G R A M                            *
 ***************************************************************************/

int main(int argc, char *argv[])
{
  /* Initialize libglobes */
  glbInit(argv[0]);

  glbRegisterProbabilityEngine(19, /*Number of parameters*/
                               &my_probability_matrix,
                               &my_set_oscillation_parameters,
                               &my_get_oscillation_parameters,
                               NULL);
  /* Define chi^2 functions for the reactor experiments */
  glbDefineChiFunction(&chiDCNorm, 5, "chiDCNorm", &DoubleChooz_true_rates);
  glbDefineChiFunction(&chiRENONorm, 7, "chiRENONorm", &RENO_true_rates);
  glbDefineChiFunction(&chiDayaBayNorm, 10, "chiDayaBayNorm", &DayaBay_true_rates);

  /* Initialize experiment(s) used in this simulation */
  glbClearExperimentList();
  glbInitExperiment("../GLB/DoubleChooz/DoubleChooz_JPfix_ND.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("../GLB/DoubleChooz/DoubleChooz_JPfix_FD.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("../GLB/RENO/RENO_JPfix_ND.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("../GLB/RENO/RENO_JPfix_FD.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("../GLB/DayaBay/DayaBay_EH1_2022.glb", &glb_experiment_list[0], &glb_num_of_exps); // 1
  glbInitExperiment("../GLB/DayaBay/DayaBay_EH2_2022.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("../GLB/DayaBay/DayaBay_EH3_2022.glb", &glb_experiment_list[0], &glb_num_of_exps);

  /* Define central values for the prior function (adopted from neutrino 2022 conference results) */
  double theta12 = asin(sqrt(0.304));
  double theta13 = asin(sqrt(0.0892)) / 2;
  double theta23 = asin(sqrt(0.573));
  double deltacp = 197 * M_PI / 180;
  double sdm = 7.42e-5;
  double ldm = 2.519e-3 + SQR(sin(theta12)) * sdm; // This parameter is determined by specific collaboration. RENO(neutrino 2022):mee= 2.74

  double theta12_error = 0.02;
  double theta13_error = 1.0 / (4.0 * sqrt(1 - 0.0852) * sqrt(0.0852)) * (0.0024); // DayaBay
  double theta23_error = 0.03;
  double deltacp_error = 0.9076;
  double sdm_error = 0.3e-5;
  double ldm_error = 0.060e-3;

  /* Initialize the parameter vector */
  glb_params central_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params minimum = glbAllocParams();
  glb_projection SM_projection = glbAllocProjection();

  /* Set the parameter vector */
  glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, ldm);
  glbDefineParams(input_errors, theta12_error, theta13_error, theta23_error, deltacp_error, sdm_error, ldm_error);
  glbSetDensityParams(central_values, 1.0, GLB_ALL);
  glbSetDensityParams(input_errors, 0.05, GLB_ALL);
  /*Lorentz parameters*/
  glbSetOscParams(central_values, 0, GLB_AS_210);
  glbSetOscParams(central_values, 0, GLB_AS_211);
  glbSetOscParams(central_values, 0, GLB_AC_210);
  glbSetOscParams(central_values, 0, GLB_AC_211);
  glbSetOscParams(central_values, 0, GLB_BS_21);
  glbSetOscParams(central_values, 0, GLB_BC_21);

  glbSetOscParams(central_values, 0, GLB_AS_310);
  glbSetOscParams(central_values, 0, GLB_AS_311);
  glbSetOscParams(central_values, 0, GLB_AC_310);
  glbSetOscParams(central_values, 0, GLB_AC_311);
  glbSetOscParams(central_values, 0, GLB_BS_31);
  glbSetOscParams(central_values, 0, GLB_BC_31);

  glbSetOscParams(central_values, 0, GLB_Lat);

  glbSetOscParams(input_errors, 0, GLB_AS_210);
  glbSetOscParams(input_errors, 0, GLB_AS_211);
  glbSetOscParams(input_errors, 0, GLB_AC_210);
  glbSetOscParams(input_errors, 0, GLB_AC_211);
  glbSetOscParams(input_errors, 0, GLB_BS_21);
  glbSetOscParams(input_errors, 0, GLB_BC_21);

  glbSetOscParams(input_errors, 0, GLB_AS_310);
  glbSetOscParams(input_errors, 0, GLB_AS_311);
  glbSetOscParams(input_errors, 0, GLB_AC_310);
  glbSetOscParams(input_errors, 0, GLB_AC_311);
  glbSetOscParams(input_errors, 0, GLB_BS_31);
  glbSetOscParams(input_errors, 0, GLB_BC_31);

  glbSetOscParams(input_errors, 0, GLB_Lat);

  /* The oscillation probabilities are computed */
  glbSetOscillationParameters(central_values);
  glbSetRates();
  /* Set central values and Input_error*/
  glbSetCentralValues(central_values);
  glbSetInputErrors(input_errors);

  /*Set up the Projection  */
  glbDefineProjection(SM_projection, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED);
  glbSetDensityProjectionFlag(SM_projection, GLB_FIXED, GLB_ALL);
  /*Lorentz parameters*/
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_AS_210);
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_AS_211);
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_AC_210);
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_AC_211);
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_BS_21);
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_BC_21);

  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_AS_310);
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_AS_311);
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_AC_310);
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_AC_311);
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_BS_31);
  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_BC_31);

  glbSetProjectionFlag(SM_projection, GLB_FIXED, GLB_Lat);

  glbSetProjection(SM_projection);

  /* Initiate a parameter vector for the scan */
  glbCopyParams(central_values, test_values);

  FILE *fp = fopen(th13file, "w");

  if (!fp)
  {
    perror("文件打开失败");
    return 1;
  }

  double DCth13_min, DCres, RNth13_min, RNres, DYBth13_min, DYBres;
  double DCchi_min = 1e7, RNchi_min = 1e7, DYBchi_min = 1e7;
  double xmin = 0.0700, xmax = 0.1000;
  int xsteps = 30;
  double ymin = 2.300e-3, ymax = 2.800e-3;
  int ysteps = 50;

  double th13;
  double x, y;

  // Th13 scan
  for (x = xmin; x <= xmax; x = x + (xmax - xmin) / xsteps)
  {
    /* Set vector of test values */
    th13 = asin(sqrt(x)) / 2; // Sin2 2theta13 to radian
    glbSetOscParams(test_values, th13, GLB_THETA_13);
    glbChiNP(test_values, minimum, GLB_ALL);
    DCres = glbChiNP(test_values, minimum, EXP_DC_FD);
    RNres = glbChiNP(test_values, minimum, EXP_RENO_FD);
    DYBres = glbChiNP(test_values, minimum, EXP_DYB_EH1);

    if (DCres < DCchi_min)
    {
      DCth13_min = th13;
      DCchi_min = DCres;
    }
    if (RNres < RNchi_min)
    {
      RNth13_min = th13;
      RNchi_min = RNres;
    }
    if (DYBres < DYBchi_min)
    {
      DYBth13_min = th13;
      DYBchi_min = DYBres;
    }
    fprintf(fp, "%f\t%f\t%f\t%f\n", th13, DCres, RNres, DYBres);
    /* Pass fit values to next iteration as starting values */
    glbCopyParams(minimum, test_values);
  }
  printf("Best-fit value: \nDoubleChooz : Sin^2 2theta13 = %g degrees and chi^2_min = %g \nRENO: Sin^2 2theta13 = %g degrees and chi^2_min = %g \nDayaBay Sin^2 2theta13 = %g degrees and chi^2_min = %g \n", SQR(sin(2 * DCth13_min)), DCchi_min, SQR(sin(2 * RNth13_min)), RNchi_min, SQR(sin(2 * DYBth13_min)), DYBchi_min);
  fclose(fp);

  // Two D sin^22th13-Dmee scan
  RNchi_min = 1e7, DYBchi_min = 1e7;
  double RNdm31_min, DYBdm31_min;
  double thedm31;
  int iter = 0;
  FILE *file = fopen(TwoDfile, "w");
  for (y = ymin; y <= ymax; y = y + (ymax - ymin) / ysteps)
  {
    thedm31 = y;
    glbSetOscParams(test_values, thedm31, GLB_DM_31);
    for (x = xmin; x <= xmax; x = x + (xmax - xmin) / xsteps)
    {
      /* Set vector of test values */
      th13 = asin(sqrt(x)) / 2; // Sin2 2theta13 to radian
      glbSetOscParams(test_values, th13, GLB_THETA_13);
      glbChiNP(test_values, minimum, GLB_ALL);
      glbCopyParams(minimum, test_values);
      /* Compute Chi^2 for all loaded experiments and all rules */
      RNres = glbChiNP(test_values, minimum, EXP_RENO_FD);
      glbSetSysStartingValuesList(EXP_RENO_FD, 0, GLB_ON, sys_dybstartval);
      DYBres = glbChiNP(test_values, minimum, EXP_DYB_EH1);
      glbSetSysStartingValuesList(EXP_DYB_EH1, 0, GLB_ON, sys_dybstartval);
      
      /* Pass fit values to next iteration as starting values */
      if (RNres < RNchi_min)
      {
        RNth13_min = th13;
        RNdm31_min = thedm31;
        RNchi_min = RNres;
      }
      if (DYBres < DYBchi_min)
      {
        DYBth13_min = th13;
        DYBdm31_min = thedm31;
        DYBchi_min = DYBres;
      }
      /* Print the data point to output */
      fprintf(file, "%f\t%f\t%f\t%f\n", x, y - SQR(sin(theta12)) * sdm, RNres, DYBres); // th13,dmee,RNres,DYBres
      // glbCopyParams(minimum, test_values);
      iter++;
      if ((iter + 1) % 100000 == 0)
        printf("[/]");
        fflush(stdout);
    }
  }
  fclose(file);

  double RENOdmee = SQR(cos(theta12)) * fabs(RNdm31_min) + SQR(sin(theta12)) * fabs(RNdm31_min - sdm);  // Normal order
  double DYBdmee = SQR(cos(theta12)) * fabs(DYBdm31_min) + SQR(sin(theta12)) * fabs(DYBdm31_min - sdm); // Normal order
  /* When the lowest chi^2 and Dmee for theta13 has been found, print the fit values */
  printf("RENO Best-fit value: Sin2 2theta13 = %g , DMee = %g and chi^2_min = %g \n\n", SQR(sin(2 * RNth13_min)), RENOdmee, RNchi_min);
  printf("DYB Best-fit value: Sin2 2theta13 = %g , DMee = %g and chi^2_min = %g \n\n", SQR(sin(2 * DYBth13_min)), DYBdmee, DYBchi_min);
  printf("RENO Best-fit value: Sin2 2theta13 = %g , DMee = %g and chi^2_min = %g \n\n", SQR(sin(2 * RNth13_min)), RNdm31_min - SQR(sin(theta12)) * sdm, RNchi_min);
  printf("DYB Best-fit value: Sin2 2theta13 = %g , DMee = %g and chi^2_min = %g \n\n", SQR(sin(2 * DYBth13_min)), DYBdm31_min - SQR(sin(theta12)) * sdm, DYBchi_min);

  // int i, n_bins;

  // double *true_rates;
  // double *bin_centers;
  // double *bin_widths;

  // double total_rate, bg_rate;

  // n_bins = glbGetNumberOfBins(0);
  // bin_centers = glbGetBinCentersListPtr(0);
  // bin_widths = glbGetBinSizeListPtr(0);

  // /*******************************
  //  * Produce DoubleChooz spectra *
  //  * *****************************/
  // glbDefineParams(central_values, theta12, DCth13_min, theta23, deltacp, sdm, ldm);
  // glbSetOscillationParameters(central_values);
  // glbSetRates();
  // /* Compute the combined event spectrum for DC_ND */
  // InitOutput(MYFILE1, "");
  // true_rates = glbGetRuleRatePtr(0, 0);
  // for (i = 0; i <= n_bins - 1; i++)
  //   AddToOutput(bin_centers[i] - 0.80e-3, true_rates[i], bin_widths[i]);

  // /* Provide the total IBD events for DC_ND */
  // total_rate = 0.0;
  // total_rate = glbTotalRuleRate(0, 0, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG);
  // printf("Total IBD rate in DoubleChooz_ND: %g\n", total_rate);

  // n_bins = glbGetNumberOfBins(1);
  // bin_centers = glbGetBinCentersListPtr(1);
  // bin_widths = glbGetBinSizeListPtr(1);

  // /* Compute the combined event spectrum for DC_FD */
  // InitOutput(MYFILE2, "");
  // true_rates = glbGetRuleRatePtr(1, 0);
  // for (i = 0; i <= n_bins - 1; i++)
  //   AddToOutput(bin_centers[i] - 0.80e-3, true_rates[i], bin_widths[i]);

  // /* Provide the total IBD events for DC_FD */
  // total_rate = 0.0;
  // total_rate = glbTotalRuleRate(1, 0, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG);
  // printf("Total IBD rate in DoubleChooz_FD: %g\n", total_rate);

  // /************************
  //  * Produce RENO spectra *
  //  * **********************/
  // glbDefineParams(central_values, theta12, RNth13_min, theta23, deltacp, sdm, RNdm31_min);
  // glbSetOscillationParameters(central_values);
  // glbSetRates();

  // n_bins = glbGetNumberOfBins(2);
  // bin_centers = glbGetBinCentersListPtr(2);
  // bin_widths = glbGetBinSizeListPtr(2);

  // /* Compute the combined event spectrum for RENO_ND */
  // InitOutput(MYFILE3, "");
  // true_rates = glbGetRuleRatePtr(2, 0);
  // for (i = 0; i <= n_bins - 1; i++)
  //   AddToOutput(bin_centers[i] - 0.78e-3, true_rates[i] * (0.1E-3 / bin_widths[i]), bin_widths[i]);

  // /* Provide the total IBD events for RENO_ND */
  // total_rate = 0.0;
  // total_rate = glbTotalRuleRate(2, 0, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG);
  // printf("Total IBD rate in RENO_ND: %g\n", total_rate);

  // n_bins = glbGetNumberOfBins(3);
  // bin_centers = glbGetBinCentersListPtr(3);
  // bin_widths = glbGetBinSizeListPtr(3);

  // /* Compute the combined event spectrum for RENO_FD */
  // InitOutput(MYFILE4, "");
  // true_rates = glbGetRuleRatePtr(3, 0);
  // for (i = 0; i <= n_bins - 1; i++)
  //   AddToOutput(bin_centers[i] - 0.78e-3, true_rates[i] * (0.1E-3 / bin_widths[i]), bin_widths[i]);

  // /* Provide the total IBD events for RENO_FD */
  // total_rate = 0.0;
  // total_rate = glbTotalRuleRate(3, 0, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG);
  // printf("Total IBD rate in RENO_FD: %g\n", total_rate);

  // /*****************************
  //  * Produce Daya Bay  spectra *
  //  * ***************************/
  // glbDefineParams(central_values, theta12, DYBth13_min, theta23, deltacp, sdm, DYBdm31_min);
  // glbSetOscillationParameters(central_values);
  // glbSetRates();



  // n_bins = glbGetNumberOfBins(4);
  // bin_centers = glbGetBinCentersListPtr(4);
  // bin_widths = glbGetBinSizeListPtr(4);

  // /* Compute the combined event spectrum for DayaBay_EH1 */
  // InitOutput(MYFILE5, "");
  // true_rates = glbGetRuleRatePtr(4, 0);
  // for (i = 0; i <= n_bins - 1; i++)
  //   AddToOutput(bin_centers[i] - 0.80e-3, true_rates[i] * (1E-3 / bin_widths[i]), bin_widths[i]);

  // /* Provide the total IBD events for DayaBay_EH1 */
  // total_rate = 0.0;
  // total_rate = glbTotalRuleRate(4, 0, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG);
  // printf("Total IBD rate in DayaBay_EH1: %g\n", total_rate);
  // // glbFreeParams(true_values);

  // n_bins = glbGetNumberOfBins(5);
  // bin_centers = glbGetBinCentersListPtr(5);
  // bin_widths = glbGetBinSizeListPtr(5);
  // /* Compute the combined event spectrum for DayaBay_EH2 */
  // InitOutput(MYFILE6, "");
  // true_rates = glbGetRuleRatePtr(5, 0);
  // for (i = 0; i <= n_bins - 1; i++)
  //   AddToOutput(bin_centers[i] - 0.8e-3, true_rates[i] * (1E-3 / bin_widths[i]), bin_widths[i]);

  // /* Provide the total IBD events for DayaBay_EH2 */
  // total_rate = 0.0;
  // total_rate = glbTotalRuleRate(5, 0, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG);
  // printf("Total IBD rate in DayaBay_EH2: %g\n", total_rate);

  // n_bins = glbGetNumberOfBins(6);
  // bin_centers = glbGetBinCentersListPtr(6);
  // bin_widths = glbGetBinSizeListPtr(6);
  // /* Compute the combined event spectrum for DayaBay_EH3 */
  // InitOutput(MYFILE7, "");
  // true_rates = glbGetRuleRatePtr(6, 0);
  // for (i = 0; i <= n_bins - 1; i++)
  //   AddToOutput(bin_centers[i] - 0.8e-3, true_rates[i] * (1E-3 / bin_widths[i]), bin_widths[i]);

  // /* Provide the total IBD events for DayaBay_EH3 */
  // total_rate = 0.0;
  // total_rate = glbTotalRuleRate(6, 0, GLB_ALL, GLB_W_EFF, GLB_WO_BG, GLB_W_COEFF, GLB_SIG);
  // printf("Total IBD rate in DayaBay_EH3: %g\n", total_rate);

  /* Destroy parameter vector(s) */
  glbFreeParams(central_values);
  glbFreeParams(input_errors);
  glbFreeParams(test_values);
  glbFreeParams(minimum);
  glbFreeProjection(SM_projection);
  exit(0);
}
