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
    This program is to make our fixed value equal to measured value.
    We are processing DoubleChooze data as a demonstration.
 * ##################################################################### *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include <globes/globes.h> /* GLoBES library */
#include "Lorentz_violation_lat.h"
#include "myio.h" /* my input-output routines */

#define SQR(x) ((x) * (x)) /* macro to calculate squares */

/* Pre-defined indexing for the reactor experiments */
#define EXP_DC_ND 0   /* Double Chooz near detector */
#define EXP_DC_FD 1   /* Double Chooz far detector */
#define EXP_RENO_FD 2 /* RENO far detector */
#define EXP_RENO_ND 3 /* RENO near detector */
#define EXP_DYB_EH1 4 /* DayaBay EH1 detector */
#define EXP_DYB_EH2 5 /* DayaBay EH2 detector */
#define EXP_DYB_EH3 6 /* DayaBay EH3 detector */

/* Gaussian likelihood function */
inline double likelihood(double true_rate, double fit_rate, double sqr_sigma)
{
  if (sqr_sigma > 0)
    return SQR(true_rate - fit_rate) / sqr_sigma;
  else
    return 0.0;
}

#include "DoubleChooz_chi.h"
#include "RENO_chi.h"
#include "DayaBay_chi.h"

int Target_Parameter = GLB_AS_210;
#include "glbmath_lat.h"

char MYFILE0[] = "../Data/all_reactor_div/LIV_BC31_theta13_ALL.dat";
char MYFILE1[] = "../Data/all_reactor_div/LIV_BC31_theta13_DC.dat";
char MYFILE2[] = "../Data/all_reactor_div/LIV_BC31_theta13_RN.dat";
char MYFILE3[] = "../Data/all_reactor_div/LIV_BC31_theta13_DYB.dat";

double xmin = 0.0600;
double xmax = 0.1400;
int xsteps = 800;

double ymin = 2.200e-3;
double ymax = 3.200e-3;
int ysteps = 1000;

int main(int argc, char *argv[])
{
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
  /* Initialize libglobes */
  glbInit(argv[0]);
  /* Initialize LIV engine*/
  glbRegisterProbabilityEngine(19, /*Number of parameters*/
                               &my_probability_matrix,
                               &my_set_oscillation_parameters,
                               &my_get_oscillation_parameters,
                               NULL);
  // glbDefineOscEngine(15, /*Number of parameters*/
  //                              &my_probability_matrix,
  //                              &my_set_oscillation_parameters,
  //                              &my_get_oscillation_parameters,"DC",
  //                              NULL);
  /* Define chi^2 functions for the reactor experiments */
  glbDefineChiFunction(&chiDCNorm, 5,
                       "chiDCNorm", &DoubleChooz_true_rates);
  glbDefineChiFunction(&chiRENONorm, 7, "chiRENONorm", &RENO_true_rates);
  glbDefineChiFunction(&chiDayaBayNorm, 10,
                       "chiDayaBayNorm", &DayaBay_true_rates);

  /* Initialize experiment(s) used in this simulation */
  glbInitExperiment("../GLB/DoubleChooz/DoubleChooz_JPfix_ND.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("../GLB/DoubleChooz/DoubleChooz_JPfix_FD.glb", &glb_experiment_list[0], &glb_num_of_exps);

  glbInitExperiment("../GLB/RENO/RENO_JPfix_FD.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("../GLB/RENO/RENO_JPfix_ND.glb", &glb_experiment_list[0], &glb_num_of_exps);

  // glbInitExperiment("../GLB/RENO/RENO_JPfix_FD.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("../GLB/DayaBay/DayaBay_EH1_2022.glb", &glb_experiment_list[0], &glb_num_of_exps); // 1
  glbInitExperiment("../GLB/DayaBay/DayaBay_EH2_2022.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbInitExperiment("../GLB/DayaBay/DayaBay_EH3_2022.glb", &glb_experiment_list[0], &glb_num_of_exps);

  /* Initialize parameter and projection vector(s) */
  glb_params central_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params minimum = glbAllocParams();
  glb_params minimum_t = glbAllocParams();
  glb_projection LIV_projection = glbAllocProjection();

  /* Set the parameter vector */
  glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, ldm);
  glbDefineParams(input_errors, theta12_error, theta13_error, theta23_error, deltacp_error, sdm_error, ldm_error);
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

  glbSetDensityParams(central_values, 1.0, GLB_ALL);
  glbSetDensityParams(input_errors, 0.05, GLB_ALL);
  /* The oscillation probabilities are computed */
  glbSetOscillationParameters(central_values);
  glbSetRates();
  /* Set central values and Input_error*/
  glbSetCentralValues(central_values);
  glbSetInputErrors(input_errors);

  /*Set up the Projection  */
  glbDefineProjection(LIV_projection, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED);
  glbSetDensityProjectionFlag(LIV_projection, GLB_FIXED, GLB_ALL);
  /*Lorentz parameters*/
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_AS_210);
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_AS_211);
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_AC_210);
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_AC_211);
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_BS_21);
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_BC_21);

  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_AS_310);
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_AS_311);
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_AC_310);
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_AC_311);
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_BS_31);
  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_BC_31);

  glbSetProjectionFlag(LIV_projection, GLB_FIXED, GLB_Lat);

  glbSetProjection(LIV_projection);

  ///////////////////////
  //  TARGET SETTING   //
  ///////////////////////
  /* Initiate a parameter vector for the scan */
  glbCopyParams(central_values, test_values);
  char *LIV_name[] = {"AS210", "AS211", "AC210", "AC211", "BS21", "BC21", "AS310", "AS311", "AC310", "AC311", "BS21", "BC31"};
  // Theta13 analysis
  double th13result_save[4];
  double sigmath13_save[2] = {1, 1};
  double LIVth13cons[2];
  
  LIV_minizer(Collab_th13Chi2, test_values, th13result_save, GLB_TH13_ONLY, GLB_AS_210);
  printf("Theta13 analysis: \nThe Best SMLIV = %g, Sin^22Theta13 = %g, Dmee = %g, Chi2 = %g\n", th13result_save[0], SQR(sin(2 * th13result_save[1])), th13result_save[2] - SQR(sin(theta12)) * sdm, th13result_save[3]);
  LIV_div_minizer(LIV_th13div, test_values, GLB_AS_210, LIVth13cons, SM_MODE);
  printf("SMLIV = %g, Chi2Cons. = %g without LIV prior\n", LIVth13cons[0], LIVth13cons[1]);

  for (int i = 6; i <= 17; i++)
  {
    LIV_sigma(Collab_th13Chi2, test_values, th13result_save, sigmath13_save, GLB_TH13_LIV, i);
    printf("The Best %s = %g +- %g, Sin^22Theta13 = %g, Dmee = %g, Chi2 = %g\n", LIV_name[i - 6], th13result_save[0], fmax(sigmath13_save[0], sigmath13_save[1]), SQR(sin(2 * th13result_save[1])), th13result_save[2] - SQR(sin(theta12)) * sdm, th13result_save[3]);
    LIV_th13bestift[i - 6] = th13result_save[0];
    LIV_th13fitsigma[i - 6] = fmax(sigmath13_save[0], sigmath13_save[1]);
    LIV_div_minizer(LIV_th13div, test_values, GLB_AS_210, LIVth13cons, SM_MODE);
    printf("SM with %s prior, Chi2Cons. = %g \n", LIV_name[i - 6], LIVth13cons[1]);
    LIV_div_minizer(LIV_th13div, test_values, i, LIVth13cons, LIV_MODE);
    printf("%s = %g, Chi2Cons. = %g \n", LIV_name[i - 6], LIVth13cons[0], LIVth13cons[1]);
  }

  // 2D analysis
  double D2result_save[4];
  double sigma2D_save[2]={1,1};
  double LIV2Dcons[2];
  LIV_minizer(Collab_2DChi2, test_values, D2result_save, GLB_2D, GLB_AS_210);
  printf("2D analysis:\nThe Best SMLIV = %g, Sin^22Theta13 = %g, Dmee = %g, Chi2 = %g\n", D2result_save[0], SQR(sin(2 * D2result_save[0])), D2result_save[2] - SQR(sin(theta12)) * sdm, D2result_save[3]);
  LIV_div_minizer(LIV_2Ddiv, test_values, GLB_AS_210, LIV2Dcons, SM_MODE);
  printf("SMLIV = %g, Chi2Cons. = %g without LIV prior\n", LIV2Dcons[0], LIV2Dcons[1]);
  for (int i = 6; i <= 17; i++)
  {
    LIV_sigma(Collab_2DChi2, test_values, D2result_save, sigma2D_save, GLB_TH13_LIV, i);
    printf("The Best %s = %g +- %g, Sin^22Theta13 = %g, Dmee = %g, Chi2 = %g\n", LIV_name[i - 6], D2result_save[0], fmax(sigma2D_save[0], sigma2D_save[1]), SQR(sin(2 * D2result_save[1])), D2result_save[2] - SQR(sin(theta12)) * sdm, D2result_save[3]);
    LIV_2Dbestift[i - 6] = D2result_save[0];
    LIV_2Dfitsigma[i - 6] = fmax(sigma2D_save[0], sigma2D_save[1]);
    LIV_div_minizer(LIV_2Ddiv, test_values, i, LIV2Dcons, SM_MODE);
    printf("SM with %s prior, Chi2Cons. = %g \n", LIV_name[i - 6], LIV2Dcons[1]);
    LIV_div_minizer(LIV_2Ddiv, test_values, i, LIV2Dcons, LIV_MODE);
    printf("%s = %g, Chi2Cons. = %g \n", LIV_name[i - 6], LIV2Dcons[0], LIV2Dcons[1]);
  }

  glbFreeParams(central_values);
  glbFreeParams(test_values);
  glbFreeParams(input_errors);
  glbFreeParams(minimum);
  glbFreeProjection(LIV_projection);

  exit(0);
}