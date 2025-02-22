#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>
#include <globes/globes.h> /* GLoBES library */


#define SQR(x) ((x) * (x)) /* macro to calculate squares */

// input locations information of each experiment
double polar_angle_DC = 1.8623;
double polar_angle_RN = 2.2546;
double polar_angle_DYB = 0.73596;
/*********************************************************
 * Chi2 function for each exp with different polar angle *
 *********************************************************/
double DC_Chi2(const gsl_vector *v, void *params)
{
  // Get value from gsl_vector
  double th13,dm31, liv_value;
  double res;
  liv_value = gsl_vector_get(v, 0);
  th13 = gsl_vector_get(v, 1);
  // Set value for a temporary glb vector
  glb_params Input_values = (glb_params)params;
  glb_params test_values = glbAllocParams();
  glbCopyParams(Input_values, test_values);
  double theta12,sdm;
  theta12 = glbGetOscParams(test_values,GLB_THETA_12);
  sdm = glbGetOscParams(test_values,GLB_DM_21);
  dm31 = 2.519e-3 + SQR(sin(theta12)) * sdm;
  glbSetOscParams(test_values, liv_value, Target_Parameter);
  glbSetOscParams(test_values, th13, GLB_THETA_13);
  glbSetOscParams(test_values, dm31, GLB_DM_31);
  glbSetOscParams(test_values, polar_angle_DC, GLB_Lat);
  glb_params minimum = glbAllocParams();
  // Obtain Chi2 of DoubleChooz
  res = glbChiNP(test_values, minimum, GLB_ALL);
  res = glbChiNP(test_values, minimum, EXP_DC_FD);
  glbSetSysStartingValuesList(EXP_DC_FD, 0, GLB_ON, sys_dcstartval);
  glbFreeParams(test_values);
  glbFreeParams(minimum);
  return res;
}
double RN_Chi2(const gsl_vector *v, void *params)
{
  double th13, dm31, liv_value;
  double res;
  liv_value = gsl_vector_get(v, 0);
  th13 = gsl_vector_get(v, 1);
  dm31 = gsl_vector_get(v, 2);

  glb_params Input_values = (glb_params)params;
  glb_params test_values = glbAllocParams();
  glbCopyParams(Input_values, test_values);
  glbSetOscParams(test_values, liv_value, Target_Parameter);
  glbSetOscParams(test_values, th13, GLB_THETA_13);
  glbSetOscParams(test_values, dm31, GLB_DM_31);
  glbSetOscParams(test_values, polar_angle_RN, GLB_Lat);
  glb_params minimum = glbAllocParams();

  res = glbChiNP(test_values, minimum, GLB_ALL);
  res = glbChiNP(test_values, minimum, EXP_RENO_FD);
  glbSetSysStartingValuesList(EXP_RENO_FD, 0, GLB_ON, sys_RNstartval);
  glbFreeParams(test_values);
  glbFreeParams(minimum);
  return res;
}
double DYB_Chi2(const gsl_vector *v, void *params)
{
  double th13, dm31, liv_value;
  double res;
  liv_value = gsl_vector_get(v, 0);
  th13 = gsl_vector_get(v, 1);
  dm31 = gsl_vector_get(v, 2);

  glb_params Input_values = (glb_params)params;
  glb_params test_values = glbAllocParams();
  glbCopyParams(Input_values, test_values);
  glbSetOscParams(test_values, liv_value, Target_Parameter);
  glbSetOscParams(test_values, th13, GLB_THETA_13);
  glbSetOscParams(test_values, dm31, GLB_DM_31);
  glbSetOscParams(test_values, polar_angle_DYB, GLB_Lat);
  glb_params minimum = glbAllocParams();

  res = glbChiNP(test_values, minimum, GLB_ALL);
  res = glbChiNP(test_values, minimum, EXP_DYB_EH1);
  glbSetSysStartingValuesList(EXP_DYB_EH1, 0, GLB_ON, sys_dybstartval);
  glbFreeParams(test_values);
  glbFreeParams(minimum);
  return res;
}

//Standard goodness of git, run Livminizer with GLB_2D 
double SGChi2(const gsl_vector *v, void *params){
  double res_DC, res_DYB, res_RN, res;
  res_DC = DC_Chi2(v, params);
  res_RN = RN_Chi2(v, params);
  res_DYB = DYB_Chi2(v, params);

  res = res_DYB + res_RN + res_DC;
  // res = res_RN  + res_DYB;
  return res;
}

/*****************************
 * LIV Minizer function      *
 * vec[0]---LIV value        *
 * vec[1]---Theta 13         *
 * vec[2]---dm31 NH          *
 * out[0-2]---vec[0-2]       *
 * out[3]---Collab_Chi2      *
 *****************************/
#define GLB_SM -1 // th13 + dmee
#define GLB_LIV -2  // liv + th13 + dmee 
void LIV_minizer(double (*chi2_func)(const gsl_vector *v, void *params),
                 glb_params test_values, double *out, int GLB_CHOICE, int LIV_target)
{
  Target_Parameter = LIV_target;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  // Oscillation Parameter scale evalution
  const double param_scales[] = {1, 0.1, 2.5e-3}; // [liv_value, th13, dm31]
  const double step_ratio = 0.1;
  const size_t total_dim = 3;

  // Minizer vector setting up
  ss = gsl_vector_alloc(total_dim);
  x = gsl_vector_alloc(total_dim);

  // Get initial value from test_value
  glb_params initial_params = glbAllocParams();
  glbCopyParams(test_values, initial_params);
  double liv_init = glbGetOscParams(initial_params, Target_Parameter);
  double th13_init = glbGetOscParams(initial_params, GLB_THETA_13);
  double dm31_init = glbGetOscParams(initial_params, GLB_DM_31);
  glbFreeParams(initial_params);

  // Switch mode
  switch (GLB_CHOICE)
  {
  case GLB_SM: // th13 + dm31
    // liv_init = 0.0;
    for (int i = 0; i < total_dim; i++)
      gsl_vector_set(ss, i, (i > 0) ? param_scales[i] * step_ratio : 0.0);
    break;

  case GLB_LIV: // ALL
    for (int i = 0; i < total_dim; i++)
      gsl_vector_set(ss, i, param_scales[i] * step_ratio);
    break;

  default:
    fprintf(stderr, "[ERROR] Undefined GLB_CHOICE=%d\n", GLB_CHOICE);
    exit(EXIT_FAILURE); // Break up
  }

  // Initial vector
  gsl_vector_set(x, 0, liv_init);
  gsl_vector_set(x, 1, th13_init);
  gsl_vector_set(x, 2, dm31_init);

  // Minizer configuration
  minex_func.n = total_dim;
  minex_func.f = chi2_func;
  minex_func.params = test_values;
  s = gsl_multimin_fminimizer_alloc(T, total_dim);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  // Minimum iteration
  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status)
      break;
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-7);
  } while (status == GSL_CONTINUE && iter < 100);

  // Output
  out[0] = gsl_vector_get(s->x, 0);
  out[1] = gsl_vector_get(s->x, 1);
  out[2] = gsl_vector_get(s->x, 2);
  out[3] = s->fval;

  // Free vector
  gsl_vector_free(ss);
  gsl_vector_free(x);
  gsl_multimin_fminimizer_free(s);
}
double PGchi2(glb_params test_values,int GLB_Parameter){
  // 2D analysis
  double D2result_save[4];
  double LIV2Dcons[2];
  double chi2DC,chi2RN,chi2DYB;
  double chi2Combined;
  double chi2PG;
  if(GLB_Parameter == GLB_SM){// GLB_AS_210 makes no sense here, Chi2_combined in this case
      LIV_minizer(SGChi2, test_values, D2result_save, GLB_SM, GLB_AS_210);
      chi2Combined = D2result_save[3];
      // printf("SGchi2 = %g\n",chi2DC);
      LIV_minizer(DC_Chi2, test_values, D2result_save, GLB_SM, GLB_AS_210);
      chi2DC = D2result_save[3];
      // printf("DCchi2 = %g\n",chi2DC);
      LIV_minizer(RN_Chi2, test_values, D2result_save, GLB_SM, GLB_AS_210);
      chi2RN = D2result_save[3];
      LIV_minizer(DYB_Chi2, test_values, D2result_save, GLB_SM, GLB_AS_210);
      chi2DYB = D2result_save[3];
      chi2PG = chi2Combined - chi2DC- chi2RN- chi2DYB;
      
  }
  else if(GLB_Parameter>=6&&GLB_Parameter<=17){
    LIV_minizer(SGChi2, test_values, D2result_save, GLB_LIV, GLB_Parameter);
    chi2Combined = D2result_save[3];
    LIV_minizer(DC_Chi2, test_values, D2result_save, GLB_LIV, GLB_Parameter);
    chi2DC = D2result_save[3];
    LIV_minizer(RN_Chi2, test_values, D2result_save, GLB_LIV, GLB_Parameter);
    chi2RN = D2result_save[3];
    LIV_minizer(DYB_Chi2, test_values, D2result_save, GLB_LIV, GLB_Parameter);
    chi2DYB = D2result_save[3];
    chi2PG = chi2Combined - chi2DC- chi2RN- chi2DYB;
  }
  return chi2PG;
}

