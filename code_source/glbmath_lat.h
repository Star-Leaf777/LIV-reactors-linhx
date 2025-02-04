#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_multimin.h>
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
  double th13, liv_value;
  double res;
  liv_value = gsl_vector_get(v, 0);
  th13 = gsl_vector_get(v, 1);
  // Set value for a temporary glb vector
  glb_params Input_values = (glb_params)params;
  glb_params test_values = glbAllocParams();
  glbCopyParams(Input_values, test_values);
  glbSetOscParams(test_values, liv_value, Target_Parameter);
  glbSetOscParams(test_values, th13, GLB_THETA_13);
  glbSetOscParams(test_values, polar_angle_DC, GLB_Lat);
  glb_params minimum = glbAllocParams();
  // Obtain Chi2 of DoubleChooz
  res = glbChiNP(test_values, minimum, GLB_ALL);
  res = glbChiNP(test_values, minimum, EXP_DC_FD);
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
  glbSetOscParams(test_values, polar_angle_RN, GLB_Lat);
  glb_params minimum = glbAllocParams();

  res = glbChiNP(test_values, minimum, GLB_ALL);
  res = glbChiNP(test_values, minimum, EXP_DYB_EH1);
  glbFreeParams(test_values);
  glbFreeParams(minimum);
  return res;
}
// Collaboration Chi2 function
double Collab_th13Chi2(const gsl_vector *v, void *params)
{
  double res_DC, res_DYB, res_RN, res;
  res_DC = DC_Chi2(v, params);
  res_RN = RN_Chi2(v, params);
  res_DYB = DYB_Chi2(v, params);

  res = res_DYB + res_RN + res_DC;
  // res = res_RN  + res_DYB;
  return res;
}
// Collaboration Chi2 function
double Collab_2DChi2(const gsl_vector *v, void *params)
{
  double res_DYB, res_RN, res;
  res_RN = RN_Chi2(v, params);
  res_DYB = DYB_Chi2(v, params);

  res = res_DYB + res_RN;
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
#define GLB_TH13_ONLY   -1  // 仅优化 th13
#define GLB_TH13_LIV    -2  // 优化 th13 + liv_value
#define GLB_2D          -3  // 优化 th13 + dm31
#define GLB_2D_LIV      -4  // 优化全部三个参数
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

  // 参数物理特性配置
  const double param_scales[] = {1, 0.1, 2.5e-3};  // [liv_value, th13, dm31]
  const double step_ratio = 0.1;
  const size_t total_dim = 3;

  // 初始化参数向量
  ss = gsl_vector_alloc(total_dim);
  x = gsl_vector_alloc(total_dim);
  
  // 从 test_values 获取初始值（后续可能覆盖）
  glb_params initial_params = glbAllocParams();
  glbCopyParams(test_values, initial_params);
  double liv_init = glbGetOscParams(initial_params, Target_Parameter);
  double th13_init = glbGetOscParams(initial_params, GLB_THETA_13);
  double dm31_init = glbGetOscParams(initial_params, GLB_DM_31);
  glbFreeParams(initial_params);

  // 根据模式调整初始值和步长
  switch (GLB_CHOICE) {
    case GLB_TH13_ONLY:  // 仅优化 th13
      liv_init = 0.0;    // 强制固定 liv_value
      for (int i=0; i<total_dim; i++)
        gsl_vector_set(ss, i, (i == 1) ? param_scales[i]*step_ratio : 0.0);
      break;

    case GLB_TH13_LIV:   // 优化 th13 + liv_value
      for (int i=0; i<total_dim; i++)
        gsl_vector_set(ss, i, (i < 2) ? param_scales[i]*step_ratio : 0.0);
      break;

    case GLB_2D:         // 优化 th13 + dm31
      liv_init = 0.0;    // 强制固定 liv_value
      for (int i=0; i<total_dim; i++)
        gsl_vector_set(ss, i, (i > 0) ? param_scales[i]*step_ratio : 0.0);
      break;

    case GLB_2D_LIV:     // 优化全部三个参数
      for (int i=0; i<total_dim; i++)
        gsl_vector_set(ss, i, param_scales[i]*step_ratio);
      break;

    default:
      fprintf(stderr, "[ERROR] Undefined GLB_CHOICE=%d\n", GLB_CHOICE);
      exit(EXIT_FAILURE); // 直接终止程序
  }

  // 设置初始参数向量（注意覆盖liv_init）
  gsl_vector_set(x, 0, liv_init);
  gsl_vector_set(x, 1, th13_init);
  gsl_vector_set(x, 2, dm31_init);

  // 配置优化器
  minex_func.n = total_dim;
  minex_func.f = chi2_func;
  minex_func.params = test_values;
  s = gsl_multimin_fminimizer_alloc(T, total_dim);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  // 优化迭代（与原始逻辑一致）
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status) break;
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-7);
  } while (status == GSL_CONTINUE && iter < 100);

  // 输出结果
  out[0] = gsl_vector_get(s->x, 0);
  out[1] = gsl_vector_get(s->x, 1);
  out[2] = gsl_vector_get(s->x, 2);
  out[3] = s->fval;

  // 释放资源
  gsl_vector_free(ss);
  gsl_vector_free(x);
  gsl_multimin_fminimizer_free(s);
}

#define TH13_MODE -5 // Chi2TH13cons_cal
#define D2_MODE -6 // Chi2D2cons_cal

double Chi2TH13cons_cal(double DC_th13, double RN_th13, double DYB_th13)
{
  double DC_th13_sigma, RN_th13_sigma, DYB_th13_sigma;
  DC_th13_sigma = 0.012;// standard deviation from reference
  RN_th13_sigma = 0.009;
  DYB_th13_sigma = 0.0024;
  double div1, div2, div3, fixth13;
  fixth13 = (DC_th13 / SQR(DC_th13_sigma) + RN_th13 / SQR(RN_th13_sigma) + DYB_th13 / SQR(DYB_th13_sigma)) / (1 / SQR(DC_th13_sigma) + 1 / SQR(RN_th13_sigma) + 1 / SQR(DYB_th13_sigma));
  div1 = SQR(DC_th13 - fixth13) / SQR(DC_th13_sigma);
  div2 = SQR(RN_th13 - fixth13) / SQR(RN_th13_sigma);
  div3 = SQR(DYB_th13 - fixth13) / SQR(DYB_th13_sigma);
  return div1 + div2 + div3;
}
double Chi2D2cons_cal(double RN_th13, double RN_dm31, double DYB_th13, double DYB_dm31)
{
  double RN_th13_sigma, DYB_th13_sigma;
  double RN_dm31_sigma, DYB_dm31_sigma;
  RN_th13_sigma = 0.009;
  RN_dm31_sigma = 0.16e-3;
  DYB_th13_sigma = 0.0024;
  DYB_dm31_sigma = 0.06e-3;
  double fixth13,fixdm31;
  double div1, div2;
  fixth13 = (RN_th13 / SQR(RN_th13_sigma) + DYB_th13 / SQR(DYB_th13_sigma)) / (1 / SQR(RN_th13_sigma) + 1 / SQR(DYB_th13_sigma));
  fixdm31 = (RN_dm31 / SQR(RN_dm31_sigma) + DYB_dm31 / SQR(DYB_dm31_sigma)) / (1 / SQR(RN_dm31_sigma) + 1 / SQR(DYB_dm31_sigma));
  div1 = SQR(RN_th13 - fixth13) / SQR(RN_th13_sigma)+SQR(RN_dm31 - fixdm31) / SQR(RN_dm31_sigma);
  div2 = SQR(DYB_th13 - fixth13) / SQR(DYB_th13_sigma)+SQR(DYB_dm31 - fixdm31) / SQR(DYB_dm31_sigma);
  return div1 + div2;
}

double LIV_th13div(const gsl_vector *v, void *params)
{
  glb_params Input_values = (glb_params)params;
  glb_params test_values = glbAllocParams();
  glbCopyParams(Input_values, test_values);
  glb_params minimum = glbAllocParams();


  double LIV_parameter;
  LIV_parameter = gsl_vector_get(v, 0);
  glbSetOscParams(test_values, LIV_parameter, Target_Parameter);

  double div_liv,result[4];
  double th13_DC, th13_RN, th13_DYB;
  double chi2_DC, chi2_RN, chi2_DYB;
  LIV_minizer(DC_Chi2, test_values, result, GLB_TH13_ONLY,Target_Parameter);
  th13_DC = result[1];chi2_DC = result[3];
  LIV_minizer(RN_Chi2, test_values, result, GLB_TH13_ONLY,Target_Parameter);
  th13_RN = result[1];chi2_RN = result[3];
  LIV_minizer(DYB_Chi2, test_values, result, GLB_TH13_ONLY,Target_Parameter);
  th13_DYB = result[1];chi2_DYB = result[3];

  div_liv = Chi2TH13cons_cal(th13_DC, th13_RN, th13_DYB);

  glbFreeParams(test_values);
  glbFreeParams(minimum);

  return div_liv;
}
double LIV_2Ddiv(const gsl_vector *v, void *params)
{
  glb_params Input_values = (glb_params)params;
  glb_params test_values = glbAllocParams();
  glbCopyParams(Input_values, test_values);
  glb_params minimum = glbAllocParams();


  double LIV_parameter;
  LIV_parameter = gsl_vector_get(v, 0);
  glbSetOscParams(test_values, LIV_parameter, Target_Parameter);

  double div_liv,result[4];
  double th13_RN, th13_DYB, dm31_RN, dm31_DYB;
  double chi2_DC, chi2_RN, chi2_DYB;
  LIV_minizer(RN_Chi2, test_values, result, GLB_2D,Target_Parameter);
  th13_RN = result[1];dm31_RN = result[2];chi2_RN = result[3];
  LIV_minizer(DYB_Chi2, test_values, result, GLB_2D,Target_Parameter);
  th13_DYB = result[1];dm31_DYB = result[2];chi2_DYB = result[3];

  div_liv = Chi2D2cons_cal(th13_RN,dm31_RN,th13_DYB,th13_DYB);

  glbFreeParams(test_values);
  glbFreeParams(minimum);

  return div_liv;
}
double LIV_div_minizer(double (*chi2_func)(const gsl_vector *v, void *params),glb_params test_values, int LIV_target, double *output)
{
  Target_Parameter = LIV_target;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  // 初始化方法和参数
  size_t n = 1;
  ss = gsl_vector_alloc(n);
  x = gsl_vector_alloc(n);
  gsl_vector_set_all(ss, 0.1);
  gsl_vector_set_all(x, 0);

  minex_func.n = n;
  minex_func.f = chi2_func;
  minex_func.params = test_values;

  s = gsl_multimin_fminimizer_alloc(T, n);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  // 迭代求解
  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status)
      break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-7);

    if (status == GSL_SUCCESS)
      printf("Minimum found at:\n");

    printf("%5ld %.10f %10.5f\n", iter,
           gsl_vector_get(s->x, 0),
           s->fval);
  } while (status == GSL_CONTINUE && iter < 100);
  double LIV_parameter = 0, div = 0;
  LIV_parameter = gsl_vector_get(s->x, 0);
  div = s->fval;
  output[0] = LIV_parameter;
  output[1] = div;
}

