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
  glbSetOscParams(test_values, polar_angle_DYB, GLB_Lat);
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
#define GLB_TH13_ONLY -1 // th13
#define GLB_TH13_LIV -2  // liv + th13
#define GLB_2D -3        // th13 + dm31
#define GLB_2D_LIV -4    // liv + th13 + dm31
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
  case GLB_TH13_ONLY: // th13
    // liv_init = 0.0;   // liv_value = 0
    for (int i = 0; i < total_dim; i++)
      gsl_vector_set(ss, i, (i == 1) ? param_scales[i] * step_ratio : 0.0);
    break;

  case GLB_TH13_LIV: // liv_value + th13
    for (int i = 0; i < total_dim; i++)
      gsl_vector_set(ss, i, (i < 2) ? param_scales[i] * step_ratio : 0.0);
    break;

  case GLB_2D: // th13 + dm31
    // liv_init = 0.0;
    for (int i = 0; i < total_dim; i++)
      gsl_vector_set(ss, i, (i > 0) ? param_scales[i] * step_ratio : 0.0);
    break;

  case GLB_2D_LIV: // ALL
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

/* 辅助结构体用于传递固定参数 */
struct root_params
{
  glb_params params;
  double (*chi2_func)(const gsl_vector *, void *);
  double target_chi2;
  int LIV_target;
  double *bestfit_LIV;
  int direction;
};

/* 修改后的约束检查函数 */
static double liv_chi2_diff(const gsl_vector *v, void *params)
{
  struct root_params *p = (struct root_params *)params;

  /* 新增：方向约束检查 */
  const double current_LIV = gsl_vector_get(v, 0);
  if ((p->direction == -1 && (current_LIV > p->bestfit_LIV[0] + 1e-8)) || // 下界搜索禁止右移
      (p->direction == 1 && (current_LIV < p->bestfit_LIV[0] - 1e-8)))    // 上界搜索禁止左移
  {
    return 1e20 * fabs(current_LIV - p->bestfit_LIV[0]); // 指数级惩罚
  }

  glb_params test_values1 = glbAllocParams();
  glbCopyParams(p->params, test_values1);

  /* 计算卡方值 */
  double current_chi2 = p->chi2_func(v, test_values1);

  /* 清理资源 */
  glbFreeParams(test_values1);

  return fabs(current_chi2 - p->target_chi2);
}
/*****************************
 * LIV Sigma function        *
 * out---LIV minizer out[3]  *
 * sigma[0]---sigma_low      *
 * sigma[1]---sigma_high     *
 *****************************/

/* 主函数：计算LIV参数的1σ边界 */
void LIV_sigma(double (*chi2_func)(const gsl_vector *, void *),
               glb_params test_values, double *bestfit,
               double *sigma, int GLB_CHOICE, int LIV_target)
{
  Target_Parameter = LIV_target;
  double result_save[4];
  LIV_minizer(chi2_func, test_values, result_save, GLB_CHOICE, Target_Parameter);
  for (int i = 0; i < 4; i++)
  {
    bestfit[i] = result_save[i];
  }

  const double param_scales[] = {1, 0.1, 2.5e-3};
  const size_t total_dim = 3;
  const int max_attempts = 3;
  const double convergence_eps = 1e-5;

  /* 存储边界结果 */
  double bounds[2] = {bestfit[0], bestfit[0]};
  double target_chi2 = bestfit[3] + 1.0;

  /* 双向搜索循环 */
  for (int direction = -1; direction <= 1; direction += 2)
  {
    if (direction == 0)
      continue;

    /* 配置方向约束参数 */
    struct root_params r_params = {
        .params = test_values,
        .chi2_func = chi2_func,
        .target_chi2 = target_chi2,
        .LIV_target = LIV_target,
        .direction = direction, // 新增方向标记
        .bestfit_LIV = bestfit  // 新增最佳值
    };

    /* 初始化优化器 */
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(
        gsl_multimin_fminimizer_nmsimplex2, total_dim);

    /* 参数初始化 */
    gsl_vector *x = gsl_vector_alloc(total_dim);
    gsl_vector_set(x, 0, bestfit[0]);
    gsl_vector_set(x, 1, bestfit[1]);
    gsl_vector_set(x, 2, bestfit[2]);

    /* 动态步长调整 */
    gsl_vector *ss = gsl_vector_alloc(total_dim);
    const double init_step = direction * fmax(fabs(bestfit[0] * 0.1), param_scales[0] * 0.1);
    gsl_vector_set(ss, 0, init_step); // 主参数步长
    gsl_vector_set(ss, 1, 0);         // th13固定
    gsl_vector_set(ss, 2, 0);         // dm31固定

    /* 配置优化参数 */
    gsl_multimin_function minex_func = {
        .n = total_dim,
        .f = liv_chi2_diff,
        .params = &r_params};

    /* 迭代优化流程 */
    int status = GSL_CONTINUE;
    size_t iter = 0;
    double size;
    double final_value = bestfit[0];
    double prev_value = bestfit[0];
    int stagnation_count = 0;

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

    /* 结果验证和处理 */
    if (status == GSL_SUCCESS)
    {
      const double found_value = gsl_vector_get(s->x, 0);
      const double delta_chi2 = liv_chi2_diff(s->x, &r_params);

      /* 二次验证结果有效性 */
      if (fabs(delta_chi2) < 0.2 &&
          ((direction == -1 && found_value < bestfit[0]) ||
           (direction == 1 && found_value > bestfit[0])))
      {
        bounds[(direction > 0) ? 1 : 0] = found_value;
      }
    }

    /* 释放当前方向资源 */
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
  }

  /* 最终结果处理 */
  sigma[0] = (bounds[0] < bestfit[0]) ? (bestfit[0] - bounds[0]) : 0.0;
  sigma[1] = (bounds[1] > bestfit[0]) ? (bounds[1] - bestfit[0]) : 0.0;
  // printf("sigma0 = %g, sigma1 = %g\n", sigma[0], sigma[1]);
  /* 错误处理：当两边搜索都失败时 */
  if (sigma[0] <= 0 && sigma[1] <= 0)
  {
    fprintf(stderr, "LIV_sigma: Both direction searches failed!\n");
    sigma[0] = sigma[1] = 0.0;
  }
}

double LIV_th13bestift[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double LIV_th13fitsigma[12] = {1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7};
double LIV_2Dbestift[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double LIV_2Dfitsigma[12] = {1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e7};

#define SM_MODE -5  // Chi2TH13cons_cal
#define LIV_MODE -6 // Chi2D2cons_cal

double Chi2TH13cons_cal(double DC_th13, double RN_th13, double DYB_th13)
{
  double DC_th13_sigma, RN_th13_sigma, DYB_th13_sigma;
  DC_th13_sigma = 0.012; // standard deviation from reference
  RN_th13_sigma = 0.009;
  DYB_th13_sigma = 0.0024;
  double div1, div2, div3, sin22fixth13;
  sin22fixth13 = (SQR(sin(2 * DC_th13)) / SQR(DC_th13_sigma) + SQR(sin(2 * RN_th13)) / SQR(RN_th13_sigma) + SQR(sin(2 * DYB_th13)) / SQR(DYB_th13_sigma)) / (1 / SQR(DC_th13_sigma) + 1 / SQR(RN_th13_sigma) + 1 / SQR(DYB_th13_sigma));
  div1 = SQR(SQR(sin(2 * DC_th13)) - sin22fixth13) / SQR(DC_th13_sigma);
  div2 = SQR(SQR(sin(2 * RN_th13)) - sin22fixth13) / SQR(RN_th13_sigma);
  div3 = SQR(SQR(sin(2 * DYB_th13)) - sin22fixth13) / SQR(DYB_th13_sigma);
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
  double sin22fixth13, fixdm31;
  double div1, div2;
  sin22fixth13 = (SQR(sin(2 * RN_th13)) / SQR(RN_th13_sigma) + SQR(sin(2 * DYB_th13)) / SQR(DYB_th13_sigma)) / (1 / SQR(RN_th13_sigma) + 1 / SQR(DYB_th13_sigma));
  fixdm31 = (RN_dm31 / SQR(RN_dm31_sigma) + DYB_dm31 / SQR(DYB_dm31_sigma)) / (1 / SQR(RN_dm31_sigma) + 1 / SQR(DYB_dm31_sigma));
  div1 = SQR(SQR(sin(2 * RN_th13)) - sin22fixth13) / SQR(RN_th13_sigma) + SQR(RN_dm31 - fixdm31) / SQR(RN_dm31_sigma);
  div2 = SQR(SQR(sin(2 * DYB_th13)) - sin22fixth13) / SQR(DYB_th13_sigma) + SQR(DYB_dm31 - fixdm31) / SQR(DYB_dm31_sigma);
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

  double div_liv, result[4];
  double th13_DC, th13_RN, th13_DYB;
  double chi2_DC, chi2_RN, chi2_DYB;
  LIV_minizer(DC_Chi2, test_values, result, GLB_TH13_ONLY, Target_Parameter);
  th13_DC = result[1];
  chi2_DC = result[3];
  // printf("DC =%g\t %g\t%g\n",result[0], result[1],result[3]);
  LIV_minizer(RN_Chi2, test_values, result, GLB_TH13_ONLY, Target_Parameter);
  th13_RN = result[1];
  chi2_RN = result[3];
  // printf("RN =%g\t %g\t%g\n",result[0], result[1],result[3]);
  LIV_minizer(DYB_Chi2, test_values, result, GLB_TH13_ONLY, Target_Parameter);
  th13_DYB = result[1];
  chi2_DYB = result[3];
  // printf("DYB =%g\t %g\t%g\n",result[0], result[1],result[3]);

  div_liv = Chi2TH13cons_cal(th13_DC, th13_RN, th13_DYB)+ SQR(LIV_parameter - LIV_th13bestift[Target_Parameter - 6]) / SQR(LIV_th13fitsigma[Target_Parameter - 6]);

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

  double div_liv, result[4];
  double th13_RN, th13_DYB, dm31_RN, dm31_DYB;
  double chi2_DC, chi2_RN, chi2_DYB;
  LIV_minizer(RN_Chi2, test_values, result, GLB_2D, Target_Parameter);
  th13_RN = result[1];
  dm31_RN = result[2];
  chi2_RN = result[3];
  LIV_minizer(DYB_Chi2, test_values, result, GLB_2D, Target_Parameter);
  th13_DYB = result[1];
  dm31_DYB = result[2];
  chi2_DYB = result[3];

  div_liv = Chi2D2cons_cal(th13_RN, dm31_RN, th13_DYB, th13_DYB) + SQR(LIV_parameter - LIV_2Dbestift[Target_Parameter - 6]) / SQR(LIV_2Dfitsigma[Target_Parameter - 6]);

  glbFreeParams(test_values);
  glbFreeParams(minimum);

  return div_liv;
}
double LIV_div_minizer(double (*chi2_func)(const gsl_vector *v, void *params), glb_params test_values, int LIV_target, double *output, int GLB_CHOICE)
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
  if (GLB_CHOICE == SM_MODE) // test_value with LIV = 0
  {
    gsl_vector_set_all(ss, 0);
    gsl_vector_set_all(x, 0);
  }
  else if (GLB_CHOICE == LIV_MODE)
  {
    gsl_vector_set_all(ss, 0.1);
    gsl_vector_set_all(x, 0);
  }
  else
  {
    printf("There is no such mode\n");
    exit(EXIT_FAILURE); // Break up;
  }

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

    // if (status == GSL_SUCCESS)
    // {
    //   printf("Minimum found at:\n");

    //   printf("%5ld %.10f %10.5f\n", iter,
    //          gsl_vector_get(s->x, 0),
    //          s->fval);
    // }
  } while (status == GSL_CONTINUE && iter < 100);
  double LIV_parameter = 0, div = 0;
  LIV_parameter = gsl_vector_get(s->x, 0);
  div = s->fval;
  output[0] = LIV_parameter;
  output[1] = div;
}
