
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "glb_error.h"
#include "glb_wrapper.h"
#include "glb_minimize.h"
#include "glb_probability.h"
//#include "DIAG_FUNC.h"

  

#define GLB_AS_210 6
#define GLB_AS_211 7
#define GLB_AC_210 8
#define GLB_AC_211 9
#define GLB_BS_21 10
#define GLB_BC_21 11

#define GLB_AS_310 12
#define GLB_AS_311 13
#define GLB_AC_310 14
#define GLB_AC_311 15
#define GLB_BS_31 16
#define GLB_BC_31 17

#define GLB_Lat 18

     

/***************************************************************************
 *     U S E R - D E F I N E D   P R O B A B I L I T Y   E N G I N E       *
 ***************************************************************************/

/* Fundamental oscillation parameters */
static double th12, th13, th23; // Mixing angles
static double delta;            // Dirac CP phase
static double mq[3];            // Squared masses static
static double sigma_E;
double deltacp;
static double sdm;
static double ldm;
static double beta21,beta31,beta211,beta311,beta212,beta312;//Lorentz isotopic

static double lat;
static double as210,as211,ac210,ac211,bs21,bc21;//Lorntz noisotopic
static double as310,as311,ac310,ac311,bs31,bc31;


/***************************************************************************
 * Store oscillation parameters in internal data structures.               *
 * For more sophisticated probability engines, this would be the right     *
 * place to pre-compute the mixing matrix and parts of the Hamiltonian in  *
 * order to speed up the calls to the actual probability matrix function.  *
 ***************************************************************************/
int my_set_oscillation_parameters(glb_params p, void *user_data)
{
  th12    = glbGetOscParams(p, GLB_THETA_12);
  th13    = glbGetOscParams(p, GLB_THETA_13);
  th23    = glbGetOscParams(p, GLB_THETA_23);
  deltacp = glbGetOscParams(p, GLB_DELTA_CP);
  sdm     = glbGetOscParams(p, GLB_DM_21) * 1.0e-18;   /* Convert to GeV^2 */
  ldm     = glbGetOscParams(p, GLB_DM_31) * 1.0e-18;   /* Convert to GeV^2 */
/* Lorentz PARAMETERS*/
  as210 =   glbGetOscParams(p, GLB_AS_210);
  as211 =   glbGetOscParams(p, GLB_AS_211);
  ac210 =   glbGetOscParams(p, GLB_AC_210);
  ac211 =   glbGetOscParams(p, GLB_AC_211);
  bs21 =   glbGetOscParams(p, GLB_BS_21);
  bc21 =   glbGetOscParams(p, GLB_BC_21);

  as310 =   glbGetOscParams(p, GLB_AS_310);
  as311 =   glbGetOscParams(p, GLB_AS_311);
  ac310 =   glbGetOscParams(p, GLB_AC_310);
  ac311 =   glbGetOscParams(p, GLB_AC_311);
  bs31 =   glbGetOscParams(p, GLB_BS_31);
  bc31 =   glbGetOscParams(p, GLB_BC_31);
/* Location Parameter */
  lat =   glbGetOscParams(p, GLB_Lat); // rad

  return 0;
}

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/
int my_get_oscillation_parameters(glb_params p, void *user_data)
{
  glbSetOscParams(p, th12, GLB_THETA_12);
  glbSetOscParams(p, th13, GLB_THETA_13);
  glbSetOscParams(p, th23, GLB_THETA_23);
  glbSetOscParams(p, deltacp, GLB_DELTA_CP);
  glbSetOscParams(p, sdm*1.0e18, GLB_DM_21);  /* Convert to eV^2 */
  glbSetOscParams(p, ldm*1.0e18, GLB_DM_31);  /* Convert to eV^2 */
/* Lorentz PARAMETERS*/
  glbSetOscParams(p, as210,  GLB_AS_210);
  glbSetOscParams(p, as211,  GLB_AS_211);
  glbSetOscParams(p, ac210,  GLB_AC_210);
  glbSetOscParams(p, ac211,  GLB_AC_211);
  glbSetOscParams(p, bs21,  GLB_BS_21);
  glbSetOscParams(p, bc21,  GLB_BC_21);

  glbSetOscParams(p, as310,  GLB_AS_310);
  glbSetOscParams(p, as311,  GLB_AS_311);
  glbSetOscParams(p, ac310,  GLB_AC_310);
  glbSetOscParams(p, ac311,  GLB_AC_311);
  glbSetOscParams(p, bs31,  GLB_BS_31);
  glbSetOscParams(p, bc31,  GLB_BC_31);

  glbSetOscParams(p, lat,  GLB_Lat); //rad
  return 0;
}


/***************************************************************************
 *                  HERE STARTS THE IMPLEMENTATION OF the CODE             *
 ***************************************************************************/

/***************************************************************************
 * Function glb_hamiltonian_cd                                             *
 ***************************************************************************
 * Calculates the Hamiltonian for neutrinos (cp_sign=1) or antineutrinos   *
 * (cp_sign=-1) with energy E, propagating in matter of density V          *
 * (> 0 even for antineutrinos) and stores the result in H.                *
 ***************************************************************************/
int glb_hamiltonian_cd_pedro(double E, double V, int cp_sign, gsl_matrix_complex *Hp)
{

  gsl_matrix_complex *Up=gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS); /* The vacuum mixing matrix                           */
  gsl_matrix_complex *T = gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS); /*  Auxiliar Matrix                               */
  double inv_E = 1.0 / E;

  double complex (*_Hp)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(Hp, 0, 0);
  double complex (*_Up)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(Up, 0, 0);
    int i, j;
  //assum that the m[0]==ldm
  if(ldm*1.0e+18<0){mq[0] =-ldm*1.0e+18;}else{mq[0] =ldm*1.0e+18;}
  if(ldm*1.0e+18<0){mq[1] =-ldm*1.0e+18 + sdm*1.0e+18;}else{mq[1] =ldm*1.0e+18 + sdm*1.0e+18;}
  if(ldm*1.0e+18<0){mq[2] =(-ldm + ldm)*1.0e+18;}else{mq[2] =(ldm + ldm)*1.0e+18;}



  if (cp_sign > 0)
  {
  _Up[0][0] = cos(th12)*cos(th13);
  _Up[0][1] = sin(th12)*cos(th13);
  _Up[0][2] = sin(th13) * cexp(-I *deltacp);

  _Up[1][0] = -sin(th12)*cos(th23) - cos(th12)*sin(th23)*sin(th13) * cexp(I*deltacp);
  _Up[1][1] =  cos(th12)*cos(th23) - sin(th12)*sin(th23)*sin(th13) * cexp(I*deltacp);
  _Up[1][2] =  sin(th23)*cos(th13);

  _Up[2][0] =  sin(th12)*sin(th23) - cos(th12)*cos(th23)*sin(th13) * cexp(I*deltacp);
  _Up[2][1] = -cos(th12)*sin(th23) - sin(th12)*cos(th23)*sin(th13) * cexp(I*deltacp);
  _Up[2][2] =  cos(th23)*cos(th13);
  }else
  {/* delta_CP -> -delta_CP */
  _Up[0][0] = cos(th12)*cos(th13);
  _Up[0][1] = sin(th12)*cos(th13);
  _Up[0][2] = sin(th13) * cexp(I *deltacp);

  _Up[1][0] = -sin(th12)*cos(th23) - cos(th12)*sin(th23)*sin(th13) * cexp(-I*deltacp);
  _Up[1][1] =  cos(th12)*cos(th23) - sin(th12)*sin(th23)*sin(th13) * cexp(-I*deltacp);
  _Up[1][2] =  sin(th23)*cos(th13);

  _Up[2][0] =  sin(th12)*sin(th23) - cos(th12)*cos(th23)*sin(th13) * cexp(-I*deltacp);
  _Up[2][1] = -cos(th12)*sin(th23) - sin(th12)*cos(th23)*sin(th13) * cexp(-I*deltacp);
  _Up[2][2] =  cos(th23)*cos(th13);
  }


/* Calculate energy independent matrix H * E */
  gsl_matrix_complex_set_zero(Hp);

  for (i=0; i < GLB_NU_FLAVOURS; i++)
    gsl_matrix_complex_set(Hp, i, i, gsl_complex_rect(0.5*mq[i], 0.0));


// The function H(E) should have unity of 10e-9 eV (neV) but E comes as GeV, thus, one shoud account a (10e9)^n+1 for a E^n power of E.

//H_LIV = (AS_0+AS_1 E) sin(lat) + (AC_0+AC_1 E) cos(lat) + BS E sin(lat) + BC E cos(lat)

  //Because the energy of reactor neutrino ~ MeV = 10^-3 GeV, we set 
  
  double AS21 = cp_sign*as210 + as211 * pow((1.0e3)*E,1);
  double AC21 = cp_sign*ac210 + ac211 * pow((1.0e3)*E,1);
  double BS21 = bs21 * pow((1.0e3)*E,1);
  double BC21 = bc21 * pow((1.0e3)*E,1);

  double AS31 = cp_sign*as310 + as311 * pow((1.0e3)*E,1);
  double AC31 = cp_sign*ac310 + ac311 * pow((1.0e3)*E,1);
  double BS31 = bs31 * pow((1.0e3)*E,1);
  double BC31 = bc31 * pow((1.0e3)*E,1);
  
  _Hp[0][0] = _Hp[0][0]*inv_E;

  _Hp[1][1] = _Hp[1][1]*inv_E +(1.0e9)*(1.0e6)*(AS21*sin(lat) + AC21*cos(lat)+ BS21*sin(2*lat)+BC21*cos(2*lat))*(1.0e-17);//10^-17 is the scale of LIV coefficiences 

  _Hp[2][2] = _Hp[2][2]*inv_E +(1.0e9)*(1.0e6)*(AS31*sin(lat) + AC31*cos(lat)+ BS31*sin(2*lat)+BC31*cos(2*lat))*(1.0e-17);

  gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, Hp, Up, /* T=H0.U^\dagger */
                 GSL_COMPLEX_ZERO, T);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, Up, T,             /* H0=U.T */
                 GSL_COMPLEX_ZERO, Hp);

   _Hp[0][0] = _Hp[0][0] + cp_sign*V;


  gsl_matrix_complex_free(T);
  gsl_matrix_complex_free(Up);


  return 0;
}



/***************************************************************************
 * Function glb_S_matrix_cd                                                *
 ***************************************************************************
 * Calculates the S matrix for neutrino oscillations in matter of constant *
 * density using a fast eigenvalue solver optimized to 3x3 matrices.       *
 ***************************************************************************
 * Parameters:                                                             *
 *   E: Neutrino energy                                                    *
 *   L: Baseline                                                           *
 *   V: Matter potential (must be > 0 even for antineutrinos)              *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 ***************************************************************************/
int glb_S_matrix_cd_pedro(double E, double L, double V, int cp_sign, gsl_matrix_complex *Sp)
{

  gsl_matrix_complex *Qp=gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS); /* The vacuum mixing matrix                           */
  gsl_matrix_complex *T0p=gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS); /* The vacuum mixing matrix                           */


   gsl_vector *lambda = NULL;

   lambda=gsl_vector_alloc (GLB_NU_FLAVOURS);

  /* Introduce some abbreviations */
  double complex (*_Sp)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(Sp,0,0);
  double complex (*_Qp)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(Qp,0,0);
  double complex (*_T0p)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(T0p,0,0);
  double *_lambda = gsl_vector_ptr(lambda,0);
  int status;
  int i, j, k;

  if(ldm*1.0e+18<0){mq[0] =-ldm*1.0e+18;}else{mq[0] =ldm*1.0e+18;}
  if(ldm*1.0e+18<0){mq[1] =-ldm*1.0e+18 + sdm*1.0e+18;}else{mq[1] =ldm*1.0e+18 + sdm*1.0e+18;}
  if(ldm*1.0e+18<0){mq[2] =(-ldm + ldm)*1.0e+18;}else{mq[2] =(ldm + ldm)*1.0e+18;}





  
    /* Calculate neutrino Hamiltonian :always calculate the correct matrix due to new parameters*/

    gsl_matrix_complex *Hp=gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);  
    
    glb_hamiltonian_cd_pedro(E,V*1.0e9, cp_sign,Hp);

    double complex (*_Hp)[GLB_NU_FLAVOURS]= (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(Hp, 0, 0);      


    /* Calculate eigenvalues of Hamiltonian */
    if ((status=zheevh3(_Hp, _Qp, _lambda)) != 0)
      return status;

     gsl_matrix_complex_free(Hp);



  //L*=1.97327e-10;

  /* Calculate S-Matrix in mass basis in matter ... */
  double phase;
  gsl_matrix_complex_set_zero(Sp);
  for (i=0; i < GLB_NU_FLAVOURS; i++)
  {
    phase    = -1e-9*L * _lambda[i];
    _Sp[i][i] = cos(phase) + I*sin(phase);
  }

  /* ... and transform it to the flavour basis */
  gsl_matrix_complex_set_zero(T0p);
  double complex *p = &_T0p[0][0];
  for (i=0; i < GLB_NU_FLAVOURS; i++)              /* T0 = S.Q^\dagger */
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      for (int k=0; k < GLB_NU_FLAVOURS; k++)
      {
        *p += ( creal(_Sp[i][k])*creal(_Qp[j][k])+cimag(_Sp[i][k])*cimag(_Qp[j][k]) )
                + I * ( cimag(_Sp[i][k])*creal(_Qp[j][k])-creal(_Sp[i][k])*cimag(_Qp[j][k]) );
      }
      p++;
    }

  gsl_matrix_complex_set_zero(Sp);
  p = &_Sp[0][0];
  for (i=0; i < GLB_NU_FLAVOURS; i++)              /* S = Q.T0 */
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      for (k=0; k < GLB_NU_FLAVOURS; k++)
      {
        *p += ( creal(_Qp[i][k])*creal(_T0p[k][j])-cimag(_Qp[i][k])*cimag(_T0p[k][j]) )
                + I * ( cimag(_Qp[i][k])*creal(_T0p[k][j])+creal(_Qp[i][k])*cimag(_T0p[k][j]) );
      }
      p++;
    }


  gsl_matrix_complex_free(T0p);
  gsl_matrix_complex_free(Qp);
  gsl_vector_free(lambda);


  return 0;
}


/***************************************************************************
 * Function glb_probability_matrix                                         *
 ***************************************************************************
 * Calculates the neutrino oscillation probability matrix.                 *
 ***************************************************************************
 * Parameters:                                                             *
 *   P:       Buffer for the storage of the matrix                         *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 *   E:       Neutrino energy (in GeV)                                     *
 *   psteps:  Number of layers in the matter density profile               *
 *   length:  Lengths of the layers in the matter density profile in km    *
 *   density: The matter densities in g/cm^3                               *
 *   filter_sigma: Width of low-pass filter or <0 for no filter            *
 *   user_data: Unused here, should be NULL                                *
 ***************************************************************************/
int my_probability_matrix(double P[3][3], int cp_sign, double E,
    int psteps, const double *length, const double *density,
    double filter_sigma, void *user_data)
{
  int status;
  int i, j,l,k;

    gsl_matrix_complex *Sp=gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS); 
    gsl_matrix_complex_set_zero(Sp);
   
   
    if (psteps > 1)
    {

      gsl_matrix_complex *S1p=gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS); 
      gsl_matrix_complex *T0p=gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS); 
      gsl_matrix_complex_set_identity(S1p);                                 /* S1 = 1 */
      gsl_matrix_complex_set_zero(T0p);
      for (i=0; i < psteps; i++)
      {
        glb_S_matrix_cd_pedro(E, GLB_KM_TO_EV(length[i]), density[i]*GLB_V_FACTOR*GLB_Ne_MANTLE, cp_sign,Sp); 

        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, Sp, S1p, GSL_COMPLEX_ZERO, T0p);   /* T0 = S.S1 */
        gsl_matrix_complex_memcpy(S1p, T0p);                                 /* S1 = T0 */
      }
      gsl_matrix_complex_memcpy(Sp, S1p);                                  /* S = S1 */
      gsl_matrix_complex_free(S1p);
      gsl_matrix_complex_free(T0p);
    }
    else
    {
     glb_S_matrix_cd_pedro(E, GLB_KM_TO_EV(length[0]), density[0]*GLB_V_FACTOR*GLB_Ne_MANTLE, cp_sign,Sp); 
    }
          
    double complex (*_Sp)[GLB_NU_FLAVOURS]= (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(Sp, 0, 0); 


    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++){
	   P[j][i] = SQR_ABS(_Sp[i][j]);
        }

  gsl_matrix_complex_free(Sp);

  return 0;
}



