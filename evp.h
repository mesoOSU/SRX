#ifndef EVP_H
#define EVP_H
#include <stdlib.h>
#include <mpi.h>
#include <fftw_mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <float.h>
#include <vector>
#include <algorithm>
#include "evpFlags.h"
#include "V3math.h"
#include "pf_xrx.h"

/**********************
  types and structures
  *********************/
typedef double real;
#define MPI_real MPI_DOUBLE
//typedef float real;
//#define MPI_real MPI_FLOAT
typedef struct{
	real x;
	real y;
	real z;
}Vec3R;
struct G_Info{
	int ID;
	real t1;
	real Phi;
	real t2;
	static bool before(const G_Info& lx, const G_Info& rx){return lx.ID<rx.ID;}
};
struct N_Info{
	int X;
	int Y;
	int Z;
	int ID;
//	static bool before(const G_Info& lx, const G_Info& rx){return lx.ID<rx.ID;}
};
typedef real ten2nd[3][3];
typedef real ten4th[3][3][3][3];
typedef real voigt[6];
typedef real voigt5[5];
typedef real voigtch[NSYSMX][5];
typedef real voigt66[6][6];


/***************************
  loops, allocators, indices
 ***************************/
#define sIDX (((px+lxs)*CellDim[1] + py)*CellDim[2] + pz)	/* index for space array */
#define pIDX ((px*CellDim[1] + py)*CellDim[2] + pz)	/* index for local array */
#define AllocMem(a, n, t) a = (t*)malloc((n)*sizeof(t))
#define B(i,j,p) (B[p-1][i-1][j-1])
#define Np_loop for(ip=0;ip<N_phases;ip++)
#define local_loop for(px=0;px<lnx;px++)for(py=0;py<CellDim[1];py++)for(pz=0;pz<CellDim[2];pz++)
#define T2_loop for(mi=0;mi<3;mi++)for(mj=0;mj<3;mj++)
#define T2p_loop for(mip=0;mip<3;mip++)for(mjp=0;mjp<3;mjp++)
#define C6_loop for(mi=0;mi<6;mi++)for(mj=0;mj<6;mj++)
#define C4_loop for(mi=0;mi<3;mi++)for(mj=0;mj<3;mj++)for(mk=0;mk<3;mk++)for(ml=0;ml<3;ml++)
#define PError(s,i) perror(#s),\
	MPI_Abort(MPI_COMM_WORLD, i)
extern int ip;		// dummy indices for global usage
extern int px,py,pz;
extern int mi, mj, mk, ml, mip, mjp;	// loop over the 3x3/3x3x3x3 tensors (stress/strain, etc)
#ifdef EVP_GLOBALS_ONCE
int ip;
int px,py,pz;
int mi, mj, mk, ml, mip, mjp;
#endif


/***************************
  units & constants
 ***************************/
#define PI 3.14159265359
#define TWO_PI 6.28318530718
#define TRIAL_FREQ 1.0E10
#define TRIAL_FREQ_DRX 1.0E10
#define EV_IN_J 1.60217733e-19
#define J_IN_EV (1.0/EV_IN_J)
#define MY_J_IN_EV 0.62415063631
#define BOLZ_IN_J__K 1.380658e-23
#define BOLZ BOLZ_IN_J__K
#define kT (T__K*BOLZ*J_IN_EV)
#define EPSILON_FFT 1E-5


/***************************
  MPI parameters
 ***************************/
extern int mpirank, NumPE;
extern int lnx, lxs, lnyt, lyst, lsize,lnx_pf,lxs_pf,lnyt_pf,lyst_pf,lsize_pf1;
extern MPI_Status status;
#ifdef EVP_GLOBALS_ONCE
int mpirank = 0, NumPE = 0;
int lnx = 0, lxs = 0, lnyt = 0, lyst = 0, lsize = 0 , lnx_pf=0, lxs_pf = 0 , lnyt_pf =0 , lyst_pf = 0, lsize_pf1=0;
MPI_Status status;
#endif


/***************************
  fft, k-space
 ***************************/
extern fftwnd_mpi_plan plan, iplan, plan_pf;
extern fftw_complex *fft_data, *fft_work;
extern Vec3R *g;
extern real *g2;
extern ten4th *GAMMA;	// Green operator in Fourier space, a function of Ref. stiffness and frequency
extern int ElastBC;	// boundar condition used when calculating GAMMA. Need to be consistent with Stress_BC_Flag
extern  ten2nd *kSig_r, *kSig_i;
#ifdef EVP_GLOBALS_ONCE
fftwnd_mpi_plan plan, iplan, plan_pf;
fftw_complex *fft_data = 0, *fft_work = 0;
Vec3R *g = 0;
real *g2 = 0;
ten4th *GAMMA = 0;
int ElastBC;
ten2nd *kSig_r=0, *kSig_i=0;
#endif


/***************************
 Computational grid, ICs, BCs 
 Simulation controllers
 ***************************/
extern int CellDim[3];
extern int Nxyz;
extern int old_id;
extern real TimeStep;
extern real TimeTot;
extern int N_steps;
extern real Err;
extern int IterMax;
extern int PrintControl[2];
extern int Update_Flag;
extern int Hard_Flag;
extern int Tex_Flag;
extern int VelGrad_BC_Flag[6];// pyz: we currently assume symmetric vel. gradient applied
extern real VelGrad_BC[6];	
extern ten2nd Udot;		// full matrix form of VelGrad_BC
extern ten2nd Udot_s, Udot_a;	// sym and ant part of Udot
extern voigt D_bar6;	// voigt notation of Udot_s
extern voigt5 D_bar5;	// devi. part of D_bar6
extern ten2nd DisGradAvg;	// macro displ. gradent
extern ten2nd dDisGradAvg;	// equivalent to "ddisgradmacro"
extern ten2nd dDisGradAvg_acum;	// equivalent to "ddisgradmacroacum"
extern int Stress_BC_Flag[6];
extern int CREEP_FLAG;
extern real Stress_BC[6];
extern ten2nd Scauchy;	// full matrix form of Stress_BC
#ifdef EVP_GLOBALS_ONCE
int CellDim[3] = {0};
int Nxyz = 0;
int old_id = 0;
real TimeStep = 0.0;
real TimeTot = 0.0;
int N_steps = 0;
real Err = 0.0;
int IterMax = 0;
int PrintControl[2] = {0};
int Update_Flag = 0;
int Hard_Flag = 0;
int Tex_Flag = 0;
int VelGrad_BC_Flag[6] = {0};
real VelGrad_BC[6] = {0.};
ten2nd Udot = {0.0};
ten2nd Udot_s = {0.};
ten2nd Udot_a = {0.};
voigt D_bar6 = {0.};
voigt5 D_bar5 = {0.};
ten2nd DisGradAvg = {0.};
ten2nd dDisGradAvg = {0.};
ten2nd dDisGradAvg_acum = {0.0};
int Stress_BC_Flag[6] = {0};
int CREEP_FLAG = 0;
real Stress_BC[6] = {0.};
ten2nd Scauchy = {0.};
#endif


/***************************
	measures and errors
 ***************************/
extern real d_vm;	//	von Mises (vm) strain rate
extern real s_vm;	//  von Mises stress
extern real e_vm;	//  von Mises strain
extern real s_vm1;	// von Mises stress in phase-1
extern real Err_e;	// error in strain
extern real Err_s;	// error in stress
#ifdef EVP_GLOBALS_ONCE
real d_vm = 0.;
real s_vm = 0.;
real e_vm = 0.;
real s_vm1 = 0.;
real Err_e = 0.;
real Err_s = 0.;
#endif


/***************************
  global files
 ***************************/
extern FILE *fp_vm;
extern FILE *fp_dd;
extern FILE *fp_err;
#ifdef EVP_GLOBALS_ONCE
FILE *fp_vm = 0;
FILE *fp_dd = 0;
FILE *fp_err = 0;
#endif

/**********************
  Material types and phases
  *********************/
/*
   N_phases -- # of phases under consideration
   Slip_Phase* -- Path of the input file containing plasticity (slip) of phase *
   */
extern int N_phases;	// <= 2 for current version
extern int Type_phases[2];
extern real ElastConst_PhaseI[4];
extern real ElastConst_PhaseII[4];
extern char Slip_PhaseI[100];	// currently assume only 2 phases at most
extern char Slip_PhaseII[100];
extern char initial_ms[100]; // file name of inital microstructure
extern ten4th Cijkl[NPHMAX];
extern real **nSRS;		// strain rate sensitvity of each slip/twin system for each phase
extern int *nSYS;		// number of slip+twin systems for each phase
extern int N_modes_max, N_modes;
extern char XtalSys[10];
extern real XtalAxis[3];
extern int *iMode;
//extern real ***bn, ***bb;
extern real bn[NPHMAX][NSYSMX][3];
extern real bb[NPHMAX][NSYSMX][3];
extern real bp[NPHMAX][NSYSMX][3];
extern voigt5 Schm_xt[NPHMAX][NSYSMX];	// Schmid tensors in xtal axes.<-> # of slip systems
extern voigtch *Schm_gr;	// Schmid tensors in sample axes <-> # of grid points
extern ten2nd B[6];		// used in chg_basis
#ifdef EVP_GLOBALS_ONCE
int N_phases = 0;
int Type_phases[2] = {0};
real ElastConst_PhaseI[4] = {0.0};
real ElastConst_PhaseII[4] = {0.0};
char Slip_PhaseI[100] = {0};
char Slip_PhaseII[100] = {0};
char initial_ms[100] = {0};
ten4th Cijkl[NPHMAX] = {0.0};
real **nSRS = 0;
int *nSYS = 0;
int N_modes_max, N_modes;
char XtalSys[10] = {0};
real XtalAxis[3] = {0.0};
int *iMode = 0;
//real ***bn, ***bb;
real bn[NPHMAX][NSYSMX][3];
real bb[NPHMAX][NSYSMX][3];
real bp[NPHMAX][NSYSMX][3];
voigt5 Schm_xt[NPHMAX][NSYSMX] = {0.0};
voigtch *Schm_gr;
ten2nd B[6] = {{0.0}};		// used in chg_basis
#endif


/**********************
	Grains
  *********************/
extern real WGT;	// weight of each grid point, i.e. 1/Nxyz
extern real Wgt_ph1;	// weight of phase#1 (solid)
extern ten4th C0;	// avg Cijkl of grain essemble
extern ten4th S0;	// avg Sijkl of grain essemble
extern voigt66 C066;
extern voigt66 caux66;// 6x6 form of C0
extern voigt66 S066;	// 6x6 form of S0
extern voigt66 *C_gr;	// Field containing the inhomogeneous Cijkl of grain essemble
extern voigt66 *FS_gr;	// Field used in DRX_ELASTIC_TREATlocal_
extern int *grain_f;
extern int *grain_f_new;	// Field to identify grain #
extern std::vector<G_Info> gID_list;// a list of grain IDs and euler angles
//extern std::vector<N_Info> nucleiID_list; //list of nuclei coordinates
extern int *phase_f;	// Field to identify phase #
extern ten2nd *TranMat_xt2sa;	// Field to store the transformation matrices (xtal -> sample)
extern int init_gid;
#ifdef EVP_GLOBALS_ONCE
real WGT = 0.0;
real Wgt_ph1 = 0.0;
ten4th C0 = {0.0};
ten4th S0 = {0.0};
voigt66 C066 = {0.};
voigt66 caux66 = {0.};
voigt66 S066 = {0.0};
voigt66 *C_gr = 0;
voigt66 *FS_gr = 0;	// Field used in DRX_ELASTIC_TREAT
int *grain_f = 0;
int *grain_f_new = 0;
std::vector<G_Info> gID_list;
int *phase_f = 0;
ten2nd *TranMat_xt2sa = 0;
int init_gid=0;
#endif


/**********************
  hardening
  *********************/
extern int nsmx, isectwx;
extern real nrsx, gamd0x, twshx, tau0xf, tau0xb, tau1x, thet0, thet1, hselfx, hlatex;
extern real ***trial_tau;	// trial crss
extern real ***crss;
extern real **xkin;
extern real *gamacum;	// accumulated shear
extern real tau[NPHMAX][NSYSMX][3];
extern real thet[NPHMAX][NSYSMX][2];
extern real gam0[NPHMAX][NSYSMX];
extern real Hard[NPHMAX][NSYSMX][NSYSMX];	// strain dependent hardening rates
#ifdef EVP_GLOBALS_ONCE
int nsmx, isectwx;
real nrsx, gamd0x, twshx, tau0xf, tau0xb, tau1x, thet0, thet1, hselfx, hlatex;
real ***trial_tau = 0;	// trial crss
real ***crss = 0;
real **xkin = 0;
real *gamacum=0;	// accumulated shear
real tau[NPHMAX][NSYSMX][3] = {0.0};
real thet[NPHMAX][NSYSMX][2] = {0.0};
real gam0[NPHMAX][NSYSMX] = {0.0};
real Hard[NPHMAX][NSYSMX][NSYSMX] = {0.0};
#endif


/**********************
  dislocation related
  *********************/
#ifdef DD_BASED_FLAG
extern real **rho_m;	// mobile dislocation density
extern real **rho_s;	// SSD density
extern real **rho_dipole;
extern real **rho_forest;
extern real **rho_inplane;
extern real **trial_rho_s;	// trial SSD density
extern real **rho_g1, **rho_g2, **rho_g3;	// GND density
extern real **rho_dot_s, **rho_dot_g1, **rho_dot_g2;  // dislocation evolution rates
extern Vec3R **gadot_grad, **ga_grad;	// gradient of shear rate and shear
extern real **trial_rho_g1, **trial_rho_g2, **trial_rho_g3;	// trial GND density
extern real *local_pxp_send, *local_pxm_send;	// data sent to neighboring PEs
extern real *local_pxp_recv, *local_pxm_recv;	// data (shear and shear rate) received from neighboring PEs
extern int *local_pxpG_send, *local_pxmG_send;	// data (grain IDs) sent to neighboring PEs
extern int *local_pxpG_recv, *local_pxmG_recv;	// data (grain IDs) received from neighboring PEs
extern int **GB_checkX, **GB_checkY, **GB_checkZ;	// check if the current point is at GB (3 directions): [0] is "+" direction, [1] is "-" direction
extern real **gam_acum;	// accumulated shear field on each slip system
extern real **trial_gam_acum;	// trial version of gam_acum
extern real **rho_F;	// Forest dislocation density
extern real **rho_P;	// Parallel dislocation density
extern real rho_F_avg, rho_P_avg, rho_s_avg, rho_m_avg, rho_g_avg, rho_dipole_avg, rho_forest_avg, rho_inplane_avg;
extern real bb_B[NPHMAX][NSYSMX];	// Burgers vectors
extern real a0[NPHMAX];		// lattice parameters
extern real L0;	// voxel length used to determine GND (and rate)
extern real GND_ScaleRatio;	// Scaling ratio used to scale L0
extern real C_1[NPHMAX],C_2[NPHMAX],C_3[NPHMAX],C_4[NPHMAX],C_5[NPHMAX],C_6[NPHMAX],C_7[NPHMAX],C_8[NPHMAX];
extern real Q_slip[NPHMAX], Q_bulk[NPHMAX];
extern real Selfinter[NPHMAX],Coplanar[NPHMAX],CrossSlip[NPHMAX],GlissileJunction[NPHMAX],HirthLock[NPHMAX],LomerCottrellLock[NPHMAX];
extern real  Selfinter0[NPHMAX],Coplanar0[NPHMAX],CrossSlip0[NPHMAX],GlissileJunction0[NPHMAX],HirthLock0[NPHMAX],LomerCottrellLock0[NPHMAX];
extern real T__K;	// temperature
extern real cos_theta_s[NPHMAX][NSYSMX][NSYSMX], sin_theta_s[NPHMAX][NSYSMX][NSYSMX];	// angle for SSD projection
extern real cos_theta_g1[NPHMAX][NSYSMX][NSYSMX], sin_theta_g1[NPHMAX][NSYSMX][NSYSMX];	// angle for GND 1st type
extern real cos_theta_g2[NPHMAX][NSYSMX][NSYSMX], sin_theta_g2[NPHMAX][NSYSMX][NSYSMX];	// angle for GND 2nd type
extern real cos_theta_g3[NPHMAX][NSYSMX][NSYSMX], sin_theta_g3[NPHMAX][NSYSMX][NSYSMX];	// angle for GND 3rd type
extern real **Chi[NPHMAX];
extern real **Chi_1[NPHMAX];// interaction strength between slip systems
extern real rho_SSD_initial;	// initial value of SSD density
extern real Shear_G[NPHMAX];
extern real *rho_tot;	// used in either record different dislocations or DRX
extern real D_0;        //diffusion coefficient in the climb equation
#ifdef EVP_GLOBALS_ONCE
real **rho_m = 0;	
real **rho_s = 0;
real **rho_dipole = 0;
real **rho_forest = 0;
real **rho_inplane = 0;
real **trial_rho_s = 0;
real **rho_g1 = 0, **rho_g2 = 0, **rho_g3 = 0;
real **rho_dot_s = 0, **rho_dot_g1 = 0, **rho_dot_g2 = 0;  // dislocation evolution rates
real **trial_rho_g1 = 0, **trial_rho_g2 = 0, **trial_rho_g3 = 0;
Vec3R **gadot_grad = 0, **ga_grad = 0;	
real *local_pxp_send = 0, *local_pxm_send = 0;	// data sent to neighboring PEs
real *local_pxp_recv = 0, *local_pxm_recv = 0;	// data received from neighboring PEs
int *local_pxpG_send = 0, *local_pxmG_send = 0;	// data received from neighboring PEs
int *local_pxpG_recv = 0, *local_pxmG_recv = 0;	// data received from neighboring PEs
int **GB_checkX = 0, **GB_checkY = 0, **GB_checkZ = 0;
real **gam_acum = 0;
real **trial_gam_acum = 0;	// updated in StrainRate_Orowan() using guess stress
real **rho_F = 0;
real **rho_P = 0;
real rho_F_avg = 0.0, rho_P_avg = 0.0, rho_s_avg = 0.0, rho_m_avg = 0.0, rho_g_avg = 0.0, rho_dipole_avg = 0.0, rho_forest_avg = 0.0, rho_inplane_avg = 0.0;
real bb_B[NPHMAX][NSYSMX] = {0.0};
real a0[NPHMAX] = {0.0};
real L0 = 0.0;
real GND_ScaleRatio = 0.0;
real C_1[NPHMAX],C_2[NPHMAX],C_3[NPHMAX],C_4[NPHMAX],C_5[NPHMAX],C_6[NPHMAX],C_7[NPHMAX],C_8[NPHMAX];
real Q_slip[NPHMAX], Q_bulk[NPHMAX];
real Selfinter[NPHMAX],Coplanar[NPHMAX],CrossSlip[NPHMAX],GlissileJunction[NPHMAX],HirthLock[NPHMAX],LomerCottrellLock[NPHMAX];
real Selfinter0[NPHMAX],Coplanar0[NPHMAX],CrossSlip0[NPHMAX],GlissileJunction0[NPHMAX],HirthLock0[NPHMAX],LomerCottrellLock0[NPHMAX];
real T__K = 0;	
real cos_theta_s[NPHMAX][NSYSMX][NSYSMX]={0.0}, sin_theta_s[NPHMAX][NSYSMX][NSYSMX]={0.0};	// angle for SSD projection
real cos_theta_g1[NPHMAX][NSYSMX][NSYSMX]={0.0}, sin_theta_g1[NPHMAX][NSYSMX][NSYSMX]={0.0};	// angle for GND 1st type
real cos_theta_g2[NPHMAX][NSYSMX][NSYSMX]={0.0}, sin_theta_g2[NPHMAX][NSYSMX][NSYSMX]={0.0};	// angle for GND 2nd type
real cos_theta_g3[NPHMAX][NSYSMX][NSYSMX]={0.0}, sin_theta_g3[NPHMAX][NSYSMX][NSYSMX]={0.0};	// angle for GND 3rd type
real **Chi[NPHMAX]={0};
real **Chi_1[NPHMAX] = {0};// interaction strength between slip systems
real rho_SSD_initial = 0.0;
real Shear_G[NPHMAX] = {0.0};
real *rho_tot=0;
real D_0=0.0;
#endif
#endif

/**********************
  stress/strain fields
  macro stress/strain
  *********************/
extern ten2nd *Sig;	
extern ten2nd *fluctuation;// stress in Eq. (2). After iteration, it should be the stress at time t+dt.
extern ten2nd *Eps;		// plastic strain field in Eq. (2). At beginning of iteration, it is the known quantity for time t.
extern ten2nd *Edot;	// plastic strain rate in Eq. (2). Note this is in correspondance to the current stress (Sig).
extern ten2nd SigAvg, EpsAvg, EdotAvg;
extern ten2nd SigDevAvg, VelGradAvg;
extern ten2nd SigAvg1;
extern ten2nd *DisGrad;		// displacement gradient field used during iteration. It gives the total strain in Eq. (15)
extern ten2nd *VelGrad;		// velocity gradient field
extern real **gamdot;	// shear rate fields
extern real **new_position;
extern real **displacement_fluct;
#ifdef EVP_GLOBALS_ONCE
ten2nd *Sig = 0;
ten2nd *fluctuation = 0;
ten2nd *Eps = 0;
ten2nd *Edot = 0;
ten2nd SigAvg = {0.};
ten2nd SigAvg1 = {0.};
ten2nd SigDevAvg = {0.0};
ten2nd EpsAvg = {0.};
ten2nd EdotAvg = {0.0};
ten2nd VelGradAvg = {0.};
ten2nd *DisGrad = 0;
ten2nd *VelGrad = 0;
real **gamdot =0;	// shear rate fields
real **new_position = 0;
real **displacement_fluct = 0;
#endif


/**********************
  parameters for DRX
  *********************/
#ifdef DD_BASED_FLAG
#ifdef PF_DRX
extern int CellDim_pf[3];
extern int lx_pf, ly_pf, lz_pf, lsize_pf;
extern real pf_Vol, pf_dt;
extern real E_gb, M_gb;
//extern real M_bar;	// the normalized GB mobility in phase-field
extern real ssdScaler_DRX;
extern real gndScaler_DRX;
extern real sigScaler_DRX;
extern real Scale_Fdeform;
extern int nucleus_radius;
extern int flag_DRX;
extern real *diff_rho;
//extern real *dd_average;
extern real *dd_average_drx;
extern real *rho_tot_drx; // in case PF grid differs from FFT grid
extern int *grain_f_drx; // in case PF grid differs from FFT grid
//extern int *gID_rex_drx; // in case PF grid differs from FFT grid
extern int *gID_new_drx; // in case PF grid differs from FFT grid
//extern int *grex_new_drx; // in case PF grid differs from FFT grid
extern int *first_pf_drx;
extern real *table_rho_drx; // in case PF grid differs from FFT grid
extern real *table_dd_average_drx;
extern int *table_pf_drx;
extern int *table_grain_drx;// in case PF grid differs from FFT grid
extern int *GB_indicator;
extern int *GB_type;
extern int *gID_rex;
extern int *growth;
extern int *swap_check;
extern int *nuclei;
extern int *nucleation_count;
extern int *gID_new;
//extern int *grex_new;
extern int *first_pf;
extern real k_c, alpha, zeta_kc, age_drx;
extern real *k_c_field;
extern real L00;  // store the initial L0
extern real TimeStep_DRX; // time step in phase-field simulation
extern int Nucl_Static_Flag;
extern int DRX_Interpolation_Flag;
extern real *kappa_drx;
extern const gsl_rng_type *RandType;
extern gsl_rng *RandInstance;
extern int RandSeed;
extern FILE *fp_drx;
extern FILE *fp1;
extern ten2nd *Eps0_last;
extern ten2nd *Tau_0;
extern ten2nd *Tau;
extern real *t_last;	// the last time DRX nucleation occurs
extern real tau_DRX;	// characteristic aging time for DRX grains
extern real DeltaQ_DRX;	// additional barrier for slip for new DRX grains
extern real DeltaTau_DRX;	// additional passing stress for new DRX grains
extern real pf_length_scale;
extern real pf_time_step;
#ifdef EVP_GLOBALS_ONCE
int CellDim_pf[3]={0};
int lx_pf=0, ly_pf=0, lz_pf=0, lsize_pf=0;
real pf_Vol=0.0, pf_dt=0.0;
real E_gb=0.0;
real M_gb=0.0;
//real M_bar = 0.0;
real ssdScaler_DRX = 0.0;
real gndScaler_DRX = 0.0;
real sigScaler_DRX = 0.0;
real Scale_Fdeform=0.0;
int nucleus_radius=0;
int flag_DRX=0;
real *diff_rho=0;
real *rho_tot_drx=0;
//real *dd_average=0;
real *dd_average_drx =0;
int *grain_f_drx=0;
real *table_rho_drx=0;
int *table_grain_drx=0;
real *table_dd_average_drx=0;
int *table_pf_drx=0;
int *gID_rex=0;
int *growth =0;
int *swap_check =0;
int *nuclei = 0;
int *nucleation_count = 0;
//int *gID_rex_drx=0;
int *gID_new=0;
int *gID_new_drx=0;
//int *grex_new=0;
int *first_pf=0;
int *first_pf_drx=0;
//int *grex_new_drx=0;
int *GB_indicator=0;
int *GB_type=0;
real k_c=0.0, alpha=0.0, zeta_kc=0.0, age_drx=0.0;
real *k_c_field=0;
real L00=0.0;  // store the initial L0
real TimeStep_DRX=0.0; // time step in phase-field simulation
int Nucl_Static_Flag=0;
int DRX_Interpolation_Flag=0;
real *kappa_drx=0;
const gsl_rng_type *RandType=0;
gsl_rng *RandInstance=0;
int RandSeed=0;
FILE *fp_drx=0;
FILE *fp1 =0;
//FILE *fp_stat=0;
ten2nd *Eps0_last=0;
ten2nd *Tau_0=0;
ten2nd *Tau=0;
real *t_last=0;		// the last time DRX nucleation occurs
real tau_DRX=0.0;
real DeltaQ_DRX=0.0;
real DeltaTau_DRX=0.0;	// additional passing stress for new DRX grains
real pf_length_scale=0.0;
real pf_time_step=0.0;
#endif
#endif
#endif

/**********************
  io.c
  *********************/
int GetNameList(char**);
void PrintNameList(void);
void OpenFiles(void);
void CloseFiles(void);
void WriteEdotMPI(char *s, int step);
void WriteEpsMPI(char *s, int step);
void WriteElsMPI(char *s, int step);
void WriteDisgradMPI(char *s, int step);
void WriteSigMPI(char *s, int step);
void WriteTextureMPI(char *s, int step);
void WriteSLIPMPI(char *s, int step);
void WriteNewPositionMPI(char *s,int step);
void WriteDDMPI(char *s, int step);
void WriteGrainMPI(char *s, int step);
void WriteGrainPFMPI(char *s, int step);
#ifdef DD_BASED_FLAG
void WriteRhoMPI(char *s, char *type, int step);
void WriteRhoDotMPI(char *s, int step);
#endif

/**********************
  init.c
  *********************/
void SetupJob(void);
void DestroyJob(void);

/**********************
  kinematics.c
  *********************/
void VoigtToFull(real*, ten2nd);
void chg_basis(voigt V6, ten2nd T2, voigt66 C2, ten4th T4, int opt);
void chg_basis5(voigt5 V6, ten2nd T2, voigt66 C2, ten4th T4, int opt);
void SymAntDecompose(ten2nd t, ten2nd s, ten2nd a);
void Ten4thTransform(ten4th c1, ten2nd a, ten4th c2, int opt);
void EulerToTransMatrix(real *t1, real *ph, real *t2, ten2nd a, int opt);
void LU_inv_66(voigt66 c);
void update_schmid(void);
real VonMises(ten2nd);


/**********************
  evolution.c
  *********************/
void Evolution(void);


/**********************
  constitutive.c
  *********************/
void StrainRate_eval(voigt,voigt,voigt66,int,int);
void get_trialtau(int, int);
void update_orient(void);
void harden(void);
void get_trialtau_anal(int, int);
void harden_anal(void);
#ifdef DD_BASED_FLAG
void StrainRate_Orowan(voigt, voigt,voigt66,int,int);
#ifdef DD_POWER_LAW
void StrainRate_Orowan_POWER(voigt, voigt,voigt66,int,int);
#endif
void DislocationEvolution(int istep);
void trial_DislocationEvolution(int, int);
void AvgDislocation(void);

void Gradient_ExchangeGrainID(void);
#ifdef DD_GND
void Gradient_ShearRate(void);
void Gradient_Shear(void);
#endif
#endif


/**********************
  integrated.c
  *********************/
#ifdef DD_BASED_FLAG
#ifdef PF_DRX
int NucleationCheckDRX(int STEP, real* kc, real alpha);
void InterpolateGID(void);
void InterpolateGrid(void);
void Update_microstructure(void);
real DrivingForceBoundaryMigration(void);
#endif
#endif

/**********************
  elasticity.c
  *********************/
void HomogeneousStressSolver(ten4th c0ijkl, ten2nd *eps0, ten2nd epsavg,
		ten2nd *eps, ten2nd *sig);
void InhomogeneousStressSolver(voigt66 *cij, ten2nd *eps0, ten2nd epsavg,
		ten4th c0ijkl, ten2nd *eps0_last,
		ten2nd *Tau_0, ten2nd *Tau,
		ten2nd *eps, ten2nd *sig);
void ElasticEVP(ten2nd epsavg);
void ElasticEVP_2(ten2nd epsavg, ten2nd *esp);

/**********************
  testf.c
  *********************/
void PrintTensor(ten2nd);
void PrintTensorNorm(ten2nd);
void PrintB(void);
void PrintVoigt(voigt);
void WriteEtaMPI(char *s);	// record microstructure
void PrintVoigt66(voigt66 m);
void PrintSchmidXt(void);
void Ti_alpha_beta(void);
void SingleXtal_Al(void);
void Bicrystal_Ti(int nx, int ny, int nz);
void Bicrystal_Ti_diffuse(int nx, int ny, int nz);
void Bicrystal_Ti_2(int nx, int ny, int nz);
void Bicrystal_Ti_3(int nx, int ny, int nz);
void Gamma_GammaPrime(int nx, int ny, int nz, real WithFrac);
void SphericalInclusion(void);
void SphericalVoid(void);

#endif	// end evp.h

