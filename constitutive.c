 #include "evp.h"
#include<set>


void StrainRate_eval(voigt stress, voigt edot, voigt66 d_edot, int idx, int jph)
{
	int i,j,k;
	int isign;
	real rss[NSYSMX];	// resolved shear stress
	real rss1[NSYSMX];
	real rss2[NSYSMX];
	voigt5 sc[NSYSMX];
	real tau_t[NSYSMX][2];	// crss
	real nsr[NSYSMX];	// strain rate sensitivity
	real xkinaux[NSYSMX];

	for(i=0;i<nSYS[jph-1];i++){
		nsr[i] = nSRS[jph-1][i];	// SRS could be different for different slip systems
		tau_t[i][0] = trial_tau[idx][i][0];
		tau_t[i][1] = trial_tau[idx][i][1];
		xkinaux[i] = xkin[idx][i];
		for(j=0;j<5;j++) sc[i][j] = Schm_gr[idx][i][j];
	}

	// calculate resolved shear stress and shear rate
	for(i=0;i<nSYS[jph-1];i++){
		rss[i] = sc[i][0]*stress[0]+sc[i][1]*stress[1]+sc[i][2]*stress[2]+
			sc[i][3]*stress[3]+sc[i][4]*stress[4];
		if((rss[i]-xkinaux[i])<0){
			isign = 1;
		}
		else{
			isign = 0;
		}
		rss[i] = (rss[i]-xkinaux[i])/tau_t[i][isign];
		rss1[i] = gam0[jph-1][i]*nsr[i]*fabs(pow(rss[i],nsr[i]-1.))/tau_t[i][isign];
		rss2[i] = gam0[jph-1][i]*fabs(pow(rss[i],nsr[i]))*((real)(2*(rss[i]>0)-1));
	
		// shear rate for each system
		gamdot[idx][i] = rss2[i];
	}

	// plastic strain rate
	for(i=0;i<5;i++){
		edot[i] = 0.;
		for(j=0;j<nSYS[jph-1];j++){
			edot[i] += sc[j][i]*rss2[j];
		}
	}
	edot[5] = 0.0;

	for(i=0;i<5;i++){
		for(j=0;j<5;j++){
			d_edot[i][j] = 0.0;
			for(k=0;k<nSYS[jph-1];k++){
				d_edot[i][j] += sc[k][i]*sc[k][j]*rss1[k];
			}
		}
	}

	for(i=0;i<6;i++){
		d_edot[i][5] = 0.0;
		d_edot[5][i] = 0.0;
	}
	
	return;
}/*end StrainRate_eval()*/

void get_trialtau(int idx, int jph)
{
	int i,j;
	real gamtot, deltgam;
	real dtau;
	real tau0, tau1, thet0, thet1;
	real voce, fact, exp_ini, exp_del;
	real TINY;

	gamtot = gamacum[idx];
	deltgam = 0.0;
	for(i=0;i<nSYS[jph-1];i++)
		deltgam += fabs(gamdot[idx][i])*TimeStep;

	for(i=0;i<nSYS[jph-1];i++){
		dtau = 0.0;
		for(j=0;j<nSYS[jph-1];j++)
			dtau += Hard[jph-1][i][j]*fabs(gamdot[idx][j])*TimeStep;
		tau0 = tau[jph-1][i][0];
		tau1 = tau[jph-1][i][2];
		thet0 = thet[jph-1][i][0];
		thet1 = thet[jph-1][i][1];
		TINY = 1.E-4*tau0;

		voce = 0.0;
		if(fabs(thet0)>TINY){
			voce = thet1*deltgam;
			if(fabs(tau1)>TINY){
				fact = fabs(thet0/tau1);
				exp_ini = exp(-1.0*gamtot*fact);
				exp_del = exp(-1.0*deltgam*fact);
				voce += -1.*(fact*tau1-thet1)/fact*exp_ini*(exp_del-1.0)-
					thet1/fact*exp_ini*(exp_del*((gamtot+deltgam)*fact+1.0)-(gamtot*fact+1.0));
			}
		}

		trial_tau[idx][i][0] = crss[idx][i][0]+dtau*voce/deltgam;
		trial_tau[idx][i][1] = crss[idx][i][1]+dtau*voce/deltgam;
	}

	return;
}/*end get_trialtau()*/

void get_trialtau_anal(int idx, int jph)
{
	/* use analytical hardening law */
	int i;
	real gamtot, deltgam;
	real dtau;
	real Koeff, N_exponent;

	Koeff = 2.0;
	N_exponent = 0.3;

	gamtot = gamacum[idx];
	deltgam = 0.0;
	for(i=0;i<nSYS[jph-1];i++)
		deltgam += fabs(gamdot[idx][i])*TimeStep;

//	dtau = Koeff*N_exponent*pow(gamtot+deltgam, N_exponent-1.0)*deltgam;
	//dtau = gamtot+deltgam;
	dtau = deltgam;
	for(i=0;i<nSYS[jph-1];i++){

		//trial_tau[idx][i][0] = crss[idx][i][0]+dtau;
		//trial_tau[idx][i][1] = crss[idx][i][1]-dtau;
		trial_tau[idx][i][0] = crss[idx][i][0]+Koeff*dtau;
		trial_tau[idx][i][1] = crss[idx][i][1]+Koeff*dtau;
	}

	return;
}/*end get_trialtau_anal()*/

static void Orientation(ten2nd tran, ten2nd rot)
{
	int i;
	real w[3];
	real w_norm, w_tan;
	ten2nd omega,omega2;
	ten2nd tran_new;

	w[0] = rot[2][1];
	w[1] = rot[0][2];
	w[2] = rot[1][0];

	w_norm = 0.0;
	for(i=0;i<3;i++)
		w_norm += w[i]*w[i];
	w_norm = sqrt(w_norm);
	w_tan = tan(w_norm/2.0);
	if(fabs(w_norm)<1E-6)
		w_norm = 1.0;
	for(i=0;i<3;i++)
		w[i] *= w_tan/w_norm;

	// construct "normalized" omega matrix
	omega[0][0]=0.; omega[0][1]=-1.*w[2]; omega[0][2]=w[1];
	omega[1][0]=w[2]; omega[1][1]=0.; omega[1][2]=-1.*w[0];
	omega[2][0]=-1.*w[1]; omega[2][1]=w[0]; omega[2][2]=0.;
	T2_loop{
		omega2[mi][mj] = 0.0;
		for(i=0;i<3;i++){
			omega2[mi][mj] += omega[mi][i]*omega[i][mj];
		}
	}

	w_tan = w_tan*w_tan;
	T2_loop{
		rot[mi][mj] = (real)(mi==mj) +
			2.*(omega[mi][mj]+omega2[mi][mj])/(1.+w_tan);
	}

	// transf. matrix at t+dt
	T2_loop{
		tran_new[mi][mj] = 0.0;
		for(i=0;i<3;i++)
			tran_new[mi][mj] += rot[mi][i]*tran[i][mj];
	}

	// record new transf. matrix
	T2_loop{
		tran[mi][mj] = tran_new[mi][mj];
	}

	return;
}/*end Orientation()*/

void update_orient(void)
{
	int is, i,j;
	int jph;
	ten2nd tranmat;	// xtal -> sample (grain)
	ten2nd RotLoc;
	ten2nd Lp, RotSlip;
	ten2nd Rot;	// total rotation
	real bb_sa[3], bn_sa[3];
	real rsl_bar, rlc_bar;
	real rsl_bar_tot, rlc_bar_tot;


	rsl_bar = 0.0;
	rlc_bar = 0.0;
	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			// local rotate rate
			T2_loop{
				RotLoc[mi][mj] = (VelGrad[pIDX][mi][mj]-VelGrad[pIDX][mj][mi])/2.;
				Lp[mi][mj] = 0.0;
				tranmat[mi][mj] = TranMat_xt2sa[pIDX][mi][mj];
			}

			// slip ration rate
			for(is=0;is<nSYS[jph-1];is++){
				for(i=0;i<3;i++){
					bb_sa[i] = 0.0;
					bn_sa[i] = 0.0;
					for(j=0;j<3;j++){
						bb_sa[i] += tranmat[i][j]*bb[jph-1][is][j];
						bn_sa[i] += tranmat[i][j]*bn[jph-1][is][j];
					}
				}
				T2_loop{
					Lp[mi][mj] += bb_sa[mi]*bn_sa[mj]*gamdot[pIDX][is];
				}
			}
			T2_loop{
				RotSlip[mi][mj] = (Lp[mi][mj]-Lp[mj][mi])/2.;
			}

			// avg rotation rate
			rsl_bar += sqrt(RotSlip[2][1]*RotSlip[2][1] +
					RotSlip[0][2]*RotSlip[0][2] + RotSlip[1][0]*RotSlip[1][0])*WGT;
			rlc_bar += sqrt(RotLoc[2][1]*RotLoc[2][1] +
					RotLoc[0][2]*RotLoc[0][2] + RotLoc[1][0]*RotLoc[1][0])*WGT;

			// total rotation
			T2_loop{
				// Udot_a is the applied rotation (anti-symm of Udot)
				// RotSlip does NOT change crystal axes orientation
				Rot[mi][mj] = (Udot_a[mi][mj]+RotLoc[mi][mj]-RotSlip[mi][mj])*TimeStep;
			}

			// reorientation
			Orientation(tranmat,Rot);

			// update Transformation matrix
			T2_loop{
				TranMat_xt2sa[pIDX][mi][mj] = tranmat[mi][mj];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&rsl_bar, &rsl_bar_tot, 1, MPI_real,
			MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rlc_bar, &rlc_bar_tot, 1, MPI_real,
			MPI_SUM, MPI_COMM_WORLD);

	if(mpirank==0){
		printf("Average plastic rotation = %e\n",rsl_bar_tot);
		printf("Average local rotation = %e\n",rlc_bar_tot);
	}

	return;
}/*end update_orient()*/

void harden(void)
{
	int i,j;
	int jph;
	real gamtot, deltgam;
	real dtau;
	real tau0, tau1, thet0, thet1;
	real voce, fact, exp_ini, exp_del;
	real TINY;

	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			gamtot = gamacum[pIDX];
			deltgam = 0.0;
			for(i=0;i<nSYS[jph-1];i++)
				deltgam += fabs(gamdot[pIDX][i])*TimeStep;

			for(i=0;i<nSYS[jph-1];i++){
				dtau = 0.;
				for(j=0;j<nSYS[jph-1];j++){
					dtau += Hard[jph-1][i][j]*fabs(gamdot[pIDX][j])*TimeStep;
				}
				tau0 = tau[jph-1][i][0];
				tau1 = tau[jph-1][i][2];
				thet0 = thet[jph-1][i][0];
				thet1 = thet[jph-1][i][1];
				TINY = 1.E-4*tau0;

				voce = 0.0;
				if(fabs(thet0)>TINY){
					voce = thet1*deltgam;
					if(fabs(tau1)>TINY){
						fact = fabs(thet0/tau1);
						exp_ini = exp(-1.0*gamtot*fact);
						exp_del = exp(-1.0*deltgam*fact);
						voce += -1.*(fact*tau1-thet1)/fact*exp_ini*(exp_del-1.0)-
							thet1/fact*exp_ini*(exp_del*((gamtot+deltgam)*fact+1.0)-(gamtot*fact+1.0));
					}
				}
				crss[pIDX][i][0] += dtau*voce/deltgam;
				crss[pIDX][i][1] += dtau*voce/deltgam;

				trial_tau[pIDX][i][0] = crss[pIDX][i][0];
				trial_tau[pIDX][i][1] = crss[pIDX][i][1];
			}
			gamacum[pIDX] = gamtot+deltgam;
		}
	}

	return;
}/*end harden()*/

void harden_anal(void)
{
	int i;
	int jph;
	real gamtot, deltgam;
	real dtau;
	real Koeff, N_exponent;

	Koeff = 2.0;
	N_exponent = 0.3;

	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			gamtot = gamacum[pIDX];
			deltgam = 0.0;
			for(i=0;i<nSYS[jph-1];i++)
				deltgam += fabs(gamdot[pIDX][i])*TimeStep;

			//dtau = Koeff*N_exponent*pow(gamtot+deltgam, N_exponent-1.0)*deltgam;
			//dtau = gamtot+deltgam;
			dtau = deltgam;
			for(i=0;i<nSYS[jph-1];i++){
				//crss[pIDX][i][0] += dtau;
				//crss[pIDX][i][1] -= dtau;
				crss[pIDX][i][0] += Koeff*dtau;
				crss[pIDX][i][1] += Koeff*dtau;

				trial_tau[pIDX][i][0] = crss[pIDX][i][0];
				trial_tau[pIDX][i][1] = crss[pIDX][i][1];
			}
			gamacum[pIDX] = gamtot+deltgam;
		}
	}


	return;
}/*end harden_anal()*/

#ifdef DD_BASED_FLAG

#ifdef DD_POWER_LAW
void StrainRate_Orowan_POWER(voigt stress, voigt edot, voigt66 d_edot, int idx, int jph)
{
	/* Using the constitutive laws for strain rate, derived through Orowan equation
	   and the Arrenius equation of thermal activated slip.
input:
	stress -- the Cauchy stress tensor for the given grid point
	idx -- the index of the grid point
	jph -- the phase identity of the grid point
output:
	edot -- the plastic strain rate for the given grid point
	d_edot -- the derivative of edot w.r.t. Cauchy stress, used in augmented Lagrangians */

	int i,j,k;
	real rss[NSYSMX];	// resolved shear stress
	real rss1[NSYSMX];	// shear rates
	real rss2[NSYSMX];	
	real nsr[NSYSMX];	// strain rate sensitivity
	voigt5 sc[NSYSMX];
	real rhoS[NSYSMX];
	real rhoG1[NSYSMX];
	real rhoG2[NSYSMX];
	real rhoG3[NSYSMX];
	real tau_pass;	// athermal resistance or passing stress for mobile dislocations
	real tau_cut;	// cutting stress due to forest dislocaions
	real tau_eff;	// effective shear stress
	real gam_o;		// pre-factor, reference shear rate
	real tmp_gam;
	real shear_G;	// shear modulus
	real c_1,c_2,c_3;
	real q_slip;

	shear_G = Shear_G[jph-1];	// It has been checked that the phase is xtal.
	c_1 = C_1[jph-1];
	c_2 = C_2[jph-1];
	c_3 = C_3[jph-1];
	q_slip = Q_slip[jph-1]*MY_J_IN_EV/kT;

	for(i=0;i<nSYS[jph-1];i++){
		rhoS[i] = trial_rho_s[idx][i];
		rhoG1[i] = trial_rho_g1[idx][i];
		rhoG2[i] = trial_rho_g2[idx][i];
		rhoG3[i] = trial_rho_g3[idx][i];
		if(rhoS[i]<0||rhoG1[i]<0||rhoG2[i]<0||rhoG3[i]<0){
			printf("pyz:rank%d:idx=%d,rhos=%le, tmp_rhos=%le\n",mpirank,idx,rho_s[idx][i],trial_rho_s[idx][i]);
		}
	}
	/* Calcuate parallel and forest dislocation densities */
	for(i=0;i<nSYS[jph-1];i++){
		rho_F[idx][i] = 0.0;
		rho_P[idx][i] = 0.0;
		for(j=0;j<nSYS[jph-1];j++){
			rho_F[idx][i] += Chi[jph-1][i][j]*(rhoS[j]*cos_theta_s[jph-1][i][j]+
					fabs(rhoG1[j]*cos_theta_g1[jph-1][i][j])+
					fabs(rhoG2[j]*cos_theta_g2[jph-1][i][j])+
					fabs(rhoG3[j]*cos_theta_g3[jph-1][i][j]));
			rho_P[idx][i] += Chi[jph-1][i][j]*(rhoS[j]*sin_theta_s[jph-1][i][j]+
					fabs(rhoG1[j]*sin_theta_g1[jph-1][i][j])+
					fabs(rhoG2[j]*sin_theta_g2[jph-1][i][j])+
					fabs(rhoG3[j]*sin_theta_g3[jph-1][i][j]));
		}
	}


	for(i=0;i<nSYS[jph-1];i++){
		nsr[i] = nSRS[jph-1][i];	// SRS could be different for different slip systems
		for(j=0;j<5;j++) sc[i][j] = Schm_gr[idx][i][j];
	}

	for(i=0;i<nSYS[jph-1];i++){
		//tau_pass = c_1*shear_G*bb_B[jph-1][i]*sqrt(rho_P[idx][i])*1.0E-3;	// MPa
		tau_pass = c_1*shear_G*bb_B[jph-1][i]*sqrt(fabs(rho_P[idx][i]))*1.0E-3;	// MPa
		/* Consider both forward and backward jumps */
		tau_cut = (q_slip*kT/MY_J_IN_EV)*sqrt(rho_F[idx][i])/c_2/c_3/bb_B[jph-1][i]/bb_B[jph-1][i]*0.1;	// MPa
    //
		// resolved shear stress
		rss[i] = sc[i][0]*stress[0]+sc[i][1]*stress[1]+sc[i][2]*stress[2]+
			sc[i][3]*stress[3]+sc[i][4]*stress[4];
		if(fabs(rss[i])>tau_pass){
			tau_eff = fabs(fabs(rss[i])-tau_pass);
		}
		else{
			tau_eff = 0.0;
		}
    tau_eff /=tau_cut;
    tau_eff = MIN(tau_eff,1.2);
//    if(grex_new[idx]==1){
//      printf("pyz:rank=%d:tau_eff=%e\n",mpirank,tau_eff);
//    }
    rss1[i] = gam0[jph-1][i]*nsr[i]*fabs(pow(tau_eff,nsr[i]-1.))/tau_cut;
    rss2[i] = gam0[jph-1][i]*fabs(pow(tau_eff,nsr[i]))*((real)(2*(rss[i]>0)-1));

		if(isnan(rss1[i])||isnan(rss2[i])){
			printf("pyz:rss is nan rank%d:idx=%d,tau_eff=%le, tau_cut=%le, rss1=%le, rss2=%le\n",mpirank,idx,tau_eff,tau_cut,
          rss1[i],rss2[i]);
		}
    if(rss1[i]>DBL_MAX){
			printf("pyz:rank%d:idx=%d,tau_eff=%le, tau_cut=%le, tau_pass=%le, rss=%le, rss1=%le, rss2=%le\n",
          mpirank,idx,tau_eff,tau_cut,tau_pass,rss[i],
          rss1[i],rss2[i]);
		}
    
		// shear rate for each system
		gamdot[idx][i] = rss2[i];
    
		// accumulated shear for each system
		trial_gam_acum[idx][i] = gam_acum[idx][i] + rss2[i]*TimeStep;
	}

	// plastic strain rate
	for(i=0;i<5;i++){
		edot[i] = 0.;
		for(j=0;j<nSYS[jph-1];j++){
			edot[i] += sc[j][i]*rss2[j];
		}
	}
	edot[5] = 0.0;

	for(i=0;i<5;i++){
		for(j=0;j<5;j++){
			d_edot[i][j] = 0.0;
			for(k=0;k<nSYS[jph-1];k++){
				d_edot[i][j] += sc[k][i]*sc[k][j]*rss1[k];
			}
		}
	}

	for(i=0;i<6;i++){
		d_edot[i][5] = 0.0;
		d_edot[5][i] = 0.0;
	}


  return;
}/*end StrainRate_Orowan_POWER()*/
#endif


void StrainRate_Orowan(voigt stress, voigt edot, voigt66 d_edot, int idx, int jph)
{
	/* Using the constitutive laws for strain rate, derived through Orowan equation
	   and the Arrenius equation of thermal activated slip.
input:
	stress -- the Cauchy stress tensor for the given grid point
	idx -- the index of the grid point
	jph -- the phase identity of the grid point
output:
	edot -- the plastic strain rate for the given grid point
	d_edot -- the derivative of edot w.r.t. Cauchy stress, used in augmented Lagrangians */

	int i,j,k;
	real rss[NSYSMX];	// resolved shear stress
	real gam[NSYSMX];	// shear rates
	real gam2[NSYSMX];	
	voigt5 sc[NSYSMX];
	real rhoS[NSYSMX];
	real rhoG1[NSYSMX];
	real rhoG2[NSYSMX];
	real rhoG3[NSYSMX];
	real rho_f[NSYSMX];
	real tau_pass;	// athermal resistance or passing stress for mobile dislocations
	real tau_cut;	// cutting stress due to forest dislocaions
	real tau_eff;	// effective shear stress
	real gam_o;		// pre-factor, reference shear rate
	real tmp_gam;
	real shear_G;	// shear modulus
	real c_1,c_2,c_3;
	real q_slip;
    real tau_0 = tau[jph-1][i][0];
	shear_G = Shear_G[jph-1];	// It has been checked that the phase is xtal.
	c_1 = C_1[jph-1];
	c_2 = C_2[jph-1];
	c_3 = C_3[jph-1];
#ifdef PF_DRX
	q_slip = Q_slip[jph-1]*MY_J_IN_EV/kT;
 // if(grex_new[idx]==1) q_slip *= 1.0;
#else
	q_slip = Q_slip[jph-1]*MY_J_IN_EV/kT;
#endif

	for(i=0;i<nSYS[jph-1];i++){
		rhoS[i] = trial_rho_s[idx][i];
		rhoG1[i] = trial_rho_g1[idx][i];
		rhoG2[i] = trial_rho_g2[idx][i];
		rhoG3[i] = trial_rho_g3[idx][i];
		if(isnan(rhoS[i])){
			printf("pyz:rank%d:idx=%d,rhos=%le, tmp_rhos=%le\n",mpirank,idx,rho_s[idx][i],trial_rho_s[idx][i]);
		}
	}
	/* Calcuate parallel and forest dislocation densities */
	for(i=0;i<nSYS[jph-1];i++){
		rho_F[idx][i] = 0.0;
		rho_P[idx][i] = 0.0;
		rho_f[i] =0.0;
		
		for(j=0;j<nSYS[jph-1];j++){
			rho_F[idx][i] += Chi_1[jph-1][i][j]*(rhoS[j]+
					fabs(rhoG1[j]*cos_theta_g1[jph-1][i][j])+
					fabs(rhoG2[j]*cos_theta_g2[jph-1][i][j])+
					fabs(rhoG3[j]*cos_theta_g3[jph-1][i][j]));
				rho_f[i] += Chi_1[jph-1][i][j]*(rhoS[j]+
					fabs(rhoG1[j]*cos_theta_g1[jph-1][i][j])+
					fabs(rhoG2[j]*cos_theta_g2[jph-1][i][j])+
					fabs(rhoG3[j]*cos_theta_g3[jph-1][i][j]));		
			rho_P[idx][i] += Chi[jph-1][i][j]*(rhoS[j]+
					fabs(rhoG1[j]*sin_theta_g1[jph-1][i][j])+
					fabs(rhoG2[j]*sin_theta_g2[jph-1][i][j])+
					fabs(rhoG3[j]*sin_theta_g3[jph-1][i][j]));
		}
	}


	for(i=0;i<nSYS[jph-1];i++){
		for(j=0;j<5;j++) sc[i][j] = Schm_gr[idx][i][j];
	}

	for(i=0;i<nSYS[jph-1];i++){
		//tau_pass = c_1*shear_G*bb_B[jph-1][i]*sqrt(rho_P[idx][i])*1.0E-3;	// MPa
		tau_pass = c_1*shear_G*bb_B[jph-1][i]*sqrt(fabs(rho_P[idx][i] + rho_m[idx][i]))*1.0E-3;	// MPa
#ifdef PF_DRX
    tau_pass += DeltaTau_DRX*exp((t_last[idx]-TimeTot)/TimeStep/tau_DRX);
#endif

		// resolved shear stress
		rss[i] = sc[i][0]*stress[0]+sc[i][1]*stress[1]+sc[i][2]*stress[2]+
			sc[i][3]*stress[3]+sc[i][4]*stress[4];
		if(fabs(rss[i])>tau_pass){
			tau_eff = (fabs(rss[i])-tau_pass);
		}
	/*	else{
			tau_eff = 0.0;
		} */

#ifdef FORWARD_ONLY
		/* Only forward jump is considered */
		tau_cut = (q_slip*kT/MY_J_IN_EV)*sqrt(rho_F[idx][i])/c_2/c_3/bb_B[jph-1][i]/bb_B[jph-1][i]*0.1;	// MPa
#ifdef PF_DRX
  //  if(grex_new[idx]){
		  gam_o = kT*TRIAL_FREQ_DRX*sqrt(rho_P[idx][i])/c_1/c_3/shear_G/bb_B[jph-1][i]/bb_B[jph-1][i]*0.160217733;	// 1/s
   // }else{
	//	  gam_o = kT*TRIAL_FREQ*sqrt(rho_P[idx][i])/c_1/c_3/shear_G/bb_B[jph-1][i]/bb_B[jph-1][i]*0.160217733;	// 1/s
   // }
#else
		gam_o = kT*TRIAL_FREQ*sqrt(rho_P[idx][i])/c_1/c_3/shear_G/bb_B[jph-1][i]/bb_B[jph-1][i]*0.160217733;	// 1/s
#endif
		rho_m[idx][i] = 2.0*kT/c_1/c_2/c_3/shear_G/(pow(bb_B[jph-1][i],3.0))*sqrt(rho_P[idx][i]*rho_F[idx][i])*160.2177;
//    if((tau_eff/tau_cut)>200.0){
//      tau_eff = 1.0*tau_cut;
//    }
		tmp_gam = gam_o*exp(-1.0*q_slip*(1.0-tau_eff/tau_cut));		// 1/s
		if(isnan(tmp_gam)){
			printf("pyz:rank#%d:rho_F[%d][%d]=%le\n",mpirank,idx,i,rho_F[idx][i]);
			printf("pyz:rank#%d:tau_eff=%le, tau_cut=%le, tau_eff/tau_cut=%le\n",mpirank,tau_eff, tau_cut,tau_eff/tau_cut);
			fflush(stdout);
		}
		assert(isnan(tmp_gam)==0);
		gam[i] = tmp_gam*((real)(2*(rss[i]>0)-1));		// 1/s
		tmp_gam *= q_slip/tau_cut;
		gam2[i] = tmp_gam*(fabs(rss[i])>tau_pass);	// be careful of the derivative of absolute value	
#else
		/* Consider both forward and backward jumps */
		tau_cut = (kT/MY_J_IN_EV)*sqrt(rho_f[i])/c_2/c_3/bb_B[jph-1][i]/bb_B[jph-1][i]*0.1;	// MPa
		//if(grex_new[idx]){
      rho_m[idx][i] = 2.0*kT/c_1/c_2/c_3/shear_G/(pow(bb_B[jph-1][i],3.0))*sqrt(rho_P[idx][i]*rho_f[i])*160.2177;
		 gam_o = 1E-3*(rho_m[idx][i])*TRIAL_FREQ_DRX*c_2*bb_B[jph-1][i]*exp(-1.0*q_slip)/sqrt(rho_f[i]);
   // if(growth[pIDX] == 1 && (tau_eff/tau_cut)>1.0){
//printf("gam_0 = %le\n",gam_o);
     // tau_eff = 1.0*tau_cut;
  //  }
		
		if(fabs(rss[i])>tau_pass){
		tmp_gam = gam_o*sinh(tau_eff/(tau_0+tau_cut));
		} else {
		    tmp_gam = 0.0;
		}

		if(isnan(tmp_gam)){
			printf("pyz:rank#%d:rho_F[%d][%d]=%le\n",mpirank,idx,i,rho_F[idx][i]);
			for(j=0;j<nSYS[jph-1];j++){
				printf("pyz:rank#%d:Chi[%d][%d][%d]=%le\n",mpirank,jph-1,i,j,Chi[jph-1][i][j]);
				printf("pyz:rank#%d:cos_theta_s[%d][%d][%d]=%le\n",mpirank,jph-1,i,j,cos_theta_s[jph-1][i][j]);
				printf("pyz:rank#%d:rhoS[%d]=%le\n",mpirank,j,rhoS[j]);
				rho_F[idx][i] += Chi[jph-1][i][j]*(rhoS[j]*cos_theta_s[jph-1][i][j]+
						fabs(rhoG1[j]*cos_theta_g1[jph-1][i][j])+
						fabs(rhoG2[j]*cos_theta_g2[jph-1][i][j])+
						fabs(rhoG3[j]*cos_theta_g3[jph-1][i][j]));
			}

			printf("pyz:rank#%d:tau_eff=%le, tau_cut=%le, tau_eff/tau_cut=%le\n",mpirank,tau_eff, tau_cut,tau_eff/tau_cut);
			fflush(stdout);
		}
		assert(isnan(tmp_gam)==0);

		gam[i] = tmp_gam*((real)(2*(rss[i]>0)-1));
		if(fabs(rss[i])>tau_pass){ 
		tmp_gam = gam_o*cosh(tau_eff/(tau_0+tau_cut))/(tau_0+tau_cut);
		} else {
		    tmp_gam = 0.0; 
		}
		if(isnan(tmp_gam)){
			printf("pyz:rank#%d:rho_F[%d][%d]=%le\n",mpirank,idx,i,rho_F[idx][i]);
			for(j=0;j<nSYS[jph-1];j++){
				printf("pyz:rank#%d:Chi[%d][%d][%d]=%le\n",mpirank,jph-1,i,j,Chi[jph-1][i][j]);
				printf("pyz:rank#%d:cos_theta_s[%d][%d][%d]=%le\n",mpirank,jph-1,i,j,cos_theta_s[jph-1][i][j]);
				printf("pyz:rank#%d:rhoS[%d]=%le\n",mpirank,j,rhoS[j]);
				rho_F[idx][i] += Chi[jph-1][i][j]*(rhoS[j]*cos_theta_s[jph-1][i][j]+
						fabs(rhoG1[j]*cos_theta_g1[jph-1][i][j])+
						fabs(rhoG2[j]*cos_theta_g2[jph-1][i][j])+
						fabs(rhoG3[j]*cos_theta_g3[jph-1][i][j]));
			}

			printf("pyz:rank#%d:tau_eff=%le, tau_cut=%le, tau_eff/tau_cut=%le\n",mpirank,tau_eff, tau_cut,tau_eff/tau_cut);
			fflush(stdout);
		}
		assert(isnan(tmp_gam)==0);
		gam2[i] = tmp_gam*(fabs(rss[i])>tau_pass);	// be careful of the derivative of absolute value	
#endif
		

		// shear rate for each system
		gamdot[idx][i] = gam[i];
		// accumulated shear for each system
		trial_gam_acum[idx][i] = gam_acum[idx][i] + gam[i]*TimeStep;
	}

	// plastic strain rate
	for(i=0;i<5;i++){
		edot[i] = 0.;
		for(j=0;j<nSYS[jph-1];j++){
			edot[i] += sc[j][i]*gam[j];
		}
		if(isnan(edot[i])){
//exit(0);
			for(j=0;j<nSYS[jph-1];j++){
				printf("pyz:rank#%d:sc[%d][%d]=%le, gam[%d]=%le\n",mpirank,j,i,sc[j][i],j,gam[j]);
				fflush(stdout);
			}
		}
		assert(isnan(edot[i])==0);
	}
	edot[5] = 0.0;
	
	// derivative of strain rate w.r.t. Cauchy stress
	for(i=0;i<5;i++){
		for(j=0;j<5;j++){
			d_edot[i][j] = 0.0;
			for(k=0;k<nSYS[jph-1];k++){
				d_edot[i][j] += sc[k][i]*sc[k][j]*gam2[k];
			}
			assert(isnan(d_edot[i][j])==0);
		}
	}

	for(i=0;i<6;i++){
		d_edot[i][5] = 0.0;
		d_edot[5][i] = 0.0;
	}

	return;
}/*end StrainRate_Orowan()*/


void Gradient_ExchangeGrainID(void)
{
	int px_p, px_m, py_p, py_m, pz_p, pz_m;
	int pIDX_xp, pIDX_xm, pIDX_yp, pIDX_ym, pIDX_zp, pIDX_zm;

	if(NumPE>1){	// multiple PEs, use MPI_Send() and MPI_Recv()
		for(px=0;px<1;px++){
			for(py=0;py<CellDim[1];py++){
				for(pz=0;pz<CellDim[2];pz++){
					local_pxpG_send[py*CellDim[2]+pz] = grain_f[pIDX];
				}
			}
		}
		for(px=lnx-1;px<lnx;px++){
			for(py=0;py<CellDim[1];py++){
				for(pz=0;pz<CellDim[2];pz++){
					local_pxmG_send[py*CellDim[2]+pz] = grain_f[pIDX];
				}
			}
		}
//		printf("allocation_end");
		if(mpirank==0){
			MPI_Sendrecv(local_pxmG_send,CellDim[1]*CellDim[2],MPI_INT,1,0,
					local_pxpG_recv,CellDim[1]*CellDim[2],MPI_INT,1,3,
					MPI_COMM_WORLD, &status);
			MPI_Sendrecv(local_pxpG_send,CellDim[1]*CellDim[2],MPI_INT,NumPE-1,1,
					local_pxmG_recv,CellDim[1]*CellDim[2],MPI_INT,NumPE-1,2*(NumPE-1),
					MPI_COMM_WORLD, &status);
						printf("send_0");
		}
	
		else if(mpirank==NumPE-1){
			MPI_Sendrecv(local_pxpG_send,CellDim[1]*CellDim[2],MPI_INT,mpirank-1,2*mpirank+1,
					local_pxmG_recv,CellDim[1]*CellDim[2],MPI_INT,mpirank-1,2*(mpirank-1),
					MPI_COMM_WORLD, &status);
			MPI_Sendrecv(local_pxmG_send,CellDim[1]*CellDim[2],MPI_INT,0,2*mpirank,
					local_pxpG_recv,CellDim[1]*CellDim[2],MPI_INT,0,1,
					MPI_COMM_WORLD, &status);
//						printf("send_N");
		}
		else{
			MPI_Sendrecv(local_pxpG_send,CellDim[1]*CellDim[2],MPI_INT,mpirank-1,2*mpirank+1,
					local_pxmG_recv,CellDim[1]*CellDim[2],MPI_INT,mpirank-1,2*(mpirank-1),
					MPI_COMM_WORLD, &status);
			MPI_Sendrecv(local_pxmG_send,CellDim[1]*CellDim[2],MPI_INT,mpirank+1,2*mpirank,
					local_pxpG_recv,CellDim[1]*CellDim[2],MPI_INT,mpirank+1,2*(mpirank+1)+1,
					MPI_COMM_WORLD, &status);
		}
	}
	else{	// Only one PE (serial run), no need for communication
		for(px=0;px<1;px++){
			for(py=0;py<CellDim[1];py++){
				for(pz=0;pz<CellDim[2];pz++){
					local_pxpG_recv[py*CellDim[2]+pz] = grain_f[pIDX];
				}
			}
		}
		for(px=lnx-1;px<lnx;px++){
			for(py=0;py<CellDim[1];py++){
				for(pz=0;pz<CellDim[2];pz++){
					local_pxmG_recv[py*CellDim[2]+pz] = grain_f[pIDX];
				}
			}
		}

	}

	for(px=0;px<lnx;px++){
		px_p = px+1; px_m = px-1;
		if(px_p>=lnx) px_p = 0;
		if(px_m<0) px_m = lnx-1;
		for(py=0;py<CellDim[1];py++){
			py_p = py+1; py_m = py-1;
			if(py_p>=CellDim[1]) py_p = 0;
			if(py_m<0) py_m = CellDim[1]-1;
			for(pz=0;pz<CellDim[2];pz++){
				pz_p = pz+1; pz_m = pz-1;
				if(pz_p>=CellDim[2]) pz_p = 0;
				if(pz_m<0) pz_m = CellDim[2]-1;
				pIDX_xp = (px_p*CellDim[1]+py)*CellDim[2]+pz;
				pIDX_xm = (px_m*CellDim[1]+py)*CellDim[2]+pz;
				pIDX_yp = (px*CellDim[1]+py_p)*CellDim[2]+pz;
				pIDX_ym = (px*CellDim[1]+py_m)*CellDim[2]+pz;
				pIDX_zp = (px*CellDim[1]+py)*CellDim[2]+pz_p;
				pIDX_zm = (px*CellDim[1]+py)*CellDim[2]+pz_m;

				if(px==0){
					GB_checkX[pIDX][0] = (grain_f[pIDX]==grain_f[pIDX_xp]);
					GB_checkX[pIDX][1] = (grain_f[pIDX]==local_pxmG_recv[py*CellDim[2]+pz]);
					GB_checkY[pIDX][0] = (grain_f[pIDX]==grain_f[pIDX_yp]);
					GB_checkY[pIDX][1] = (grain_f[pIDX]==grain_f[pIDX_ym]);
					GB_checkZ[pIDX][0] = (grain_f[pIDX]==grain_f[pIDX_zp]);
					GB_checkZ[pIDX][1] = (grain_f[pIDX]==grain_f[pIDX_ym]);

          /* 09/18/2015 --- PYZ */
          int nn_gList[] = {grain_f[pIDX],grain_f[pIDX_xp],local_pxmG_recv[py*CellDim[2]+pz],
            grain_f[pIDX_yp],grain_f[pIDX_ym],grain_f[pIDX_zp],grain_f[pIDX_zm]};
          std::set<int> uqList;
          for(int iloc=0; iloc<7; iloc++){
            std::set<int>::iterator tmpit = uqList.find(nn_gList[iloc]);
            if(tmpit==uqList.end()) uqList.insert(nn_gList[iloc]);
          }
          GB_type[pIDX] = uqList.size()-1;  // 0-bulk; 1-normal; 2-triple; 3-quadruple; ...

				}
				else if(px==lnx-1){
					GB_checkX[pIDX][0] = (grain_f[pIDX]==local_pxpG_recv[py*CellDim[2]+pz]);
					GB_checkX[pIDX][1] = (grain_f[pIDX]==grain_f[pIDX_xm]);
					GB_checkY[pIDX][0] = (grain_f[pIDX]==grain_f[pIDX_yp]);
					GB_checkY[pIDX][1] = (grain_f[pIDX]==grain_f[pIDX_ym]);
					GB_checkZ[pIDX][0] = (grain_f[pIDX]==grain_f[pIDX_zp]);
					GB_checkZ[pIDX][1] = (grain_f[pIDX]==grain_f[pIDX_ym]);

          /* 09/18/2015 --- PYZ */
          int nn_gList[] = {grain_f[pIDX],local_pxpG_recv[py*CellDim[2]+pz],grain_f[pIDX_xm],
            grain_f[pIDX_yp],grain_f[pIDX_ym],grain_f[pIDX_zp],grain_f[pIDX_zm]};
          std::set<int> uqList;
          for(int iloc=0; iloc<7; iloc++){
            std::set<int>::iterator tmpit = uqList.find(nn_gList[iloc]);
            if(tmpit==uqList.end()) uqList.insert(nn_gList[iloc]);
          }
          GB_type[pIDX] = uqList.size()-1;
				}
				else{
					GB_checkX[pIDX][0] = (grain_f[pIDX]==grain_f[pIDX_xp]);
					GB_checkX[pIDX][1] = (grain_f[pIDX]==grain_f[pIDX_xm]);
					GB_checkY[pIDX][0] = (grain_f[pIDX]==grain_f[pIDX_yp]);
					GB_checkY[pIDX][1] = (grain_f[pIDX]==grain_f[pIDX_ym]);
					GB_checkZ[pIDX][0] = (grain_f[pIDX]==grain_f[pIDX_zp]);
					GB_checkZ[pIDX][1] = (grain_f[pIDX]==grain_f[pIDX_ym]);

          /* 09/18/2015 --- PYZ */
          int nn_gList[] = {grain_f[pIDX],grain_f[pIDX_xp],grain_f[pIDX_xm],
            grain_f[pIDX_yp],grain_f[pIDX_ym],grain_f[pIDX_zp],grain_f[pIDX_zm]};
          std::set<int> uqList;
          for(int iloc=0; iloc<7; iloc++){
            std::set<int>::iterator tmpit = uqList.find(nn_gList[iloc]);
            if(tmpit==uqList.end()) uqList.insert(nn_gList[iloc]);
          }
          GB_type[pIDX] = uqList.size()-1;
				}

			}
		}
	}
				
//printf("exchange_end");
	return;

}/*end Gradient_ExchangeGrainID()*/
#ifdef DD_GND
void Gradient_ShearRate(void)
{
	int i;
	int px_p, px_m, py_p, py_m, pz_p, pz_m;
	int pIDX_xp, pIDX_xm, pIDX_yp, pIDX_ym, pIDX_zp, pIDX_zm;

	for(i=0;i<NSYSMX;i++){

		/* Becuase of the slab decomposition, calculating the derivatives
		   requires the data from neighboring PEs */
		if(NumPE>1){	// multiple PEs, use MPI_Send() and MPI_Recv()
			for(px=0;px<1;px++){
				for(py=0;py<CellDim[1];py++){
					for(pz=0;pz<CellDim[2];pz++){
						local_pxp_send[py*CellDim[2]+pz] = gamdot[pIDX][i];
					}
				}
			}
			for(px=lnx-1;px<lnx;px++){
				for(py=0;py<CellDim[1];py++){
					for(pz=0;pz<CellDim[2];pz++){
						local_pxm_send[py*CellDim[2]+pz] = gamdot[pIDX][i];
					}
				}
			}
			if(mpirank==0){
				MPI_Sendrecv(local_pxm_send,CellDim[1]*CellDim[2],MPI_real,1,0,
						local_pxp_recv,CellDim[1]*CellDim[2],MPI_real,1,3,
						MPI_COMM_WORLD, &status);
				MPI_Sendrecv(local_pxp_send,CellDim[1]*CellDim[2],MPI_real,NumPE-1,1,
						local_pxm_recv,CellDim[1]*CellDim[2],MPI_real,NumPE-1,2*(NumPE-1),
						MPI_COMM_WORLD, &status);
			}
			else if(mpirank==NumPE-1){
				MPI_Sendrecv(local_pxp_send,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*mpirank+1,
						local_pxm_recv,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*(mpirank-1),
						MPI_COMM_WORLD, &status);
				MPI_Sendrecv(local_pxm_send,CellDim[1]*CellDim[2],MPI_real,0,2*mpirank,
						local_pxp_recv,CellDim[1]*CellDim[2],MPI_real,0,1,
						MPI_COMM_WORLD, &status);
			}
			else{
				MPI_Sendrecv(local_pxp_send,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*mpirank+1,
						local_pxm_recv,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*(mpirank-1),
						MPI_COMM_WORLD, &status);
				MPI_Sendrecv(local_pxm_send,CellDim[1]*CellDim[2],MPI_real,mpirank+1,2*mpirank,
						local_pxp_recv,CellDim[1]*CellDim[2],MPI_real,mpirank+1,2*(mpirank+1)+1,
						MPI_COMM_WORLD, &status);
			}
		}
		else{	// Only one PE (serial run), no need for communication
			for(px=0;px<1;px++){
				for(py=0;py<CellDim[1];py++){
					for(pz=0;pz<CellDim[2];pz++){
						local_pxp_recv[py*CellDim[2]+pz] = gamdot[pIDX][i];
					}
				}
			}
			for(px=lnx-1;px<lnx;px++){
				for(py=0;py<CellDim[1];py++){
					for(pz=0;pz<CellDim[2];pz++){
						local_pxm_recv[py*CellDim[2]+pz] = gamdot[pIDX][i];
					}
				}
			}

		}

		for(px=0;px<lnx;px++){
			px_p = px+1; px_m = px-1;
			if(px_p>=lnx) px_p = 0;
			if(px_m<0) px_m = lnx-1;
			for(py=0;py<CellDim[1];py++){
				py_p = py+1; py_m = py-1;
				if(py_p>=CellDim[1]) py_p = 0;
				if(py_m<0) py_m = CellDim[1]-1;
				for(pz=0;pz<CellDim[2];pz++){
					pz_p = pz+1; pz_m = pz-1;
					if(pz_p>=CellDim[2]) pz_p = 0;
					if(pz_m<0) pz_m = CellDim[2]-1;
					pIDX_xp = (px_p*CellDim[1]+py)*CellDim[2]+pz;
					pIDX_xm = (px_m*CellDim[1]+py)*CellDim[2]+pz;
					pIDX_yp = (px*CellDim[1]+py_p)*CellDim[2]+pz;
					pIDX_ym = (px*CellDim[1]+py_m)*CellDim[2]+pz;
					pIDX_zp = (px*CellDim[1]+py)*CellDim[2]+pz_p;
					pIDX_zm = (px*CellDim[1]+py)*CellDim[2]+pz_m;

					if(px==0){
						/* Gradient along x-axis */
						if((GB_checkX[pIDX][0]==0)&&(GB_checkX[pIDX][1]==1)){
							gadot_grad[pIDX][i].x = (gamdot[pIDX][i]-local_pxm_recv[py*CellDim[2]+pz]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==0)){
							gadot_grad[pIDX][i].x = (gamdot[pIDX_xp][i]-gamdot[pIDX][i]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==1)){
							gadot_grad[pIDX][i].x = (gamdot[pIDX_xp][i]-local_pxm_recv[py*CellDim[2]+pz])/2.0;
						}
						else{
							// Warning: only one sampling point along X-axis for the current grain!
							gadot_grad[pIDX][i].x = 0.0;
						}
							
						/* Gradient along y-axis */
						if((GB_checkY[pIDX][0]==0)&&(GB_checkY[pIDX][1]==1)){
							gadot_grad[pIDX][i].y = (gamdot[pIDX][i]-gamdot[pIDX_ym][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==0)){
							gadot_grad[pIDX][i].y = (gamdot[pIDX_yp][i]-gamdot[pIDX][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==1)){
							gadot_grad[pIDX][i].y = (gamdot[pIDX_yp][i]-gamdot[pIDX_ym][i])/2.0;
						}
						else{
							// Warning: only one sampling point along Y-axis for the current grain!
							gadot_grad[pIDX][i].y = 0.0;
						}

						/* Gradient along z-axis */
						if((GB_checkZ[pIDX][0]==0)&&(GB_checkZ[pIDX][1]==1)){
							gadot_grad[pIDX][i].z = (gamdot[pIDX][i]-gamdot[pIDX_zm][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==0)){
							gadot_grad[pIDX][i].z = (gamdot[pIDX_zp][i]-gamdot[pIDX][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==1)){
							gadot_grad[pIDX][i].z = (gamdot[pIDX_zp][i]-gamdot[pIDX_zm][i])/2.0;
						}
						else{
							// Warning: on.z one sampling point along Z-axis for the current grain!
							gadot_grad[pIDX][i].z = 0.0;
						}

					}
					else if(px==lnx-1){
						/* Gradient along x-axis */
						if((GB_checkX[pIDX][0]==0)&&(GB_checkX[pIDX][1]==1)){
							gadot_grad[pIDX][i].x = (gamdot[pIDX][i]-gamdot[pIDX_xm][i]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==0)){
							gadot_grad[pIDX][i].x = (local_pxp_recv[py*CellDim[2]+pz]-gamdot[pIDX][i]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==1)){
							gadot_grad[pIDX][i].x = (local_pxp_recv[py*CellDim[2]+pz]-gamdot[pIDX_xm][i])/2.0;
						}
						else{
							// Warning: only one sampling point along X-axis for the current grain!
							gadot_grad[pIDX][i].x = 0.0;
						}
							
						/* Gradient along y-axis */
						if((GB_checkY[pIDX][0]==0)&&(GB_checkY[pIDX][1]==1)){
							gadot_grad[pIDX][i].y = (gamdot[pIDX][i]-gamdot[pIDX_ym][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==0)){
							gadot_grad[pIDX][i].y = (gamdot[pIDX_yp][i]-gamdot[pIDX][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==1)){
							gadot_grad[pIDX][i].y = (gamdot[pIDX_yp][i]-gamdot[pIDX_ym][i])/2.0;
						}
						else{
							// Warning: only one sampling point along Y-axis for the current grain!
							gadot_grad[pIDX][i].y = 0.0;
						}

						/* Gradient along z-axis */
						if((GB_checkZ[pIDX][0]==0)&&(GB_checkZ[pIDX][1]==1)){
							gadot_grad[pIDX][i].z = (gamdot[pIDX][i]-gamdot[pIDX_zm][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==0)){
							gadot_grad[pIDX][i].z = (gamdot[pIDX_zp][i]-gamdot[pIDX][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==1)){
							gadot_grad[pIDX][i].z = (gamdot[pIDX_zp][i]-gamdot[pIDX_zm][i])/2.0;
						}
						else{
							// Warning: on.z one sampling point along Z-axis for the current grain!
							gadot_grad[pIDX][i].z = 0.0;
						}
					}
					else{
						/* Gradient along x-axis */
						if((GB_checkX[pIDX][0]==0)&&(GB_checkX[pIDX][1]==1)){
							gadot_grad[pIDX][i].x = (gamdot[pIDX][i]-gamdot[pIDX_xm][i]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==0)){
							gadot_grad[pIDX][i].x = (gamdot[pIDX_xp][i]-gamdot[pIDX][i]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==1)){
							gadot_grad[pIDX][i].x = (gamdot[pIDX_xp][i]-gamdot[pIDX_xm][i])/2.0;
						}
						else{
							// Warning: only one sampling point along X-axis for the current grain!
							gadot_grad[pIDX][i].x = 0.0;
						}
							
						/* Gradient along y-axis */
						if((GB_checkY[pIDX][0]==0)&&(GB_checkY[pIDX][1]==1)){
							gadot_grad[pIDX][i].y = (gamdot[pIDX][i]-gamdot[pIDX_ym][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==0)){
							gadot_grad[pIDX][i].y = (gamdot[pIDX_yp][i]-gamdot[pIDX][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==1)){
							gadot_grad[pIDX][i].y = (gamdot[pIDX_yp][i]-gamdot[pIDX_ym][i])/2.0;
						}
						else{
							// Warning: only one sampling point along Y-axis for the current grain!
							gadot_grad[pIDX][i].y = 0.0;
						}

						/* Gradient along z-axis */
						if((GB_checkZ[pIDX][0]==0)&&(GB_checkZ[pIDX][1]==1)){
							gadot_grad[pIDX][i].z = (gamdot[pIDX][i]-gamdot[pIDX_zm][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==0)){
							gadot_grad[pIDX][i].z = (gamdot[pIDX_zp][i]-gamdot[pIDX][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==1)){
							gadot_grad[pIDX][i].z = (gamdot[pIDX_zp][i]-gamdot[pIDX_zm][i])/2.0;
						}
						else{
							// Warning: on.z one sampling point along Z-axis for the current grain!
							gadot_grad[pIDX][i].z = 0.0;
						}
					}
				}
			}
		}
	}


	return;
}/*end Gradient_ShearRate()*/

void Gradient_Shear(void)
{
	int i;
	int px_p, px_m, py_p, py_m, pz_p, pz_m;
	int pIDX_xp, pIDX_xm, pIDX_yp, pIDX_ym, pIDX_zp, pIDX_zm;

	for(i=0;i<NSYSMX;i++){

		/* Becuase of the slab decomposition, calculating the derivatives
		   requires the data from neighboring PEs */
		if(NumPE>1){	// multiple PEs, use MPI_Send() and MPI_Recv()
			for(px=0;px<1;px++){
				for(py=0;py<CellDim[1];py++){
					for(pz=0;pz<CellDim[2];pz++){
						local_pxp_send[py*CellDim[2]+pz] = trial_gam_acum[pIDX][i];
					}
				}
			}
			for(px=lnx-1;px<lnx;px++){
				for(py=0;py<CellDim[1];py++){
					for(pz=0;pz<CellDim[2];pz++){
						local_pxm_send[py*CellDim[2]+pz] = trial_gam_acum[pIDX][i];
					}
				}
			}
			if(mpirank==0){
				MPI_Sendrecv(local_pxm_send,CellDim[1]*CellDim[2],MPI_real,1,0,
						local_pxp_recv,CellDim[1]*CellDim[2],MPI_real,1,3,
						MPI_COMM_WORLD, &status);
				MPI_Sendrecv(local_pxp_send,CellDim[1]*CellDim[2],MPI_real,NumPE-1,1,
						local_pxm_recv,CellDim[1]*CellDim[2],MPI_real,NumPE-1,2*(NumPE-1),
						MPI_COMM_WORLD, &status);
			}
			else if(mpirank==NumPE-1){
				MPI_Sendrecv(local_pxp_send,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*mpirank+1,
						local_pxm_recv,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*(mpirank-1),
						MPI_COMM_WORLD, &status);
				MPI_Sendrecv(local_pxm_send,CellDim[1]*CellDim[2],MPI_real,0,2*mpirank,
						local_pxp_recv,CellDim[1]*CellDim[2],MPI_real,0,1,
						MPI_COMM_WORLD, &status);
			}
			else{
				MPI_Sendrecv(local_pxp_send,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*mpirank+1,
						local_pxm_recv,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*(mpirank-1),
						MPI_COMM_WORLD, &status);
				MPI_Sendrecv(local_pxm_send,CellDim[1]*CellDim[2],MPI_real,mpirank+1,2*mpirank,
						local_pxp_recv,CellDim[1]*CellDim[2],MPI_real,mpirank+1,2*(mpirank+1)+1,
						MPI_COMM_WORLD, &status);
			}
		}
		else{	// Only one PE (serial run), no need for communication
			for(px=0;px<1;px++){
				for(py=0;py<CellDim[1];py++){
					for(pz=0;pz<CellDim[2];pz++){
						local_pxp_recv[py*CellDim[2]+pz] = trial_gam_acum[pIDX][i];
					}
				}
			}
			for(px=lnx-1;px<lnx;px++){
				for(py=0;py<CellDim[1];py++){
					for(pz=0;pz<CellDim[2];pz++){
						local_pxm_recv[py*CellDim[2]+pz] = trial_gam_acum[pIDX][i];
					}
				}
			}

		}

		for(px=0;px<lnx;px++){
			px_p = px+1; px_m = px-1;
			if(px_p>=lnx) px_p = 0;
			if(px_m<0) px_m = lnx-1;
			for(py=0;py<CellDim[1];py++){
				py_p = py+1; py_m = py-1;
				if(py_p>=CellDim[1]) py_p = 0;
				if(py_m<0) py_m = CellDim[1]-1;
				for(pz=0;pz<CellDim[2];pz++){
					pz_p = pz+1; pz_m = pz-1;
					if(pz_p>=CellDim[2]) pz_p = 0;
					if(pz_m<0) pz_m = CellDim[2]-1;
					pIDX_xp = (px_p*CellDim[1]+py)*CellDim[2]+pz;
					pIDX_xm = (px_m*CellDim[1]+py)*CellDim[2]+pz;
					pIDX_yp = (px*CellDim[1]+py_p)*CellDim[2]+pz;
					pIDX_ym = (px*CellDim[1]+py_m)*CellDim[2]+pz;
					pIDX_zp = (px*CellDim[1]+py)*CellDim[2]+pz_p;
					pIDX_zm = (px*CellDim[1]+py)*CellDim[2]+pz_m;


					if(px==0){
						/* Gradient along x-axis */
						if((GB_checkX[pIDX][0]==0)&&(GB_checkX[pIDX][1]==1)){
							ga_grad[pIDX][i].x = (trial_gam_acum[pIDX][i]-local_pxm_recv[py*CellDim[2]+pz]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==0)){
							ga_grad[pIDX][i].x = (trial_gam_acum[pIDX_xp][i]-trial_gam_acum[pIDX][i]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==1)){
							ga_grad[pIDX][i].x = (trial_gam_acum[pIDX_xp][i]-local_pxm_recv[py*CellDim[2]+pz])/2.0;
						}
						else{
							// Warning: only one sampling point along X-axis for the current grain!
							ga_grad[pIDX][i].x = 0.0;
						}
							
						/* Gradient along y-axis */
						if((GB_checkY[pIDX][0]==0)&&(GB_checkY[pIDX][1]==1)){
							ga_grad[pIDX][i].y = (trial_gam_acum[pIDX][i]-trial_gam_acum[pIDX_ym][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==0)){
							ga_grad[pIDX][i].y = (trial_gam_acum[pIDX_yp][i]-trial_gam_acum[pIDX][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==1)){
							ga_grad[pIDX][i].y = (trial_gam_acum[pIDX_yp][i]-trial_gam_acum[pIDX_ym][i])/2.0;
						}
						else{
							// Warning: only one sampling point along Y-axis for the current grain!
							ga_grad[pIDX][i].y = 0.0;
						}

						/* Gradient along z-axis */
						if((GB_checkZ[pIDX][0]==0)&&(GB_checkZ[pIDX][1]==1)){
							ga_grad[pIDX][i].z = (trial_gam_acum[pIDX][i]-trial_gam_acum[pIDX_zm][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==0)){
							ga_grad[pIDX][i].z = (trial_gam_acum[pIDX_zp][i]-trial_gam_acum[pIDX][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==1)){
							ga_grad[pIDX][i].z = (trial_gam_acum[pIDX_zp][i]-trial_gam_acum[pIDX_zm][i])/2.0;
						}
						else{
							// Warning: on.z one sampling point along Z-axis for the current grain!
							ga_grad[pIDX][i].z = 0.0;
						}

					}
					else if(px==lnx-1){
						/* Gradient along x-axis */
						if((GB_checkX[pIDX][0]==0)&&(GB_checkX[pIDX][1]==1)){
							ga_grad[pIDX][i].x = (trial_gam_acum[pIDX][i]-trial_gam_acum[pIDX_xm][i]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==0)){
							ga_grad[pIDX][i].x = (local_pxp_recv[py*CellDim[2]+pz]-trial_gam_acum[pIDX][i]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==1)){
							ga_grad[pIDX][i].x = (local_pxp_recv[py*CellDim[2]+pz]-trial_gam_acum[pIDX_xm][i])/2.0;
						}
						else{
							// Warning: only one sampling point along X-axis for the current grain!
							ga_grad[pIDX][i].x = 0.0;
						}
							
						/* Gradient along y-axis */
						if((GB_checkY[pIDX][0]==0)&&(GB_checkY[pIDX][1]==1)){
							ga_grad[pIDX][i].y = (trial_gam_acum[pIDX][i]-trial_gam_acum[pIDX_ym][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==0)){
							ga_grad[pIDX][i].y = (trial_gam_acum[pIDX_yp][i]-trial_gam_acum[pIDX][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==1)){
							ga_grad[pIDX][i].y = (trial_gam_acum[pIDX_yp][i]-trial_gam_acum[pIDX_ym][i])/2.0;
						}
						else{
							// Warning: only one sampling point along Y-axis for the current grain!
							ga_grad[pIDX][i].y = 0.0;
						}

						/* Gradient along z-axis */
						if((GB_checkZ[pIDX][0]==0)&&(GB_checkZ[pIDX][1]==1)){
							ga_grad[pIDX][i].z = (trial_gam_acum[pIDX][i]-trial_gam_acum[pIDX_zm][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==0)){
							ga_grad[pIDX][i].z = (trial_gam_acum[pIDX_zp][i]-trial_gam_acum[pIDX][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==1)){
							ga_grad[pIDX][i].z = (trial_gam_acum[pIDX_zp][i]-trial_gam_acum[pIDX_zm][i])/2.0;
						}
						else{
							// Warning: on.z one sampling point along Z-axis for the current grain!
							ga_grad[pIDX][i].z = 0.0;
						}
					}
					else{
						/* Gradient along x-axis */
						if((GB_checkX[pIDX][0]==0)&&(GB_checkX[pIDX][1]==1)){
							ga_grad[pIDX][i].x = (trial_gam_acum[pIDX][i]-trial_gam_acum[pIDX_xm][i]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==0)){
							ga_grad[pIDX][i].x = (trial_gam_acum[pIDX_xp][i]-trial_gam_acum[pIDX][i]);
						}
						else if((GB_checkX[pIDX][0]==1)&&(GB_checkX[pIDX][1]==1)){
							ga_grad[pIDX][i].x = (trial_gam_acum[pIDX_xp][i]-trial_gam_acum[pIDX_xm][i])/2.0;
						}
						else{
							// Warning: only one sampling point along X-axis for the current grain!
							ga_grad[pIDX][i].x = 0.0;
						}
							
						/* Gradient along y-axis */
						if((GB_checkY[pIDX][0]==0)&&(GB_checkY[pIDX][1]==1)){
							ga_grad[pIDX][i].y = (trial_gam_acum[pIDX][i]-trial_gam_acum[pIDX_ym][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==0)){
							ga_grad[pIDX][i].y = (trial_gam_acum[pIDX_yp][i]-trial_gam_acum[pIDX][i]);
						}
						else if((GB_checkY[pIDX][0]==1)&&(GB_checkY[pIDX][1]==1)){
							ga_grad[pIDX][i].y = (trial_gam_acum[pIDX_yp][i]-trial_gam_acum[pIDX_ym][i])/2.0;
						}
						else{
							// Warning: only one sampling point along Y-axis for the current grain!
							ga_grad[pIDX][i].y = 0.0;
						}

						/* Gradient along z-axis */
						if((GB_checkZ[pIDX][0]==0)&&(GB_checkZ[pIDX][1]==1)){
							ga_grad[pIDX][i].z = (trial_gam_acum[pIDX][i]-trial_gam_acum[pIDX_zm][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==0)){
							ga_grad[pIDX][i].z = (trial_gam_acum[pIDX_zp][i]-trial_gam_acum[pIDX][i]);
						}
						else if((GB_checkZ[pIDX][0]==1)&&(GB_checkZ[pIDX][1]==1)){
							ga_grad[pIDX][i].z = (trial_gam_acum[pIDX_zp][i]-trial_gam_acum[pIDX_zm][i])/2.0;
						}
						else{
							// Warning: on.z one sampling point along Z-axis for the current grain!
							ga_grad[pIDX][i].z = 0.0;
						}
					}
				}
			}
		}
	}

	return;
}/*end Gradient_Shear()*/
#endif

void DislocationEvolution(int istep)
{
	int i,j;
	int jph;
	real rss[NSYSMX],hydro[NSYSMX];	// resolved shear stress
	voigt5 sc[NSYSMX];
	ten2nd sig_aux;
	voigt sig6;
	ten4th aux3333;
	voigt66 aux66;
	real dot_SSD;	// the change rate of SSDs
	real dot_GND;	// the change rate of GNDs
	real shear_G;	// shear modulus
	real Poisson_nu;	// Poisson's ratio
	real c_4, c_5, c_6, c_7, c_8;
	real q_bulk;
	real d_dipole;

	real stress_climb, strainrate_climb;

	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			c_4 = C_4[jph-1];
			c_5 = C_5[jph-1];
			c_6 = C_6[jph-1];
			c_7 = C_7[jph-1];
			c_8 = C_8[jph-1];
			q_bulk = Q_bulk[jph-1]*MY_J_IN_EV/kT;	// reduced
			shear_G = Shear_G[jph-1];	// It has been checked that the phase is xtal.
			Poisson_nu = ElastConst_PhaseI[3];

			/* resolve shear stress */
			T2_loop{
				sig_aux[mi][mj] = Sig[pIDX][mi][mj];
			}
			chg_basis(sig6,sig_aux,aux66,aux3333,2);

			for(i=0;i<nSYS[jph-1];i++){
				for(j=0;j<5;j++) sc[i][j] = Schm_gr[pIDX][i][j];
			}
			for(i=0;i<nSYS[jph-1];i++){
			rss[i] = sc[i][0]*sig6[0]+sc[i][1]*sig6[1]+sc[i][2]*sig6[2]+
			sc[i][3]*sig6[3]+sc[i][4]*sig6[4];
			//	hydro[i] = (sig6[0]+sig6[1]+sig6[2])/3.0;		
			}

			for(i=0;i<nSYS[jph-1];i++){
				/* SSD */
				stress_climb = VonMises(Sig[pIDX])*3.0/2.0;
				strainrate_climb = VonMises(VelGrad[pIDX]);
			//	d_dipole = 1000/sqrt(rho_s[pIDX][i]);	//in nm
				d_dipole = sqrt(3.0)*shear_G*bb_B[jph-1][i]/16.0/PI/(1-Poisson_nu)/(fabs(rss[i])+EPSILON);	//in nm
				hydro[i] = shear_G*bb_B[jph-1][i]/PI/(1-Poisson_nu)/d_dipole;
				if(istep<=2000000) {
				  //  rho_m[pIDX][i] += (fabs(gamdot[pIDX][i])*(c_4*sqrt(rho_F[pIDX][i])) - fabs(gamdot[pIDX][i])*(c_6*d_dipole*rho_m[pIDX][i]))*TimeStep;
				    rho_dipole[pIDX][i] += (1E-17*rho_tot[pIDX]*0.5*1E24*(d_dipole)*(D_0/bb_B[jph-1][i])*pow(rho_s[pIDX][i],2.0)*(exp(fabs(hydro[i])*bb_B[jph-1][i]*bb_B[jph-1][i]*bb_B[jph-1][i]*100/(1.38*T__K)) -1.0)*exp(-1.0*q_bulk))*TimeStep;
				    rho_forest[pIDX][i] += fabs(gamdot[pIDX][i])*(c_4*sqrt(rho_F[pIDX][i]))*TimeStep;
				    rho_inplane[pIDX][i] +=  fabs(gamdot[pIDX][i])*(c_5*rho_s[pIDX][i])*TimeStep;
				dot_SSD = fabs(gamdot[pIDX][i])*(c_4*sqrt(rho_F[pIDX][i]) +
					c_6*(d_dipole)*rho_m[pIDX][i] - c_5*rho_s[pIDX][i]) - 1E-17*rho_tot[pIDX]*0.5*1E24*(d_dipole)*(D_0/bb_B[jph-1][i])*pow(rho_s[pIDX][i],2.0)*(exp(fabs(hydro[i])*bb_B[jph-1][i]*bb_B[jph-1][i]*bb_B[jph-1][i]*100/(1.38*T__K))-1.0)*exp(-1.0*q_bulk);
					//c_7*fabs(stress_climb)*pow(rho_s[pIDX][i],2.0)*pow(fabs(strainrate_climb),c_8)*exp(-1.0*q_bulk)/kT*MY_J_IN_EV;	// 1E12/m^2/s
				rho_s[pIDX][i] += dot_SSD*TimeStep;
				}
				if(istep>20000000) {
				dot_SSD =  - 4*1E14*D_0*pow(rho_s[pIDX][i],2.0)*sinh(fabs(rss[i])*bb_B[jph-1][i]*bb_B[jph-1][i]*bb_B[jph-1][i]*100/(1.38*823))*exp(-1.0*Q_bulk[jph-1]*MY_J_IN_EV/(823*BOLZ*J_IN_EV));
			//	printf("climb_rate = %le\n",dot_SSD);
					//c_7*fabs(stress_climb)*pow(rho_s[pIDX][i],2.0)*pow(fabs(strainrate_climb),c_8)*exp(-1.0*q_bulk)/kT*MY_J_IN_EV;	// 1E12/m^2/s
				rho_s[pIDX][i] += dot_SSD*TimeStep*50;
				}
				
				trial_rho_s[pIDX][i] = rho_s[pIDX][i];
        rho_dot_s[pIDX][i] = dot_SSD;
        	gam_acum[pIDX][i] = trial_gam_acum[pIDX][i];
if(istep<=200000000) {
#ifdef DD_GND
				/* GND */
				real gnd1;
				Vec3R gnd_v1;

				// Screw type
				VSet(gnd_v1, bp[jph-1][i][0], bp[jph-1][i][1], bp[jph-1][i][2]);
				gnd1 = VDot(gadot_grad[pIDX][i],gnd_v1);
#ifdef DD_GND_FIRST_ORDER
				dot_GND = gnd1/bb_B[jph-1][i]/L0;
#else
				real gnd2, gnd3;
				Vec3R gnd_v2, gnd_v3, gnd_vt1, gnd_vt2, gnd_vt3, gnd_vt4, gnd_vt5;
				
				VZero(gnd_v2);
				gnd_v3 = 0.0;
				for(j=0;j<nSYS[jph-1];j++){
					if(j!=i){
						VSet(gnd_vt1, bn[jph-1][j][0], bn[jph-1][j][1], bn[jph-1][j][2]);
						VSet(gnd_vt2, bb[jph-1][i][0], bb[jph-1][i][1], bb[jph-1][i][2]);
						VCross(gnd_vt3, gnd_vt1, gnd_vt2);
						VSet(gnd_vt4, bn[jph-1][i][0], bn[jph-1][i][1], bn[jph-1][i][2]);
						VSet(gnd_vt5, bb[jph-1][j][0], bb[jph-1][j][1], bb[jph-1][j][2]);
						VSAdd(gnd_v2, gnd_v2, trial_gam_acum[pIDX][j]*VDot(gnd_vt4,gnd_vt5),gnd_vt3);

						gnd_v3 += VDot(gnd_vt4,gnd_vt5)*VDot(ga_grad[pIDX][j],gnd_vt3);
					}
				}
				gnd2 = VDot(gadot_grad[pIDX][i],gnd_v2);
				gnd3 = gamdot[pIDX][i]*gnd_v3;
				dot_GND = (gnd1-gnd2-gnd3)/bb_B[jph-1][i]/L0;
#endif

				rho_g1[pIDX][i] += dot_GND*TimeStep;
				trial_rho_g1[pIDX][i] = rho_g1[pIDX][i];
        rho_dot_g1[pIDX][i] = dot_GND;


				// Edge type-I
				VSet(gnd_v1, bb[jph-1][i][0], bb[jph-1][i][1], bb[jph-1][i][2]);
				gnd1 = VDot(gadot_grad[pIDX][i],gnd_v1);
#ifdef DD_GND_FIRST_ORDER
				dot_GND = -1.0*gnd1/bb_B[jph-1][i]/L0;
#else
				VZero(gnd_v2);
				gnd_v3 = 0.0;
				for(j=0;j<nSYS[jph-1];j++){
					if(j!=i){
						VSet(gnd_vt1, bn[jph-1][j][0], bn[jph-1][j][1], bn[jph-1][j][2]);
						VSet(gnd_vt2, bp[jph-1][i][0], bp[jph-1][i][1], bp[jph-1][i][2]);
						VCross(gnd_vt3, gnd_vt1, gnd_vt2);
						VSet(gnd_vt4, bn[jph-1][i][0], bn[jph-1][i][1], bn[jph-1][i][2]);
						VSet(gnd_vt5, bb[jph-1][j][0], bb[jph-1][j][1], bb[jph-1][j][2]);
						VSAdd(gnd_v2, gnd_v2, trial_gam_acum[pIDX][j]*VDot(gnd_vt4,gnd_vt5),gnd_vt3);

						gnd_v3 += VDot(gnd_vt4,gnd_vt5)*VDot(ga_grad[pIDX][j],gnd_vt3);
					}
				}
				gnd2 = VDot(gadot_grad[pIDX][i],gnd_v2);
				gnd3 = gamdot[pIDX][i]*gnd_v3;
				dot_GND = (-1.0*gnd1-gnd2-gnd3)/bb_B[jph-1][i]/L0;
#endif
				rho_g2[pIDX][i] += dot_GND*TimeStep;
				trial_rho_g2[pIDX][i] = rho_g2[pIDX][i];
        rho_dot_g2[pIDX][i] = dot_GND;

				// Edge type-II
				gnd1 = 0.0;
#ifdef DD_GND_FIRST_ORDER
				dot_GND = gnd1;
#else
				VZero(gnd_v2);
				gnd_v3 = 0.0;
				for(j=0;j<nSYS[jph-1];j++){
					if(j!=i){
						VSet(gnd_vt1, bn[jph-1][j][0], bn[jph-1][j][1], bn[jph-1][j][2]);
						VSet(gnd_vt2, bn[jph-1][i][0], bn[jph-1][i][1], bn[jph-1][i][2]);
						VCross(gnd_vt3, gnd_vt1, gnd_vt2);
						VSet(gnd_vt4, bn[jph-1][i][0], bn[jph-1][i][1], bn[jph-1][i][2]);
						VSet(gnd_vt5, bb[jph-1][j][0], bb[jph-1][j][1], bb[jph-1][j][2]);
						VSAdd(gnd_v2, gnd_v2, trial_gam_acum[pIDX][j]*VDot(gnd_vt4,gnd_vt5),gnd_vt3);

						gnd_v3 += VDot(gnd_vt4,gnd_vt5)*VDot(ga_grad[pIDX][j],gnd_vt3);
					}
				}
				gnd2 = VDot(gadot_grad[pIDX][i],gnd_v2);
				gnd3 = gamdot[pIDX][i]*gnd_v3;
				dot_GND = (-1.0*gnd2-gnd3)/bb_B[jph-1][i]/L0;
#endif
				rho_g3[pIDX][i] += dot_GND*TimeStep;
				trial_rho_g3[pIDX][i] = rho_g3[pIDX][i];

				// Record accumulated shear
			
#endif
                 }
			}
		}
	}

	return;
}/*end DislocationEvolution()*/

void trial_DislocationEvolution(int idx, int jph)
{
	int i,j;
	real rss[NSYSMX],hydro[NSYSMX];	// resolved shear stress
	voigt5 sc[NSYSMX];
	ten2nd sig_aux;
	voigt sig6;
	ten4th aux3333;
	voigt66 aux66;
	real dot_SSD;	// the change rate of SSDs
	//real dot_GND;	// the change rate of GNDs
	//real gnd1, gnd2, gnd3;	// three terms in rate of GNDs
	//Vec3R gnd_v1, gnd_v2, gnd_vt1, gnd_vt2, gnd_vt3, gnd_vt4, gnd_vt5;	// 2nd terms in GND rate
	//real gnd_v3;	// 3rd term in GND rate
	real shear_G;	// shear modulus
	real Poisson_nu;	// Poisson's ratio
	real c_4, c_5, c_6, c_7, c_8;
	real q_bulk;
	real d_dipole;

	c_4 = C_4[jph-1];
	c_5 = C_5[jph-1];
	c_6 = C_6[jph-1];
	c_7 = C_7[jph-1];
	c_8 = C_8[jph-1];
	q_bulk = Q_bulk[jph-1]*MY_J_IN_EV/kT;	// reduced
	shear_G = Shear_G[jph-1];
	Poisson_nu = ElastConst_PhaseI[3];
	real stress_climb, strainrate_climb;

	/* resolve shear stress */
	T2_loop{
		sig_aux[mi][mj] = Sig[idx][mi][mj];
	}
	chg_basis(sig6,sig_aux,aux66,aux3333,2);

	for(i=0;i<nSYS[jph-1];i++){
		for(j=0;j<5;j++) sc[i][j] = Schm_gr[idx][i][j];
	}
	for(i=0;i<nSYS[jph-1];i++){
		rss[i] = sc[i][0]*sig6[0]+sc[i][1]*sig6[1]+sc[i][2]*sig6[2]+
			sc[i][3]*sig6[3]+sc[i][4]*sig6[4];
	     	hydro[i] = (sig6[0]+sig6[1]+sig6[2])/3.0;
	}

	for(i=0;i<nSYS[jph-1];i++){
		stress_climb = VonMises(Sig[idx])*3.0/2.0;
		strainrate_climb = VonMises(VelGrad[idx]);
		//d_dipole = sqrt(3.0)*shear_G*bb_B[jph-1][i]/16.0/PI/(1-Poisson_nu)/(fabs(rss[i]));	//in nm
		d_dipole = sqrt(3.0)*shear_G*bb_B[jph-1][i]/16.0/PI/(1-Poisson_nu)/(fabs(rss[i]+EPSILON));	//in nm
		if(d_dipole <= 1.6) {
		    d_dipole = 1.6;
		}
	//	rho_dipole[pIDX][i] += (fabs(gamdot[pIDX][i])*(c_6*d_dipole*rho_m[pIDX][i] - c_5*rho_s[pIDX][i]) - c_7*1.6*1E-9*1E12*rho_s[pIDX][i]*rho_s[pIDX][i])*TimeStep;
		dot_SSD = fabs(gamdot[idx][i])*(c_4*sqrt(rho_F[idx][i]) +
			c_6*(d_dipole)*rho_m[idx][i] - c_5*rho_s[idx][i]) - c_7*1.6*1E-9*1E12*rho_s[pIDX][i]*rho_s[pIDX][i] - 4*1E12*D_0*pow(rho_s[pIDX][i],2.0)*sinh(fabs(hydro[i])*bb_B[jph-1][i]*bb_B[jph-1][i]*bb_B[jph-1][i]*100/(1.38*T__K))*exp(-1.0*q_bulk);
			//c_7*fabs(stress_climb)*pow(rho_s[idx][i],2.0)*pow(fabs(strainrate_climb),c_8)*exp(-1.0*q_bulk)/kT*MY_J_IN_EV;	// 1E1/m^2/s
		if(isnan(dot_SSD)){
			printf("pyz:WTF:rank%d:rho_F[%d][%d]=%le,gamdot=%le,rss=%le,d_dipole=%le,rho_m=%le,rho_s=%le,stress_climb=%le,strainrate_climb=%le\n",
					mpirank,idx,i,rho_F[idx][i],gamdot[idx][i],rss[i],d_dipole,rho_m[idx][i],rho_s[idx][i],stress_climb,strainrate_climb);
		}
		trial_rho_s[idx][i] = rho_s[idx][i]+dot_SSD*TimeStep;
//		if(trial_rho_s[idx][i]<0){
//			trial_rho_s[idx][i] = fabs(trial_rho_s[idx][i]);
//		}
		
//
//#ifdef DD_GND
//		// Screw type
//		VSet(gnd_v1, bp[jph-1][i][0], bp[jph-1][i][1], bp[jph-1][i][2]);
//		gnd1 = VDot(gadot_grad[idx][i],gnd_v1);
//#ifdef DD_GND_FIRST_ORDER
//		dot_GND = gnd1/bb_B[jph-1][i]/L0;
//#else
//		VZero(gnd_v2);
//		gnd_v3 = 0.0;
//		for(j=0;j<nSYS[jph-1];j++){
//			if(j!=i){
//				VSet(gnd_vt1, bn[jph-1][j][0], bn[jph-1][j][1], bn[jph-1][j][2]);
//				VSet(gnd_vt2, bb[jph-1][i][0], bb[jph-1][i][1], bb[jph-1][i][2]);
//				VCross(gnd_vt3, gnd_vt1, gnd_vt2);
//				VSet(gnd_vt4, bn[jph-1][i][0], bn[jph-1][i][1], bn[jph-1][i][2]);
//				VSet(gnd_vt5, bb[jph-1][j][0], bb[jph-1][j][1], bb[jph-1][j][2]);
//				VSAdd(gnd_v2, gnd_v2, trial_gam_acum[idx][j]*VDot(gnd_vt4,gnd_vt5),gnd_vt3);
//
//				gnd_v3 += VDot(gnd_vt4,gnd_vt5)*VDot(ga_grad[idx][j],gnd_vt3);
//			}
//		}
//		gnd2 = VDot(gadot_grad[idx][i],gnd_v2);
//		gnd3 = gamdot[idx][i]*gnd_v3;
//		dot_GND = (gnd1-gnd2-gnd3)/bb_B[jph-1][i]/L0;
//#endif
//		trial_rho_g1[idx][i] = rho_g1[idx][i] + dot_GND*TimeStep;
//
//		// Edge type-I
//		VSet(gnd_v1, bb[jph-1][i][0], bb[jph-1][i][1], bb[jph-1][i][2]);
//		gnd1 = VDot(gadot_grad[idx][i],gnd_v1);
//#ifdef DD_GND_FIRST_ORDER
//		dot_GND = -1.0*gnd1/bb_B[jph-1][i]/L0;
//#else
//		VZero(gnd_v2);
//		gnd_v3 = 0.0;
//		for(j=0;j<nSYS[jph-1];j++){
//			if(j!=i){
//				VSet(gnd_vt1, bn[jph-1][j][0], bn[jph-1][j][1], bn[jph-1][j][2]);
//				VSet(gnd_vt2, bp[jph-1][i][0], bp[jph-1][i][1], bp[jph-1][i][2]);
//				VCross(gnd_vt3, gnd_vt1, gnd_vt2);
//				VSet(gnd_vt4, bn[jph-1][i][0], bn[jph-1][i][1], bn[jph-1][i][2]);
//				VSet(gnd_vt5, bb[jph-1][j][0], bb[jph-1][j][1], bb[jph-1][j][2]);
//				VSAdd(gnd_v2, gnd_v2, trial_gam_acum[idx][j]*VDot(gnd_vt4,gnd_vt5),gnd_vt3);
//
//				gnd_v3 += VDot(gnd_vt4,gnd_vt5)*VDot(ga_grad[idx][j],gnd_vt3);
//			}
//		}
//		gnd2 = VDot(gadot_grad[idx][i],gnd_v2);
//		gnd3 = gamdot[idx][i]*gnd_v3;
//		dot_GND = (-1.0*gnd1-gnd2-gnd3)/bb_B[jph-1][i]/L0;
//#endif
//		trial_rho_g2[idx][i] = rho_g2[idx][i] + dot_GND*TimeStep;
//
//		// Edge type-II
//		gnd1 = 0.0;
//#ifdef DD_GND_FIRST_ORDER
//		dot_GND = gnd1;
//#else
//		VZero(gnd_v2);
//		gnd_v3 = 0.0;
//		for(j=0;j<nSYS[jph-1];j++){
//			if(j!=i){
//				VSet(gnd_vt1, bn[jph-1][j][0], bn[jph-1][j][1], bn[jph-1][j][2]);
//				VSet(gnd_vt2, bn[jph-1][i][0], bn[jph-1][i][1], bn[jph-1][i][2]);
//				VCross(gnd_vt3, gnd_vt1, gnd_vt2);
//				VSet(gnd_vt4, bn[jph-1][i][0], bn[jph-1][i][1], bn[jph-1][i][2]);
//				VSet(gnd_vt5, bb[jph-1][j][0], bb[jph-1][j][1], bb[jph-1][j][2]);
//				VSAdd(gnd_v2, gnd_v2, trial_gam_acum[idx][j]*VDot(gnd_vt4,gnd_vt5),gnd_vt3);
//
//				gnd_v3 += VDot(gnd_vt4,gnd_vt5)*VDot(ga_grad[idx][j],gnd_vt3);
//			}
//		}
//		gnd2 = VDot(gadot_grad[idx][i],gnd_v2);
//		gnd3 = gamdot[idx][i]*gnd_v3;
//		dot_GND = (-1.0*gnd2-gnd3)/bb_B[jph-1][i]/L0;
//#endif
//		trial_rho_g3[idx][i] = rho_g3[idx][i] + dot_GND*TimeStep;
//#endif
	}

	return;
}/*end trial_DislocationEvolution()*/

void AvgDislocation(void)
{
	real rho_local;
	int i, jph;

	// Forest dislocation densities
	rho_local=0.0;
	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			for(i=0;i<nSYS[jph-1];i++){
				rho_local += rho_F[pIDX][i];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&rho_local, &rho_F_avg, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
	rho_F_avg *= 1.0/Nxyz;

	// Parallel dislocation densities
	rho_local=0.0;
	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			for(i=0;i<nSYS[jph-1];i++){
				rho_local += rho_P[pIDX][i];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&rho_local, &rho_P_avg, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
	rho_P_avg *= 1.0/Nxyz;

	// SSD  dislocation densities
	rho_local=0.0;
	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			for(i=0;i<nSYS[jph-1];i++){
				rho_local += rho_s[pIDX][i];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&rho_local, &rho_s_avg, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
	rho_s_avg *= 1.0/Nxyz;

	// Mobile dislocation densities
	rho_local=0.0;
	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			for(i=0;i<nSYS[jph-1];i++){
				rho_local += rho_m[pIDX][i];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&rho_local, &rho_m_avg, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
	rho_m_avg *= 1.0/Nxyz;

	// GNDs
	rho_local=0.0;
	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			for(i=0;i<nSYS[jph-1];i++){
				rho_local += sqrt(rho_g1[pIDX][i]*rho_g1[pIDX][i]+rho_g2[pIDX][i]*rho_g2[pIDX][i]+rho_g3[pIDX][i]*rho_g3[pIDX][i]);
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&rho_local, &rho_g_avg, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
	rho_g_avg *= 1.0/Nxyz;
    
    	// dipole dislocation densities
	rho_local=0.0;
	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			for(i=0;i<nSYS[jph-1];i++){
				rho_local += rho_dipole[pIDX][i];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&rho_local, &rho_dipole_avg, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
	rho_dipole_avg *= 1.0/Nxyz;

    	// forest generation dislocation densities
	rho_local=0.0;
	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			for(i=0;i<nSYS[jph-1];i++){
				rho_local += rho_forest[pIDX][i];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&rho_local, &rho_forest_avg, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
	rho_forest_avg *= 1.0/Nxyz;
	
	// inplane annihilation
		rho_local=0.0;
	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			for(i=0;i<nSYS[jph-1];i++){
				rho_local += rho_inplane[pIDX][i];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&rho_local, &rho_inplane_avg, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
	rho_inplane_avg *= 1.0/Nxyz;
	

	return;
}/*end AvgDislocation() */

#endif
