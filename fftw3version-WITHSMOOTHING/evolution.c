#include "evp.h"

static void InteractionStrengthMatrix(char *type, real self, real coplane,
		real cross, real glissile, real hirth, real lomer, real **chi)
{
	/* The assignment of interaction strength matrix dependes on the label
	   of the slip system. Currently only FCC is considered according to
	  "Arsenlis and Parks, JMPS, 50(2002), 1979" */
	int i,j;
	if(strcmp(type, "FCC")==0||strcmp(type, "fcc")==0){
		chi[0][0] = self;
		chi[0][1] = coplane;
		chi[0][2] = coplane;
		chi[0][3] = hirth;
		chi[0][4] = glissile;
		chi[0][5] = lomer;
		chi[0][6] = cross;
		chi[0][7] = glissile;
		chi[0][8] = glissile;
		chi[0][9] = hirth;
		chi[0][10] = lomer;
		chi[0][11] = glissile;

		chi[1][1] = self;
		chi[1][2] = coplane;
		chi[1][3] = glissile;
		chi[1][4] = cross;
		chi[1][5] = glissile;
		chi[1][6] = glissile;
		chi[1][7] = hirth;
		chi[1][8] = lomer;
		chi[1][9] = lomer;
		chi[1][10] = hirth;
		chi[1][11] = glissile;

		chi[2][2] = self;
		chi[2][3] = lomer;
		chi[2][4] = glissile;
		chi[2][5] = hirth;
		chi[2][6] = glissile;
		chi[2][7] = lomer;
		chi[2][8] = hirth;
		chi[2][9] = glissile;
		chi[2][10] = glissile;
		chi[2][11] = cross;

		chi[3][3] = self;
		chi[3][4] = coplane;
		chi[3][5] = coplane;
		chi[3][6] = hirth;
		chi[3][7] = lomer;
		chi[3][8] = glissile;
		chi[3][9] = cross;
		chi[3][10] = glissile;
		chi[3][11] = glissile;

		chi[4][4] = self;
		chi[4][5] = coplane;
		chi[4][6] = lomer;
		chi[4][7] = hirth;
		chi[4][8] = glissile;
		chi[4][9] = glissile;
		chi[4][10] = hirth;
		chi[4][11] = lomer;

		chi[5][5] = self;
		chi[5][6] = glissile;
		chi[5][7] = glissile;
		chi[5][8] = cross;
		chi[5][9] = glissile;
		chi[5][10] = lomer;
		chi[5][11] = hirth;

		chi[6][6] = self;
		chi[6][7] = coplane;
		chi[6][8] = coplane;
		chi[6][9] = hirth;
		chi[6][10] = glissile;
		chi[6][11] = lomer;

		chi[7][7] = self;
		chi[7][8] = coplane;
		chi[7][9] = glissile;
		chi[7][10] = cross;
		chi[7][11] = glissile;

		chi[8][8] = self;
		chi[8][9] = lomer;
		chi[8][10] = glissile;
		chi[8][11] = hirth;

		chi[9][9] = self;
		chi[9][10] = coplane;
		chi[9][11] = coplane;

		chi[10][10] = self;
		chi[10][11] = coplane;

		chi[11][11] = self;

		/* symmetric matrix */
		for(i=1;i<12;i++){
			for(j=0;j<i;j++){
				chi[i][j] = chi[j][i];
			}
		}
	}
	else{
		PError("Interaction strength matrix is only supported for FCC currently!!", 1110);	
	}
	return;
}/*end InteractionStrengthMatrix()*/


static real NR_quality(voigt sig, voigt xlambda, voigt eps, voigt strain, voigt66 sg, int idx,int jphi)
{
	voigt edotp;
	voigt66 d_edotp;
	voigt tot_eps;	// total strain
  voigt res;
  real f;
  int i,j;
#ifdef DD_BASED_FLAG
#ifdef DD_POWER_LAW
				StrainRate_Orowan_POWER(sig,edotp,d_edotp,idx,jphi);
#else
				StrainRate_Orowan(sig,edotp,d_edotp,idx,jphi);
#endif
#else
				StrainRate_eval(sig,edotp, d_edotp,idx, jphi);
#endif

				/* tot_eps is the total  strain and eps6 is the 
				   current plastic strain, i.e. Eq. 4 */
				for(i=0;i<6;i++){
					tot_eps[i] = eps[i] + edotp[i]*TimeStep;
					for(j=0;j<6;j++){
						tot_eps[i] += sg[i][j]*sig[j];
					}
				}

				// calculate the residual R, Eq. 16
				for(i=0;i<6;i++){
					res[i] = sig[i] - xlambda[i];
					for(j=0;j<6;j++){
						res[i] += C066[i][j]*(tot_eps[j]-strain[j]);
					}
				}
        for(f=0.0,i=0;i<6;i++){
          f += res[i]*res[i];
        }

        return f/2.0;
}/*end NR_quality()*/


static void UpdateStress(int istep, real *Err_e_local, real *Err_s_local, int update_flag)
{
	/* Solve the stress using the updated strain during the iteration for
	   each single mechanical step (istep), i.e. Eq. 16 */ 
	int jph;
	voigt66 sg66;
	ten2nd xlambda_aux, sig_aux, eps_aux, strain_aux;
	voigt xlambda6, sig6, eps6, strain6;
	voigt edotp6;
	ten2nd edotp_aux;
	voigt66 d_edotp66;
	voigt tot_eps;	// total strain
	voigt sig6_old;
	ten4th aux3333;
	voigt66 aux66;
	real signorm,enorm;
	real erroral, erral;
	real conv_NR, conv_istep_e, conv_istep_s;

	voigt res;	// residual R to be nullified
	voigt66	jacob_inv;	// Jacobian of R
	int itmaxal, iterl;
	int i,j,k;

	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			C6_loop{
				sg66[mi][mj] = C_gr[pIDX][mi][mj];
			}
			LU_inv_66(sg66);
			
			T2_loop{
				/* the iteration starts with the current stress field */
				xlambda_aux[mi][mj] = Sig[pIDX][mi][mj];
				sig_aux[mi][mj] = Sig[pIDX][mi][mj];
				/* plastic strain at time t. This is used in Eq. 4 to update total strain
				 from calculated strain rate*/
				eps_aux[mi][mj] = Eps[pIDX][mi][mj];
				/* DisGrad stores the updated displacement gradient obtained from Eq. 15.
				   So here strain_aux/strain6 stores the updated total strain*/
				strain_aux[mi][mj] = (DisGrad[pIDX][mi][mj]+DisGrad[pIDX][mj][mi])/2.0;
			}
			chg_basis(xlambda6,xlambda_aux,aux66,aux3333,2);
			chg_basis(sig6,sig_aux,aux66,aux3333,2);
			chg_basis(eps6,eps_aux,aux66,aux3333,2);
			chg_basis(strain6,strain_aux,aux66,aux3333,2);

			signorm = 0.0;
			T2_loop{
				signorm += xlambda_aux[mi][mj]*xlambda_aux[mi][mj];
			}
			signorm = sqrt(signorm);
            
            enorm = 0.0;
			T2_loop{
				enorm += eps_aux[mi][mj]*eps_aux[mi][mj];
			}
			enorm = sqrt(enorm);
            
			erroral = 1E-10;
			itmaxal = 500;
			iterl = 0;
			erral = 10*erroral;

			/* Newton-Raphson method to solve augmented Lagrangians */
			while((iterl<itmaxal)&&(fabs(erral)>fabs(erroral))){
				iterl++;

				for(i=0;i<6;i++) sig6_old[i] = sig6[i];

				/* Update the plastic strain rate based on current stress */
#ifdef DD_BASED_FLAG
#ifdef DD_POWER_LAW
				StrainRate_Orowan_POWER(sig6,edotp6,d_edotp66,pIDX,jph);
#else
				StrainRate_Orowan(sig6,edotp6,d_edotp66,pIDX,jph);
//exit(0);
#endif
#else
				StrainRate_eval(sig6,edotp6, d_edotp66,pIDX, jph);
#endif

				/* tot_eps is the total  strain and eps6 is the 
				   current plastic strain, i.e. Eq. 4 */
				for(i=0;i<6;i++){
					tot_eps[i] = eps6[i] + edotp6[i]*TimeStep;
					for(j=0;j<6;j++){
						tot_eps[i] += sg66[i][j]*sig6[j];
					}
				}

				// calculate the residual R, Eq. 16
				for(i=0;i<6;i++){
					res[i] = sig6[i] - xlambda6[i];
					for(j=0;j<6;j++){
						res[i] += C066[i][j]*(tot_eps[j]-strain6[j]);
					}
				}
				// calculate the Jacobian of R
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
						jacob_inv[i][j] = (real)(i==j);
						for(k=0;k<6;k++){
							// Eq. 18
							jacob_inv[i][j] += C066[i][k]*(sg66[k][j]+d_edotp66[k][j]*TimeStep);
						}
					}
				}

#ifdef NR_MODIFIED
	      voigt66	jacob;	// Jacobian of R, not the inverse!
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
            jacob[i][j] = jacob_inv[i][j];
          }
        }
#endif
        
				LU_inv_66(jacob_inv);
//printf("residual = %le and jacob = %le\n",res[1],jacob_inv[5][1]);
				// Newton-Raphson update
#ifdef NR_MODIFIED
        /* Using line searches and backtracking to ensure global
         * convergence */

        // current f=F*F/2 value
        real f_old = 0.0;
        for(i=0;i<6;i++){
          f_old += res[i]*res[i];
        }
        f_old *= 0.5;

        // current gradient of f
        real gradf[6];
				for(i=0;i<6;i++){
          gradf[i] = 0.0;
					for(j=0;j<6;j++){
            gradf[i] += res[j]*jacob[j][i];
          }
        }

        // full NR step
        real NR_step[6];
				for(i=0;i<6;i++){
          NR_step[i] = 0.0;
					for(j=0;j<6;j++){
						NR_step[i] -= jacob_inv[i][j]*res[j];
					}
				}

        real stpmax = 10.0;  // limit of the length
        real sum, slope, test, temp, tmplam;
        real a_coef, alam, alam2, alamin, b_coef, disc, f2, rhs1, rhs2;
        real ALF = 1.E-4;
        real TOLX = 1.E-7;
        for(sum=0.0,i=0;i<6;i++) sum += NR_step[i]*NR_step[i];
        sum = sqrt(sum);
        if(sum>stpmax){
          for(i=0;i<6;i++) NR_step[i] *= stpmax/sum;
        }
        for(slope=0.0,i=0;i<6;i++){
          slope += gradf[i]*NR_step[i];
        }

       // printf("slope is = %le\n",slope);
        if(slope>=0){
           printf("NR_step is = %le %le %le %le %le %le\n", NR_step[0],NR_step[1],NR_step[2],NR_step[3],NR_step[4],NR_step[5]);
           printf("gradf is =  %le %le %le %le %le %le\n", gradf[0],gradf[1],gradf[2],gradf[3],gradf[4],gradf[5]);
          PError("Roundoff problem during line searching in NR.", 1203);
        }
        test = 0.0;
        for(i=0;i<6;i++){
          temp = fabs(NR_step[i])/MAX(fabs(sig6_old[i]),1.0);
          if(temp>test) test = temp;
        }
        alamin = TOLX/test;
        alam = 1.0;
        while(1){
          // always try full NR step first
          for(i=0;i<6;i++) sig6[i] = sig6_old[i] + alam*NR_step[i];
          real f_new = NR_quality(sig6, xlambda6, eps6, strain6, sg66, pIDX, jph);
          if(alam<alamin){  // convergence on delta_x
            for(i=0;i<6;i++) sig6[i] = sig6_old[i];
            break;
          }else if(f_new <= f_old+ALF*alam*slope){  // sufficient function decrease
            break;
          }
          else{ // backtrack
            if(fabs(alam-1.0)<1.E-4){ // first time
              tmplam = -slope/(2.0*(f_new-f_old-slope));
            }else{
              rhs1 = f_new-f_old-alam*slope;
              rhs2 = f2-f_old-alam2*slope;
              a_coef = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
              b_coef = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
              if(fabs(a_coef)<1E-4) tmplam = -slope/(2.0*b_coef);
              else{
                disc = b_coef*b_coef-3.0*a_coef*slope;
                if(disc<0.0) tmplam = 0.5*alam;
                else if(b_coef<=0.0) tmplam = (-b_coef+sqrt(disc))/(3.0*a_coef);
                else tmplam = -slope/(b_coef+sqrt(disc));
              }
              if(tmplam>0.5*alam)
                tmplam = 0.5*alam;
            }
          }
          alam2=alam;
          alam = MAX(tmplam,0.1*alam);
        }
        
#else
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
						sig6[i] -= jacob_inv[i][j]*res[j];
					}
				}
#endif

				// convergence check
				conv_NR= 0.0;
				conv_istep_e= 0.0;
				conv_istep_s = 0.0;
				for(i=0;i<6;i++){
					conv_NR += (sig6[i]-sig6_old[i])*(sig6[i]-sig6_old[i]);
					conv_istep_e += (tot_eps[i]-strain6[i])*(tot_eps[i]-strain6[i]);
					conv_istep_s += (sig6[i]-xlambda6[i])*(sig6[i]-xlambda6[i]);
				}
				erral = sqrt(conv_NR)/signorm;
                conv_istep_s = sqrt(conv_istep_s);
                conv_istep_e = sqrt(conv_istep_e);
				// update crss
				if((Hard_Flag==1)){
#ifdef DD_BASED_FLAG
					//trial_DislocationEvolution(pIDX,jph);
#else
					/* because the shear rate that is about to use
					   is calcualted (in StrainRate_eval()) based on
					   trial stress field, we use the get_trialtau()
					   subroutine which is a trial versio nof harden() */
					//get_trialtau(pIDX,jph);
					//get_trialtau_anal(pIDX,jph);
#endif
				}
			}	// end of while() loop


			chg_basis(sig6,sig_aux,aux66,aux3333,1);
			chg_basis(edotp6,edotp_aux,aux66,aux3333,1);
			// update stress and strain rate fields
			T2_loop{
				Sig[pIDX][mi][mj] = sig_aux[mi][mj];
				Edot[pIDX][mi][mj] = edotp_aux[mi][mj];
			}
			*Err_s_local += conv_istep_s;
			*Err_e_local += conv_istep_e;
		}
		else{
			T2_loop{
				Sig[pIDX][mi][mj] = 0.0;
				Edot[pIDX][mi][mj] = 0.0;
			}
		}
	}

	return;
}/*end UpdateStress()*/

static void UpdateStress_DRX(int istep, real *Err_e_local, real *Err_s_local, int update_flag)
{
	/* Solve the stress using the updated strain during the iteration for
	   each single mechanical step (istep), i.e. Eq. 16 */ 
	int jph;
	voigt66 sg66;
	ten2nd xlambda_aux, sig_aux, eps_aux, strain_aux;
	voigt xlambda6, sig6, eps6, strain6;
	voigt edotp6;
	ten2nd edotp_aux;
	voigt66 d_edotp66;
	voigt tot_eps;	// total strain
	voigt sig6_old;
	ten4th aux3333;
	voigt66 aux66;
	real signorm;
	real erroral, erral;
	real conv_NR, conv_istep_e, conv_istep_s;

	voigt res;	// residual R to be nullified
	voigt66	jacob_inv;	// Jacobian of R
	int itmaxal, iterl;
	int i,j,k;

	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			C6_loop{
				sg66[mi][mj] = C_gr[pIDX][mi][mj];
			}
			LU_inv_66(sg66);
			
			T2_loop{
				/* the iteration starts with the current stress field */
				xlambda_aux[mi][mj] = Sig[pIDX][mi][mj];
				sig_aux[mi][mj] = Sig[pIDX][mi][mj];
				/* plastic strain at time t. This is used in Eq. 4 to update total strain
				 from calculated strain rate*/
				eps_aux[mi][mj] = Eps[pIDX][mi][mj];
				/* DisGrad stores the updated displacement gradient obtained from Eq. 15.
				   So here strain_aux/strain6 stores the updated total strain*/
				strain_aux[mi][mj] = (DisGrad[pIDX][mi][mj]+DisGrad[pIDX][mj][mi])/2.0;
			}
			chg_basis(xlambda6,xlambda_aux,aux66,aux3333,2);
			chg_basis(sig6,sig_aux,aux66,aux3333,2);
			chg_basis(eps6,eps_aux,aux66,aux3333,2);
			chg_basis(strain6,strain_aux,aux66,aux3333,2);

			signorm = 0.0;
			T2_loop{
				signorm += xlambda_aux[mi][mj]*xlambda_aux[mi][mj];
			}
			signorm = sqrt(signorm);

			erroral = 1E-7;
			itmaxal = 100;
			iterl = 0;
			erral = 10*erroral;

			/* Newton-Raphson method to solve augmented Lagrangians */
			while((iterl<itmaxal)&&(fabs(erral)>fabs(erroral))){
				iterl++;

				for(i=0;i<6;i++) sig6_old[i] = sig6[i];

//				/* Update the plastic strain rate based on current stress */
//#ifdef DD_BASED_FLAG
//#ifdef DD_POWER_LAW
//				StrainRate_Orowan_POWER(sig6,edotp6,d_edotp66,pIDX,jph);
//#else
//				StrainRate_Orowan(sig6,edotp6,d_edotp66,pIDX,jph);
//#endif
//#else
//				StrainRate_eval(sig6,edotp6, d_edotp66,pIDX, jph);
//#endif

				/* tot_eps is the total  strain and eps6 is the 
				   current plastic strain, i.e. Eq. 4 */
				for(i=0;i<6;i++){
					//tot_eps[i] = eps6[i] + edotp6[i]*TimeStep;
					tot_eps[i] = eps6[i];
					for(j=0;j<6;j++){
						tot_eps[i] += sg66[i][j]*sig6[j];
					}
				}

				// calculate the residual R, Eq. 16
				for(i=0;i<6;i++){
					res[i] = sig6[i] - xlambda6[i];
					for(j=0;j<6;j++){
						res[i] += C066[i][j]*(tot_eps[j]-strain6[j]);
					}
				}
				// calculate the Jacobian of R
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
						jacob_inv[i][j] = (real)(i==j);
						for(k=0;k<6;k++){
							// Eq. 18
							//jacob_inv[i][j] += C066[i][k]*(sg66[k][j]+d_edotp66[k][j]*TimeStep);
							jacob_inv[i][j] += C066[i][k]*(sg66[k][j]);
						}
					}
				}
				LU_inv_66(jacob_inv);

				// Newton-Raphson update
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
						sig6[i] -= jacob_inv[i][j]*res[j];
					}
				}

				// convergence check
				conv_NR= 0.0;
				conv_istep_e= 0.0;
				conv_istep_s = 0.0;
				for(i=0;i<6;i++){
					conv_NR += (sig6[i]-sig6_old[i])*(sig6[i]-sig6_old[i]);
					conv_istep_e += (tot_eps[i]-strain6[i])*(tot_eps[i]-strain6[i]);
					conv_istep_s += (sig6[i]-xlambda6[i])*(sig6[i]-xlambda6[i]);
				}
				erral = conv_NR/signorm;

				// update crss
				if((Hard_Flag==1)&&(istep>2)&&(update_flag==1)){
#ifdef DD_BASED_FLAG
					trial_DislocationEvolution(pIDX,jph);
#else
					/* because the shear rate that is about to use
					   is calcualted (in StrainRate_eval()) based on
					   trial stress field, we use the get_trialtau()
					   subroutine which is a trial versio nof harden() */
					//get_trialtau(pIDX,jph);
					//get_trialtau_anal(pIDX,jph);
#endif
				}
			}	// end of while() loop


			chg_basis(sig6,sig_aux,aux66,aux3333,1);
			chg_basis(edotp6,edotp_aux,aux66,aux3333,1);
			// update stress and strain rate fields
			T2_loop{
				Sig[pIDX][mi][mj] = sig_aux[mi][mj];
				Edot[pIDX][mi][mj] = edotp_aux[mi][mj];
			}
			*Err_s_local += conv_istep_s;
			*Err_e_local += conv_istep_e;
		}
		else{
			T2_loop{
				Sig[pIDX][mi][mj] = 0.0;
				Edot[pIDX][mi][mj] = 0.0;
			}
		}
	}

	return;
}/*end UpdateStress_DRX()*/

static void get_smacro(void)
{
	real local_sigavg, local_sigavg1;
	voigt sav6;
	voigt5 sav5;
	voigt66 aux66;
	ten4th aux3333;
	int i, ii,jj, k, kk,ll;

	int IJV[2][6] = {{0,1,2,1,0,0,},{0,1,2,2,2,1}};

	// overal stress
	T2_loop{
		local_sigavg = 0.0;
		local_sigavg1 = 0.0;
		local_loop{
			local_sigavg += Sig[pIDX][mi][mj];
			if(phase_f[pIDX]==1)
				local_sigavg1 += Sig[pIDX][mi][mj];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&local_sigavg, &SigAvg[mi][mj], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&local_sigavg1, &SigAvg1[mi][mj], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		SigAvg[mi][mj] *= WGT;
		SigAvg1[mi][mj] *= WGT/Wgt_ph1;
	}

	for(i=0;i<6;i++){
		ii=IJV[0][i];
		jj=IJV[1][i];
		dDisGradAvg[ii][jj] = 0.0;
		if(VelGrad_BC_Flag[i]==0){	// the component controlled by stress, implying a disp. variation
			for(k=0;k<6;k++){
				kk=IJV[0][k];
				ll=IJV[1][k];
				dDisGradAvg[ii][jj] += S0[ii][jj][kk][ll]*Stress_BC_Flag[k]*
					(Scauchy[kk][ll]-SigAvg[kk][ll]);
			}
		}
	}
	T2_loop{
		dDisGradAvg_acum[mi][mj] += dDisGradAvg[mi][mj];
	}
	if(mpirank==0){
	//	printf("dDisGradAvg(1,1),(2,2) = %e,%e\n",dDisGradAvg[0][0],dDisGradAvg[1][1]);
	}
	chg_basis(sav6,SigAvg,aux66,aux3333,2);
	for(i=0;i<5;i++) sav5[i] = sav6[i];
	chg_basis5(sav5,SigDevAvg,aux66,aux3333,1);

	s_vm = 0.0;
	T2_loop{
	    if(mi!=mj){
		s_vm += SigDevAvg[mi][mj]*SigDevAvg[mi][mj];
	    }
	    if(mi==mj){
	        s_vm += (SigDevAvg[mi][mj] - (SigDevAvg[1][1] + SigDevAvg[2][2] + SigDevAvg[3][3])/3.0)*(SigDevAvg[mi][mj] - (SigDevAvg[1][1] + SigDevAvg[2][2] + SigDevAvg[3][3])/3.0);
	    }
	}
	s_vm = sqrt(3./2.*s_vm);	// for stress, it is 3/2

	chg_basis(sav6,SigAvg1,aux66,aux3333,2);
	for(i=0;i<5;i++) sav5[i] = sav6[i];
	chg_basis5(sav5,SigDevAvg,aux66,aux3333,1);

	s_vm1 = 0.0;
	T2_loop{
		s_vm1 += SigDevAvg[mi][mj]*SigDevAvg[mi][mj];
	}
	s_vm1 = sqrt(3./2.*s_vm1);
				
	return;
}/*end get_smacro()*/


void Evolution(void)
{
   real tmp_smooth[512];
   ten2nd *sig_smooth;
	int s;
	int step_drx,stop_check,count,tmp_flag;
    real dd_temp,dd_ave;
    int mesh_count,g_count,i,jph;
    double rho_avg[60000];
    double ph,th,om;
    double tmp[3];
    real tmp1;
    real tmp_re[3];
    real tmp_im[3];
	int istep;		// step # of deformation test
	int iter;
        int ntimes;// step # of iteration (within a given istep)
	ten2nd DisGradAvg_t = {0.0};	// store the actual macro disp. grad. at time t
	ten2nd DisGradAvg_actual = {0.0};
	voigt66 c066_local;
	ten2nd sym_du_r, sym_du_i;
	
    //int nph1, nph1_all;
  	voigt aux6;
  	ten2nd aux33;
	ten2nd sa2xt;
	 std::vector<G_Info>::iterator it, end;
	 //voigt66 c066_local = {0.0};
    //int nph1, nph1_all;
  	//voigt aux6;
	//ten2nd aux33;
    step_drx = 0;
    double Identity[3][3];
    count = 0;
    ntimes = 1;
	
  for(iter=0;iter<60000;iter++) {
rho_avg[iter] = 0.0;
}  
  /* 9/16/15 --- PYZ 
   * Add a list of isteps (involving PF relaxation) to be recorded for further analysis */
  //int RC_LIST[14] = {2259,3497,3684,4499,4908,5480,6499,7480,8502,9500,10400,7879,7910,5276};
  int RC_LIST[6] = {9604,9605,9606,9607,9608,9610};
  int RC_length = 6;

	ten2nd EpsAvg_local, EdotAvg_local;
	real evmp, dvmp;
	real Err_e_local, Err_s_local;


	if(mpirank==0){
		printf("\n\n======================================\n");
		printf("-------------Simulation starts----------------\n");
		
	}


local_loop{
	    new_position[pIDX][0] = px+lxs;
	    new_position[pIDX][1] = py;
	    new_position[pIDX][2] = pz;
	}

	for(istep=0;istep<N_steps;istep++){
		if(mpirank==0){
			printf("\n****************************************\n");
			printf("STEP = %d\n",istep);
			if(N_steps!=1){
				fprintf(fp_err,"STEP = %d\n",istep);
			}
		}
		if(CREEP_FLAG==1&&e_vm>0.3){
			break;
		}

		

		/* "dDisGradAvg" is the local strain deviation part from
		   the average, i.e. "E-<e>". It is calculated in get_smacro()
		   according to Eq. 23 (the second term on the l.h.s). It is
		   calculated for each mechanical step iteration, and "dDisGradAvg_acum"
		   is the accumulated value after all the iterations (the following
		   while loop) for the current mechanical step (istep). */
 		T2_loop{
			if(mi==mj){
    Identity[mi][mj] = 1.0;
} else{
    Identity[mi][mj] = 0.0;
}
 		}
// 
// 		/***********************************
// 		  iteration to update stress at t+dt
// 		  ***********************************/
// 		iter = 0;
// 		Err_e = 2.*Err;
// 		Err_s = 2.*Err;
// 		
// 		while((iter<IterMax)&&((fabs(Err_s)>fabs(Err))||(fabs(Err_e)>fabs(Err)))){
// 			iter++;
// 			if(mpirank==0){
// 				printf("\nITER = %d\n",iter);
// 				printf("Forward FFT of stress field\n");
// 			}
// 			// k-space stress field
// 			T2_loop{
// 				local_loop{
// 					fft_data[pIDX].re = Sig[pIDX][mi][mj];
// 					fft_data[pIDX].im = 0.0;
// 				}
// 				MPI_Barrier(MPI_COMM_WORLD);
// 				fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
// 				local_loop{
// 					kSig_r[pIDX][mi][mj] = fft_data[pIDX].re;
// 					kSig_i[pIDX][mi][mj] = fft_data[pIDX].im;
// 				}
// 			}
// 			local_loop{
// 				T2_loop{
// 					sym_du_r[mi][mj] = 0.0;;
// 					sym_du_i[mi][mj] = 0.0;;
// 					T2p_loop{
// 						sym_du_r[mi][mj] += GAMMA[pIDX][mi][mj][mip][mjp]*kSig_r[pIDX][mip][mjp];
// 						sym_du_i[mi][mj] += GAMMA[pIDX][mi][mj][mip][mjp]*kSig_i[pIDX][mip][mjp];
// 					}
// 				}
// 				T2_loop{
// 					kSig_r[pIDX][mi][mj] = sym_du_r[mi][mj];
// 					kSig_i[pIDX][mi][mj] = sym_du_i[mi][mj];
// 				}
// 
// 			}
// 
// 			if(mpirank==0){
// 				printf("Inverse FFT to get strain field\n");
// 			}
// 			// update strain field in real space
// 			T2_loop{
// 				local_loop{
// 					fft_data[pIDX].re = kSig_r[pIDX][mi][mj];
// 					fft_data[pIDX].im = kSig_i[pIDX][mi][mj];
// 				}
// 				MPI_Barrier(MPI_COMM_WORLD);
// 				fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
// 				local_loop{
// 					DisGrad[pIDX][mi][mj] += dDisGradAvg[mi][mj] + fft_data[pIDX].re/Nxyz;	// Eq. 15
// 				}
// 			}
// 
// 			if(mpirank==0){
// 				printf("Update stress field\n");
// 			}
// 			Err_e_local = 0.0;
// 			Err_s_local = 0.0;
// //#ifdef DD_BASED_FLAG
// //#ifdef DD_GND
// //			// graidents of shear rate and shear
// //			Gradient_ShearRate();
// //			Gradient_Shear();
// //#endif
// //#endif
// 			// update stress, which requires Newton-Raphson method
// 			// NOTE: N-R is run locally. So PEs could run for different iteration steps
// 			UpdateStress(istep, &Err_e_local, &Err_s_local, 1);	
// 			// collect errors
// 			MPI_Barrier(MPI_COMM_WORLD);
// 			MPI_Allreduce(&Err_e_local, &Err_e, 1, MPI_real,
// 					MPI_SUM, MPI_COMM_WORLD);
// 			MPI_Allreduce(&Err_s_local, &Err_s, 1, MPI_real,
// 					MPI_SUM, MPI_COMM_WORLD);
// 			Err_e *= WGT;
// 			Err_s *= WGT;
// 
// 			if(mpirank==0){
// 				printf("ERRE = %e\n",Err_e);
// 				printf("ERRS = %e\n",Err_s);
// 			}
// 
// 			// update the average quantities.
// 			get_smacro();
// 
// 			Err_e /= e_vm;	
// 			Err_s /= s_vm;
// 			if(mpirank==0){
// 				printf("Strain field error = %f\n", Err_e);
// 				printf("Stress field error = %e\n", Err_s);
// 				fprintf(fp_err,"%d\t%e\t%e\t%e\n",iter,Err_e,Err_s,s_vm);
//         fflush(fp_err);
// 			}
// 
// 		}/* while() loop ends and stress converges */
// 
// 		/* update the velocity gradient field. VelGrad on the r.h.s.
// 		   stores the previous displ. gradient fiels */
// 		local_loop{
// 			T2_loop{
// 				VelGrad[pIDX][mi][mj] = (DisGrad[pIDX][mi][mj] - VelGrad[pIDX][mi][mj])/TimeStep;
// 			}
// 		}
// 		T2_loop{
// 			// dDisGradAvg_acum is updated in get_smacro()
// 			/* DisGradAvg is always the initially applied strain based on linear ealstic assumption
// 			   for each mechanical step (istep), and dDisGradAvg_acum is the resulted adjustement
// 			   after the iteration, which depends on the boundary conditon. */
// 			DisGradAvg_actual[mi][mj] = DisGradAvg[mi][mj] + dDisGradAvg_acum[mi][mj];
// 			VelGradAvg[mi][mj] = (DisGradAvg_actual[mi][mj]-DisGradAvg_t[mi][mj])/TimeStep;
// 		}
// 
// 		if(mpirank==0){
// 			printf("DisGradAvg(1,1),DisGradAvg(2,2),DisGradAvg(3,3)\n");
// 			printf("%e,%e,%e\n",DisGradAvg_actual[0][0], DisGradAvg_actual[1][1], DisGradAvg_actual[2][2]);
// 			printf("DisGradAvg(1,1)/DisGradAvg(3,3)\n");
// 			printf("%e\n",DisGradAvg_actual[0][0]/DisGradAvg_actual[2][2]);
// 			printf("SigAvg(1,1),SigAvg(2,2),SigAvg(3,3)\n");
// 			printf("%e,%e,%e\n", SigAvg[0][0],SigAvg[1][1],SigAvg[2][2]);
// 		}
// 		e_vm = VonMises(DisGradAvg_actual);
// 		d_vm = VonMises(VelGradAvg);
// 		TimeTot += TimeStep;
// 		T2_loop{
// 			DisGradAvg_t[mi][mj] = DisGradAvg_actual[mi][mj];
// 		}
//      
// 		// update strain field
// 		local_loop{
// 			T2_loop{
// 				Eps[pIDX][mi][mj] += Edot[pIDX][mi][mj]*TimeStep;	// Edot was updated in UpdateStress()
// 			}
// 		}
// 
// 
// 		
// 
//     	// Plastic VM
// 		T2_loop{
// 			EpsAvg_local[mi][mj] = 0.0;
// 			EdotAvg_local[mi][mj] = 0.0;
// 		}
// 		local_loop{
// 			T2_loop{
// 				EpsAvg_local[mi][mj] += Eps[pIDX][mi][mj]*WGT;
// 				EdotAvg_local[mi][mj] += Edot[pIDX][mi][mj]*WGT;
// 			}
// 		}
// 		T2_loop{
// 			MPI_Barrier(MPI_COMM_WORLD);
// 			MPI_Allreduce(&EpsAvg_local[mi][mj], &EpsAvg[mi][mj], 1, MPI_real,
// 					MPI_SUM, MPI_COMM_WORLD);
// 			MPI_Allreduce(&EdotAvg_local[mi][mj], &EdotAvg[mi][mj], 1, MPI_real,
// 					MPI_SUM, MPI_COMM_WORLD);
// 		}
// 
// 		evmp = 0.0;
// 		dvmp = 0.0;
// 		T2_loop{
// 			evmp += EpsAvg[mi][mj]*EpsAvg[mi][mj];
// 			dvmp += EdotAvg[mi][mj]*EdotAvg[mi][mj];
// 		}
// 		evmp = sqrt(2./3.*evmp);
// 		dvmp = sqrt(2./3.*dvmp);
// 
// if((Update_Flag==1)&&(istep>1)){
// 			// grain reorientation
// 			update_orient();
// 		}

		/* Integrated modeling with phase-field */
#ifdef DD_BASED_FLAG
#ifdef PF_DRX
		/* updte the GB_indicator */
		local_loop{
			int chk = (GB_checkX[pIDX][0]+GB_checkX[pIDX][1]+
					GB_checkY[pIDX][0]+GB_checkY[pIDX][1]+
					GB_checkZ[pIDX][0]+GB_checkZ[pIDX][1]);
      /* 09/16/15 --- PYZ */
      /* >0: at GB; =0: in bulk */
      GB_indicator[pIDX] = 6-chk;
			if(chk<6){
				GB_indicator[pIDX] = 1;	// voxel is adjacent to GB
			}
			else{
				GB_indicator[pIDX] = 0;	
			}
		} 
		/* Implement the DRX nucleation based
		   on the dislocation density difference. This
		   will updat the "diff_rho" */
		   
	/*	  if(gID_rex[pIDX] == 1){
		       if(gsl_rng_uniform(RandInstance) >= 0.5) {
		           gID_rex[pIDX] = 0;
		           if(mpirank==0){
		               printf("gID_rex = 0.0");
		           
		       }
		   }
		   }*/

               //  printf("kappa is %le\n",kappa_drx[0]);
    

//cp
/*	sort(gID_list.begin(),gID_list.end(),G_Info::before);
dd_size=gID_list.size();
	double *dd_id = (double*)malloc(dd_size*sizeof(double));
	 dd_index=0;
	 	for(it=gID_list.begin();it!=end;it++){
	 	   id1=gID_list.ID;
	 	   dd=0;
	 	  struct dd_info temp_dd={id1,dd};
	 	   dd_list.pushback(temp_dd);
	 	}*/
	 //	dd_size=gID_list.size();
	 //	double *dd_average = (double*)malloc(dd_size*sizeof(double));
	 //	dd_index=0;
	     
	     
	//std::vector<G_Info>::iterator it2, end2;
	 	//end2=gID_list.end();
	//for(it2=gID_list.begin();it2!=end2;it2++){
	    
		//dd_temp = 0.0;
		//g_count = 0;
		//dd_ave=0.0;
		
		//local_loop{
		 //   if(grain_f[pIDX]==(*it2).ID){
		    //    dd_temp += rho_tot[pIDX];
		    //    g_count ++ ;
                       /* if((*it2).ID==17){
		        	printf("%lf \n",rho_tot[pIDX]);
		      //  d_index[pIDX]=dd_index;

} */		  //  }
		   
		//}
	//printf("dd_average=%lf %d\n",dd_temp,g_count);
				//MPI_Barrier(MPI_COMM_WORLD);
		//MPI_Allreduce(&dd_temp, &dd_ave, 1, MPI_real,
				//MPI_SUM, MPI_COMM_WORLD);
		//MPI_Allreduce(&g_count, &mesh_count, 1, MPI_INT,
				//MPI_SUM, MPI_COMM_WORLD);
				// printf("dd_average=%lf %d\n",dd_ave,mesh_count);
	//dd_average[(*it2).ID - 1]= real (dd_ave/mesh_count);
				// printf("dd_average=%lf %le\n",dd_ave,dd_average[(*it2).ID -1]);
			//	 exit(0);
			//	index= gID_list.ID;
				/*	for(it2=dd_list.begin();it2!=end;it2++){
					    if(index==(*it2).ID){
					        (*it2).avg_dd= dd_id;
					       
					    }
					}*/
				
					
					/*local_loop{
				
		    if(grain_f[pIDX]==(*it2).ID){
		        dd_average[(*it2.ID] =dd_ave;
		           
			//		dd_index++;
	}
						} */
					
						
	//}

local_loop{
    if(first_pf[pIDX] == 1) {
            jph = phase_f[pIDX];
           rho_tot[pIDX] = 0.0;
                    for(i=0;i<nSYS[jph-1];i++){
                    rho_tot[pIDX] += 2.0 + 3.0 + 3.0 ; // + pow(ssdScaler_DRX,ntimes)*rho_s[pIDX][i]; // pow(gndScaler_DRX,ntimes)*rho_m[pIDX][i]
			//	+  pow(gndScaler_DRX,ntimes)*sqrt(pow(rho_g1[pIDX][i],2.0)+pow(rho_g2[pIDX][i],2.0)+pow(rho_g3[pIDX][i],2.0));
		}
        } else {
        
         jph = phase_f[pIDX];
           rho_tot[pIDX] = 0.0;
                    for(i=0;i<nSYS[jph-1];i++){
                    rho_tot[pIDX] += rho_s[pIDX][i]+ rho_m[pIDX][i]
				+  sqrt(pow(rho_g1[pIDX][i],2.0)+pow(rho_g2[pIDX][i],2.0)+pow(rho_g3[pIDX][i],2.0));
		}   
            
        }
}


#ifdef DRX_ELASTIC_TREAT
			/* relax stress after DRX
			   We treat DRX as an instantaneous phase transformation and solve
			   the resulted elasticity issue*/
//			WriteSigMPI("sig_before",0);
//			WriteDisgradMPI("disg_before",0);
//			InhomogeneousStressSolver(C_gr, Eps, DisGradAvg_actual,
//					C0, Eps0_last,
//					Tau_0, Tau,
//					DisGrad, Sig);
			ten2nd epsavg;
			T2_loop{
				epsavg[mi][mj] = 0.0;
				T2p_loop{
					epsavg[mi][mj] = S0[mi][mj][mip][mjp]*SigAvg[mip][mjp];
				}
			}
//			if(mpirank==0){
//				printf("DisGradAvg:\n");
//				PrintTensor(DisGradAvg_actual);
//				printf("EpsAvg:\n");
//				PrintTensor(EpsAvg);
//				printf("epsavg:\n");
//				PrintTensor(epsavg);
//			}
			ElasticEVP(epsavg);
			//ElasticEVP_2(epsavg, Eps);
//			WriteSigMPI("sig_after",0);
//			WriteDisgradMPI("disg_after",0);
//			PError("pyz",1203);
#endif
			/* relax stress after DRX */

			// keep current DisGrad
	
			/* "dDisGradAvg" is the local strain deviation part from
			   the average, i.e. "E-<e>". It is calculated in get_smacro()
			   according to Eq. 23 (the second term on the l.h.s). It is
			   calculated for each mechanical step iteration, and "dDisGradAvg_acum"
			   is the accumulated value after all the iterations (the following
			   while loop) for the current mechanical step (istep). */
     
// 	 if(growth[pIDX]==1 && rho_tot[pIDX] > 50.0 ){
//               int jph = phase_f[pIDX];
// 			for(int i=0;i<nSYS[jph-1];i++){
// 				rho_s[pIDX][i] = rho_s[pIDX][i]*ssdScaler_DRX;
// 				trial_rho_s[pIDX][i] = rho_s[pIDX][i]*ssdScaler_DRX;
// 				rho_g1[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 				rho_g2[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 				rho_g3[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 				trial_gam_acum[pIDX][i] *=  gndScaler_DRX;	// accumulated shear on each slip system, trial version
// 				trial_rho_g1[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 				trial_rho_g2[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 				trial_rho_g3[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 			}
// 				T2_loop{
// 				/* we treat the DRX as a phase-transformation that
// 				   leads to a plastic-strain-free state */
// 				   Sig[pIDX][mi][mj] *= sigScaler_DRX;
//                                   
//      
// 			}
//             }

                               local_loop{
			T2_loop{
				// store displacement gradient of current time into velgrad
				VelGrad[pIDX][mi][mj] = DisGrad[pIDX][mi][mj];
				// enforce macro deformation
				if(CREEP_FLAG==1){
					DisGrad[pIDX][mi][mj] = DisGradAvg_actual[mi][mj];
				}
				else{
					DisGrad[pIDX][mi][mj] += Udot[mi][mj]*TimeStep;
				}

			}
		}

		if(istep==1 || Update_Flag==1){	// Note the logic operator is "OR"
			update_schmid();
		}
			
			T2_loop{
				dDisGradAvg[mi][mj] = 0.0;
				dDisGradAvg_acum[mi][mj] = 0.0;
			}
				/* re-calculate the stress field */
			iter = 0;
			Err_e = 2.*Err;
			Err_s = 2.*Err;

			/* initial guess */
			local_loop{
				T2_loop{
					Sig[pIDX][mi][mj] *= 1.0;
				}
			}

			while((iter<IterMax)&& (MAX(fabs(Err_s),fabs(Err_e))>Err)){
				iter++;
				if(mpirank==0){
				//	printf("\nITER = %d\n",iter);
				//	printf("Forward FFT of stress field\n");
				}
				// k-space stress field
				T2_loop{
					local_loop{
						//fft_data[pIDX].re = Sig[pIDX][mi][mj];
						//fft_data[pIDX].im = 0.0;
						  fft_data[pIDX][0] = Sig[pIDX][mi][mj];
					fft_data[pIDX][1] = 0.0;
					}
					MPI_Barrier(MPI_COMM_WORLD);
					//fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
fftw_mpi_execute_dft(plan,fft_data,fft_fourier);					
local_loop{
						//kSig_r[pIDX][mi][mj] = fft_data[pIDX].re;
						//kSig_i[pIDX][mi][mj] = fft_data[pIDX].im;
						 kSig_r[pIDX][mi][mj] = fft_fourier[pIDX][0];
					kSig_i[pIDX][mi][mj] = fft_fourier[pIDX][1];
					}
				}
				local_loop{
					T2_loop{
						sym_du_r[mi][mj] = 0.0;;
						sym_du_i[mi][mj] = 0.0;;
						T2p_loop{
							sym_du_r[mi][mj] += GAMMA[pIDX][mi][mj][mip][mjp]*kSig_r[pIDX][mip][mjp];
							sym_du_i[mi][mj] += GAMMA[pIDX][mi][mj][mip][mjp]*kSig_i[pIDX][mip][mjp];
						}
					}
					T2_loop{
						kSig_r[pIDX][mi][mj] = sym_du_r[mi][mj];
						kSig_i[pIDX][mi][mj] = sym_du_i[mi][mj];
					}
	
				}
	
				if(mpirank==0){
				//	printf("Inverse FFT to get strain field\n");
				}
				// update strain field in real space
				T2_loop{
					local_loop{
						//fft_data[pIDX].re = kSig_r[pIDX][mi][mj];
					//fft_data[pIDX].im = kSig_i[pIDX][mi][mj];
                    fft_fourier[pIDX][0] = kSig_r[pIDX][mi][mj];
					fft_fourier[pIDX][1] = kSig_i[pIDX][mi][mj];
					}
					MPI_Barrier(MPI_COMM_WORLD);
					//fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
                fftw_mpi_execute_dft(iplan,fft_fourier,fft_data);
					local_loop{
						//DisGrad[pIDX][mi][mj] += dDisGradAvg[mi][mj] + fft_data[pIDX].re/Nxyz;	// Eq. 15
                    DisGrad[pIDX][mi][mj] += dDisGradAvg[mi][mj] + fft_data[pIDX][0]/Nxyz;	// Eq. 15
						fluctuation[pIDX][mi][mj] +=  fft_data[pIDX][0]/Nxyz;
					}
				
				}
				//	printf("fluctuation is %le\n",fluctuation[0][1][1]);
					
                
                
				
	
				if(mpirank==0){
				//	printf("Update stress field\n");
				}
				Err_e_local = 0.0;
				Err_s_local = 0.0;
	//#ifdef DD_BASED_FLAG
	//#ifdef DD_GND
	//			// graidents of shear rate and shear
	//			Gradient_ShearRate();
	//			Gradient_Shear();
	//#endif
	//#endif
				// update stress, which requires Newton-Raphson method
				// NOTE: N-R is run locally. So PEs could run for different iteration steps
				//UpdateStress_DRX(istep, &Err_e_local, &Err_s_local, 0);	
				UpdateStress(istep, &Err_e_local, &Err_s_local, 1);
//exit(0);	
				// collect errors
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Allreduce(&Err_e_local, &Err_e, 1, MPI_real,
						MPI_SUM, MPI_COMM_WORLD);
				MPI_Allreduce(&Err_s_local, &Err_s, 1, MPI_real,
						MPI_SUM, MPI_COMM_WORLD);
				
				Err_e *= WGT;
				Err_s *= WGT;
				
			//	Err_e = sqrt(Err_e);
			//	Err_s = sqrt(Err_s);
	
				if(mpirank==0){
				//	printf("ERRE = %e\n",Err_e);
				//	printf("ERRS = %e\n",Err_s);
				}
	
				// update the average quantities.
				get_smacro();
	            if(fabs(e_vm)<1E-5){
			  Err_e = 0.0;	
      } else {
          Err_e /= e_vm;
      }
			//	Err_e /= e_vm;	
				Err_s /= s_vm;
				if(mpirank==0){
					printf("Strain field error = %f\n", Err_e);
					printf("Stress field error = %e\n", Err_s);
				//	fprintf(fp_err,"%d\t%e\t%e\t%e\n",iter,Err_e,Err_s,s_vm);
				}
				fflush(fp_err);
	
			}/* while() loop ends and stress converges */
			
			if(mpirank==0){
			//	printf("Stress is relaxed after DRX\n");
			}
			
	/* update the velocity gradient field. VelGrad on the r.h.s.
		   stores the previous displ. gradient fiels */
		local_loop{
			T2_loop{
				VelGrad[pIDX][mi][mj] = (DisGrad[pIDX][mi][mj] - VelGrad[pIDX][mi][mj])/TimeStep;
			}
		}
		T2_loop{
			// dDisGradAvg_acum is updated in get_smacro()
			/* DisGradAvg is always the initially applied strain based on linear ealstic assumption
			   for each mechanical step (istep), and dDisGradAvg_acum is the resulted adjustement
			   after the iteration, which depends on the boundary conditon. */
			DisGradAvg_actual[mi][mj] = DisGradAvg[mi][mj] + dDisGradAvg_acum[mi][mj];
			VelGradAvg[mi][mj] = (DisGradAvg_actual[mi][mj]-DisGradAvg_t[mi][mj])/TimeStep;
		}

		if(mpirank==0){
			printf("VelGrad(1,1),velGrad(2,2),velGrad(3,3)\n");
			printf("%e,%e,%e\n",VelGrad[0][0][0], VelGrad[0][1][1], VelGrad[0][2][2]);
			printf("velGrad(1,2),velGrad(1,3),velGrad(2,1),velgrad(2,3),velgrad(3,1),velgrad(3,2)\n");
			printf("%e %e %e %e %e %e\n",VelGrad[0][0][1], VelGrad[0][0][2], VelGrad[0][1][0],VelGrad[0][1][2],VelGrad[0][2][0],VelGrad[0][2][1]);
			
		}
		e_vm = VonMises(DisGradAvg_actual);
		d_vm = VonMises(VelGradAvg);
		TimeTot += TimeStep;
		T2_loop{
			DisGradAvg_t[mi][mj] = DisGradAvg_actual[mi][mj];
		}
     
		// update strain field
		local_loop{
			T2_loop{
				Eps[pIDX][mi][mj] += Edot[pIDX][mi][mj]*TimeStep;	// Edot was updated in UpdateStress()
			}
		}


		

    	// Plastic VM
		T2_loop{
			EpsAvg_local[mi][mj] = 0.0;
			EdotAvg_local[mi][mj] = 0.0;
		}
		local_loop{
			T2_loop{
				EpsAvg_local[mi][mj] += Eps[pIDX][mi][mj]*WGT;
				EdotAvg_local[mi][mj] += Edot[pIDX][mi][mj]*WGT;
			}
		}
		T2_loop{
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allreduce(&EpsAvg_local[mi][mj], &EpsAvg[mi][mj], 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&EdotAvg_local[mi][mj], &EdotAvg[mi][mj], 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
		}

		evmp = 0.0;
		dvmp = 0.0;
		T2_loop{
			evmp += EpsAvg[mi][mj]*EpsAvg[mi][mj];
			dvmp += EdotAvg[mi][mj]*EdotAvg[mi][mj];
		}
		evmp = sqrt(2./3.*evmp);
		dvmp = sqrt(2./3.*dvmp);

		
		
	//	free(dd_average);

		/* record only the grain IDs, without texture info */
		
		if((Update_Flag==1)&&(istep>1)){
 			// grain reorientation
			update_orient();
		}

#endif
#endif

		if((Hard_Flag==1)&&(istep>1)){
#ifdef DD_BASED_FLAG
#ifdef DD_GND
			// graidents of shear rate and shear
			Gradient_ShearRate();
			Gradient_Shear();
		//	L0 *= (1.0+e_vm/GND_ScaleRatio);
			//L0 *= (1.0+e_vm*e_vm/GND_ScaleRatio);
#endif
#endif
}

			/* dislocation evolution */
			DislocationEvolution(istep);
			AvgDislocation();
			if(mpirank==0){
				fprintf(fp_dd,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
					TimeTot,e_vm,rho_s_avg,rho_m_avg,rho_dipole_avg,rho_inplane_avg,rho_forest_avg);
        fflush(fp_dd);
			}
			/*
				for(jph=0;jph<NPHMAX;jph++){
				  //  Selfinter[jph] = Selfinter0[jph]*pow((log(1/(2.56E-10*sqrt(Selfinter0[jph]*rho_F_avg*1E12)))/log(1/(2.56E-10*sqrt(Selfinter0[jph]*rho_SSD_initial*12*1E12)))),2.0);
				   // Coplanar[jph] =  Coplanar0[jph]*pow((log(1/(2.56E-10*sqrt(Coplanar0[jph]*rho_F_avg*1E12)))/log(1/(2.56E-10*sqrt(Coplanar0[jph]*rho_SSD_initial*12*1E12)))),2.0);
			HirthLock[jph] = HirthLock0[jph]*pow((log(1/(2.56E-10*sqrt(HirthLock0[jph]*rho_F_avg*1E12)))/log(1/(2.56E-10*sqrt(HirthLock0[jph]*rho_SSD_initial*12*1E12)))),2.0);
			GlissileJunction[jph] = GlissileJunction0[jph]*pow((log(1/(2.56E-10*sqrt(GlissileJunction0[jph]*rho_F_avg*1E12)))/log(1/(2.56E-10*sqrt(GlissileJunction0[jph]*rho_SSD_initial*12*1E12)))),2.0);
			LomerCottrellLock[jph] = LomerCottrellLock0[jph]*pow((log(1/(2.56E-10*sqrt(LomerCottrellLock0[jph]*rho_F_avg*1E12)))/log(1/(2.56E-10*sqrt(LomerCottrellLock0[jph]*rho_SSD_initial*12*1E12)))),2.0);
			CrossSlip[jph] = CrossSlip0[jph]*pow((log(1/(2.56E-10*sqrt(CrossSlip0[jph]*rho_F_avg*1E12)))/log(1/(2.56E-10*sqrt(CrossSlip0[jph]*rho_SSD_initial*12*1E12)))),2.0);
			
			InteractionStrengthMatrix("FCC",0.0,0.0,
				CrossSlip[jph], GlissileJunction[jph],
				HirthLock[jph], LomerCottrellLock[jph],
				Chi_1[jph]);
				
				InteractionStrengthMatrix("FCC",Selfinter0[jph],Coplanar0[jph],
				0.0, 0.0,
			 0.0, 0.0,
				Chi[jph]);
				}
				*/
			
			
			if(mpirank==0){
			    printf("hirth=%le, lomer = %le, cross = %le\n",HirthLock[0],LomerCottrellLock[0],CrossSlip[0]);
			}
			
//#else
			// hardening
			// update "crss" and "gamacum"
			//harden();
			//harden_anal();
//#endif
	//	}
    for(mi=0;mi<3;mi++){
				    tmp_re[mi]=0.0;
				    tmp_im[mi]=0.0;
				    for(mj=0;mj<3;mj++){
				
					local_loop{
//CHANGES HERE
						fft_data[pIDX][0] = fluctuation[pIDX][mi][mj];
						fft_data[pIDX][1] = 0.0;
					}
					MPI_Barrier(MPI_COMM_WORLD);
					//fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
  fftw_mpi_execute_dft(plan,fft_data,fft_fourier);
				local_loop{
                     tmp[0] = g[pIDX].x;
                      tmp[1] = g[pIDX].y;
                      tmp[2] = g[pIDX].z;  
                     tmp1 = fft_data[pIDX][0];
           tmp_re[mi] += (-1)*fft_data[pIDX][1]*tmp[mj];
           tmp_im[mi] += tmp1*tmp[mj];
}
}
            fft_data[pIDX][0] = tmp_re[mi];
            fft_data[pIDX][1] = tmp_im[mi];
			MPI_Barrier(MPI_COMM_WORLD);
           // fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
 fftw_mpi_execute_dft(iplan,fft_fourier,fft_data);
            local_loop{
            displacement_fluct[pIDX][mi] = - fft_data[pIDX][0]/Nxyz;
				}
							
				}
		//			printf("displacement is %le\n",displacement_fluct[0][0]);	
    
	 local_loop{
            new_position[pIDX][0] += (new_position[pIDX][0]*VelGrad[pIDX][0][0] + new_position[pIDX][1]*VelGrad[pIDX][0][1] + new_position[pIDX][2]*VelGrad[pIDX][0][2])*TimeStep + displacement_fluct[pIDX][0];
		    new_position[pIDX][1] += (new_position[pIDX][0]*VelGrad[pIDX][1][0] + new_position[pIDX][1]*VelGrad[pIDX][1][1] + new_position[pIDX][2]*VelGrad[pIDX][1][2])*TimeStep + displacement_fluct[pIDX][1];
		    new_position[pIDX][2] += (new_position[pIDX][0]*VelGrad[pIDX][2][0] + new_position[pIDX][1]*VelGrad[pIDX][2][1] + new_position[pIDX][2]*VelGrad[pIDX][2][2])*TimeStep + displacement_fluct[pIDX][2];
    } 
    
   /* local_loop{
        new_position[pIDX][0] = (1.0 + DisGrad[pIDX][0][0])*(px+lxs -32) + DisGrad[pIDX][0][1]*(py-32) + DisGrad[pIDX][0][2]*(pz-32); 
        new_position[pIDX][1] = (1.0 + DisGrad[pIDX][1][1])*(py -32) + DisGrad[pIDX][1][0]*(px+lxs-32) + DisGrad[pIDX][1][2]*(pz-32); 
        new_position[pIDX][2] = (1.0 + DisGrad[pIDX][2][2])*(pz -32) + DisGrad[pIDX][2][0]*(px+lxs-32) + DisGrad[pIDX][2][1]*(py-32); 
    } */

		if(mpirank==0){
			fprintf(fp_vm,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
				TimeTot,e_vm,evmp,d_vm,dvmp,s_vm,s_vm1);
      fflush(fp_vm);
		}
    
		// record stress/strain/strain rate fields
	//	if((PrintControl[0]==1) && (istep>2000) && (count%40 == 0) && (istep<3000)){
	//	    WriteEdotMPI("edot",istep);
	//	    WriteElsMPI("els",istep);
	//	    WriteSigMPI("sig",istep);
	//	} 
		if((PrintControl[0]==1) && (step_drx == 0) && (istep%PrintControl[1] == 0)){
			WriteEpsMPI("eps",istep); // plastic strain field
			WriteElsMPI("els",istep); // elastic strain field
			
		//	WriteEdotMPI("edot",istep);
		//  WriteTextureMPI("tex", istep);
                   /*        WriteEpsMPI("eps_beforeDRX",istep); // plastic strain field
			    WriteSigMPI("sig_beforeDRX",istep);
			    WriteSigMPI("edot_beforeDRX",istep);
		      WriteRhoMPI("rho_beforeDRX","SSD",istep);
		      WriteRhoMPI("rho_beforeDRX","Mobile",istep);
		      WriteRhoMPI("rho_beforeDRX","GND",istep);
			    WriteRhoDotMPI("rhodot_beforeDRX",istep);
		      WriteTextureMPI("tex_beforeDRX", istep);*/
                      //    WriteEpsMPI("eps_afterDRX",istep); // plastic strain field
			    WriteSigMPI("sig_afterDRX",istep);
			 //   WriteSigMPI("edot_afterDRX",istep);
		     WriteRhoMPI("rho_afterDRX","SSD",istep);
		    //  WriteRhoMPI("rho_afterDRX","Mobile",istep);
		    //  WriteRhoMPI("rho_afterDRX","GND",istep);
		     // WriteRhoMPI("rho_afterDRX","ddw",istep);
		    //  WriteRhoMPI("rho_afterDRX","cell",istep);
			  //  WriteRhoDotMPI("rhodot_afterDRX",istep);
		      WriteTextureMPI("tex_afterDRX", istep);
		      WriteSLIPMPI("slip_afterDRX",istep);
		      WriteNewPositionMPI("new_position_afterDRX",istep);
		       WriteDDMPI("dislocation_afterDRX",istep);  
		}


    
#ifdef PF_DRX
if((PrintControl[0]==1)&&(istep%PrintControl[1]==0)){
   // WriteGrainPFMPI("gID_rex",istep);
}
    /* 9/16/15 --- PYZ */
 /*   int RC_flag = 0;
    for(int iix = 0; iix<RC_length; ++iix){
      if(istep==RC_LIST[iix]){
        RC_flag=1;
          break;
      }
    }
    if(RC_flag){
			 //   WriteEpsMPI("eps_afterDRX",istep); // plastic strain field
		//	    WriteSigMPI("sig_afterDRX",istep);
		//	    WriteSigMPI("edot_afterDRX",istep);
		      WriteRhoMPI("rho_afterDRX","SSD",istep);
		      WriteRhoMPI("rho_afterDRX","Mobile",istep);
		      WriteRhoMPI("rho_afterDRX","GND",istep);
		//	    WriteRhoDotMPI("rhodot_afterDRX",istep);
		      WriteTextureMPI("tex_afterDRX", istep);
    }*/
#endif


		// Initial guess of macro disp. gradient at t+dt, always elastic
		T2_loop{
			DisGradAvg[mi][mj] = DisGradAvg_t[mi][mj] + Udot[mi][mj]*TimeStep;
		}
		
		
//	exit(0);

	}	// end of mechanical test


T2_loop{
 	local_loop{
		
		fft_data[pIDX][1] = 0.0;
		fft_data[pIDX][0] = Sig[pIDX][mi][mj]; //sig here??
			tmp_smooth[pIDX] = (pow(g[pIDX].x,2)+pow(g[pIDX].y,2)+pow(g[pIDX].z,2))*s+1.0; //s not declared??
	}
	MPI_Barrier(MPI_COMM_WORLD);
	fftw_mpi_execute_dft(plan,fft_data,fft_fourier);
		local_loop{
			
			fft_fourier[pIDX][1] = tmp_smooth[pIDX]*fft_fourier[pIDX][1];//tmp_smooth not declared
			fft_fourier[pIDX][0] = tmp_smooth[pIDX]*fft_fourier[pIDX][0];
		}
		MPI_Barrier(MPI_COMM_WORLD);
	fftw_mpi_execute_dft(iplan,fft_data,fft_fourier);
		local_loop{
			sig_smooth[pIDX][mi][mj] = fft_data[pIDX][0]/Nxyz;
		}
}

	return;
}/*end Evolution()*/
