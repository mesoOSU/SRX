#include "evp.h"

#ifdef PF_DRX
void HomogeneousStressSolver(ten4th c0ijkl, ten2nd *eps0, ten2nd epsavg,
		ten2nd *eps, ten2nd *sig)
{
	ten2nd sym_du_r, sym_du_i;

	T2_loop{
		
		// "stress-free" stress defined through SFTS
		local_loop{
			fft_data[pIDX].re = eps0[pIDX][mi][mj];
			fft_data[pIDX].im = 0.0;
		}
		// forward FFT
		MPI_Barrier(MPI_COMM_WORLD);
		fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
		local_loop{
			kSig_r[pIDX][mi][mj] = fft_data[pIDX].re;
			kSig_i[pIDX][mi][mj] = fft_data[pIDX].im;
		}
	}

	local_loop{
		T2_loop{
			sym_du_r[mi][mj] = 0.0;;
			sym_du_i[mi][mj] = 0.0;;
			T2p_loop{
				sym_du_r[mi][mj] += c0ijkl[mi][mj][mip][mjp]*kSig_r[pIDX][mip][mjp];
				sym_du_i[mi][mj] += c0ijkl[mi][mj][mip][mjp]*kSig_i[pIDX][mip][mjp];
			}
		}
		T2_loop{
			kSig_r[pIDX][mi][mj] = sym_du_r[mi][mj];
			kSig_i[pIDX][mi][mj] = sym_du_i[mi][mj];
		}
	}

	T2_loop{
		// strain field fluctuation in k-space

		local_loop{
			fft_data[pIDX].re = 0.0;
			fft_data[pIDX].im = 0.0;
			T2p_loop{
				fft_data[pIDX].re += (GAMMA[pIDX][mi][mj][mip][mjp]+GAMMA[pIDX][mj][mi][mip][mjp])*kSig_r[pIDX][mip][mjp];
				fft_data[pIDX].im += (GAMMA[pIDX][mi][mj][mip][mjp]+GAMMA[pIDX][mj][mi][mip][mjp])*kSig_i[pIDX][mip][mjp];
			}
			fft_data[pIDX].re *= -0.5;
			fft_data[pIDX].im *= -0.5;
		}
		// inverse FFT to get total strain field
		MPI_Barrier(MPI_COMM_WORLD);
		fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
		local_loop{
			eps[pIDX][mi][mj] = fft_data[pIDX].re/Nxyz + epsavg[mi][mj];	// total strain
		}
	}

	// stress field
	local_loop{
		T2_loop{
			sig[pIDX][mi][mj]= 0.0;
			T2p_loop{
				sig[pIDX][mi][mj] += c0ijkl[mi][mj][mip][mjp]*(eps[pIDX][mip][mjp]-eps0[pIDX][mip][mjp]);
			}
		}
	}

		
	return;
}/*end HomogeneousStressSolver()*/

void InhomogeneousStressSolver(voigt66 *cij, ten2nd *eps0, ten2nd epsavg,
		ten4th c0ijkl, ten2nd *eps0_last,
		ten2nd *Tau_0, ten2nd *Tau,
		ten2nd *eps, ten2nd *sig)
{
	int step;
	real convg, convg_local, normepsavg;
	voigt66 s066,c066;
	voigt aux6;
	ten2nd aux33;
	ten4th cc;
	ten4th s0ijkl;

	chg_basis(aux6,aux33,c066,c0ijkl,4);
	C6_loop{
		s066[mi][mj] = c066[mi][mj];
	}
	LU_inv_66(s066);
	chg_basis(aux6,aux33,s066,s0ijkl,3);

	/* The first term of r.h.s. of Eq. (28) for JL note */
	local_loop{
		chg_basis(aux6,aux33,cij[pIDX],cc,3);
		T2_loop{
			Tau_0[pIDX][mi][mj] = 0.0;
			T2p_loop{
				Tau_0[pIDX][mi][mj] += cc[mi][mj][mip][mjp]*eps0[pIDX][mip][mjp];
			}
		}
	}

	/* Initial guess of the virtual strain field */
	local_loop{
		T2_loop{
			eps[pIDX][mi][mj] = eps0[pIDX][mi][mj]+epsavg[mi][mj];
		}
	}

	normepsavg = VonMises(epsavg);
	step = 0;
	while(1){
		step++;
		/* record the reference system */
		local_loop{
			T2_loop{
				eps0_last[pIDX][mi][mj] = eps[pIDX][mi][mj];
			}
		}

		/* solve the reference homogeneous system */
		HomogeneousStressSolver(c0ijkl, eps0_last, epsavg,
			eps, sig);

		/* total stress */
		local_loop{
			chg_basis(aux6,aux33,cij[pIDX],cc,3);
			T2_loop{
				Tau[pIDX][mi][mj] = 0.0;
				T2p_loop{
					/* The 2nd term of r.h.s in Eq. 28 */
					Tau[pIDX][mi][mj] += (c0ijkl[mi][mj][mip][mjp]-cc[mi][mj][mip][mjp])*eps[pIDX][mip][mjp];
				}
				Tau[pIDX][mi][mj] += Tau_0[pIDX][mi][mj];
			}
		}
		
		/* updated virtual strain field */
		local_loop{
			T2_loop{
				eps[pIDX][mi][mj] = 0.0;
				T2p_loop{
					eps[pIDX][mi][mj] += s0ijkl[mi][mj][mip][mjp]*Tau[pIDX][mip][mjp];
				}
			}
		}

		/* check convergence */
		convg_local = 0.0;
		local_loop{
			T2_loop{
				convg_local += pow((eps[pIDX][mi][mj]-eps0_last[pIDX][mi][mj]),2.0);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&convg_local, &convg, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		convg /= ((real)(CellDim[0])*(real)(CellDim[1])*(real)(CellDim[2]));
		convg = sqrt(convg);
		if(mpirank==0){
			printf("Elasticity solver step = %d: convergence = %e\n",step, convg);
		}
		if((fabs(convg) <((1E-3)*normepsavg))||(fabs(convg)<EPSILON_FFT)){
			break;
		}
		if(step>1000){
			break;
		}
	}

	return;
}/*end InhomogeneousStressSolver()*/

static real Norm_t2nd(ten2nd tt)
{
	real norm = 0.0;

	T2_loop{
		norm += tt[mi][mj]*tt[mi][mj];
	}

	return sqrt(norm);
}/* end Norm_t2nd()*/

void ElasticEVP(ten2nd epsavg)
{
	/* Follow Ricardo's elastic-FFT code */
	int iter=0;
	int itermax = 100;
	real err2mod = 2*Err;
	int jph = 0;
	ten4th cc = {0.0};	// store local C_gr
	ten4th FSloc = {0.0}; // store local FS_gr
	ten2nd sigavg = {0.0};
	ten2nd sigavg1 = {0.0};
	ten2nd epstot = {0.0};
	ten2nd e2nd = {0.0};
	ten2nd sigdiff = {0.0};
	ten2nd epsdiff = {0.0};
	ten2nd Xloc = {0.0};
	voigt aux6;
	ten2nd aux33;

	real local_sigavg = 0.0;
	real local_sigavg1 = 0.0;
	real errd_local = 0.0;
	real errs_local = 0.0;
	real errd, errs;
	/* initialize strain field fluctuation field, which
	   is stored in Eps0_last */
	local_loop{
		T2_loop{
			Eps0_last[pIDX][mi][mj] = 0.0;
		}
	}
	/* initialize stress field */
	local_loop{
		jph = phase_f[pIDX];
		chg_basis(aux6,aux33,C_gr[pIDX],cc,3);
		T2_loop{
			Sig[pIDX][mi][mj] = 0.0;
			if(!Type_phases[jph-1]){
				T2p_loop{
					Sig[pIDX][mi][mj] += cc[mi][mj][mip][mjp]*(epsavg[mip][mjp]+
							Eps0_last[pIDX][mip][mjp]);
				}
			}
		}
	}
	T2_loop{
		local_sigavg = 0.0;
		local_sigavg1 = 0.0;
		local_loop{
			local_sigavg += Sig[pIDX][mi][mj];
			if(phase_f[pIDX]==1)
				local_sigavg1 += Sig[pIDX][mi][mj];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&local_sigavg, &sigavg[mi][mj], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&local_sigavg1, &sigavg1[mi][mj], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		sigavg[mi][mj] *= WGT;
		sigavg1[mi][mj] *= WGT/Wgt_ph1;
	}

	while((iter<itermax)&&(err2mod>Err)){
		iter++;
		if(mpirank==0){
			printf("\nITER = %d\n",iter);
			printf("Direct FFT of polarization and Lagrange multiplier fields\n");
		}

		/* stress polarization field, stored in Tau, in real space */
		local_loop{
			T2_loop{
				Tau[pIDX][mi][mj] = Sig[pIDX][mi][mj];	// Sig is updated every while step
				T2p_loop{
					Tau[pIDX][mi][mj] -= C0[mi][mj][mip][mjp]*Eps0_last[pIDX][mip][mjp];
				}
			}
		}

		/* forward FFT to get stress polarization field in k-space */
		T2_loop{
			local_loop{
				fft_data[pIDX].re = Tau[pIDX][mi][mj];
				fft_data[pIDX].im = 0.0;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
			local_loop{
				kSig_r[pIDX][mi][mj] = fft_data[pIDX].re;
				kSig_i[pIDX][mi][mj] = fft_data[pIDX].im;
			}
		}

		/* multiply Green operator to get strain polarization field in k-space */
		T2_loop{
			local_loop{
				fft_data[pIDX].re = 0.0;
				fft_data[pIDX].im = 0.0;
				T2p_loop{
					fft_data[pIDX].re += GAMMA[pIDX][mi][mj][mip][mjp]*kSig_r[pIDX][mip][mjp];
					fft_data[pIDX].im += GAMMA[pIDX][mi][mj][mip][mjp]*kSig_i[pIDX][mip][mjp];
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
			local_loop{
				Eps0_last[pIDX][mi][mj] = fft_data[pIDX].re/Nxyz;	// dis_grad fluctuation field
			}
		}
		local_loop{
			T2_loop{
				aux33[mi][mj] = (Eps0_last[pIDX][mi][mj] + Eps0_last[pIDX][mj][mi])/2.0;
			}
			T2_loop{
				Eps0_last[pIDX][mi][mj] = aux33[mi][mj];	// strain fluctuation field
			}
		}

		/* update Lagrange multiplier */
		errd_local = 0.0;
		errs_local = 0.0;
		local_loop{
			jph = phase_f[pIDX];
			if(!Type_phases[jph-1]){
				T2_loop{
					epstot[mi][mj] = epsavg[mi][mj]+Eps0_last[pIDX][mi][mj];	// total strain field
					DisGrad[pIDX][mi][mj] = epstot[mi][mj];
				}
				T2_loop{
					Xloc[mi][mj] = Sig[pIDX][mi][mj];
					T2p_loop{
						Xloc[mi][mj] += C0[mi][mj][mip][mjp]*epstot[mip][mjp];
					}
				}
				chg_basis(aux6,aux33,FS_gr[pIDX],FSloc,3);
				T2_loop{
					e2nd[mi][mj] = 0.0;
					T2p_loop{
						// The 2nd term on the r.h.s of Eq. 16
						e2nd[mi][mj] += FSloc[mi][mj][mip][mjp]*Xloc[mip][mjp];
					}
				}
				T2_loop{
					epsdiff[mi][mj] = epstot[mi][mj]-e2nd[mi][mj];
					sigdiff[mi][mj] = 0.0;
					T2p_loop{
						sigdiff[mi][mj] += C0[mi][mj][mip][mjp]*(epstot[mip][mjp]-e2nd[mip][mjp]);
					}
					Sig[pIDX][mi][mj] += sigdiff[mi][mj];	// update full stress field
				}
				errd_local += Norm_t2nd(epsdiff);
				errs_local += Norm_t2nd(sigdiff);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&errd_local, &errd, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&errs_local, &errs, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		errd *= WGT;
		errs *= WGT;
		

		/* update stress avg */
		T2_loop{
			local_sigavg = 0.0;
			local_sigavg1 = 0.0;
			local_loop{
				local_sigavg += Sig[pIDX][mi][mj];
				if(phase_f[pIDX]==1)
					local_sigavg1 += Sig[pIDX][mi][mj];
			}
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allreduce(&local_sigavg, &sigavg[mi][mj], 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&local_sigavg1, &sigavg1[mi][mj], 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
			sigavg[mi][mj] *= WGT;
			sigavg1[mi][mj] *= WGT/Wgt_ph1;
		}

		errs /= Norm_t2nd(sigavg);
		errd /= Norm_t2nd(epsavg);
		err2mod = (errs+errd)/2.0;

		if(mpirank==0){
			printf("Strain field error = %f\n", errd);
			printf("Stress field error = %e\n", errs);
		}
	}

	return;
}/*end ElasticEVP()*/

void ElasticEVP_2(ten2nd epsavg0, ten2nd *eps)
{
	/* Follow Ricardo's elastic-FFT code */
	int iter=0;
	int itermax = 100;
	real err2mod = 2*Err;
	int jph = 0;
	ten4th cc = {0.0};	// store local C_gr
	ten4th FSloc = {0.0}; // store local FS_gr
	ten2nd sigavg = {0.0};
	ten2nd sigavg1 = {0.0};
	ten2nd epsavg = {0.0};
	ten2nd epstot = {0.0};
	ten2nd e2nd = {0.0};
	ten2nd sigdiff = {0.0};
	ten2nd epsdiff = {0.0};
	ten2nd Xloc = {0.0};
	voigt aux6;
	ten2nd aux33;

	real local_sigavg = 0.0;
	real local_sigavg1 = 0.0;
	real local_epsavg = 0.0;
	real errd_local = 0.0;
	real errs_local = 0.0;
	real errd, errs;
	real epsnorm, signorm;
	/* initialize strain field fluctuation field, which
	   is stored in Eps0_last */
	local_loop{
		T2_loop{
			Eps0_last[pIDX][mi][mj] = 0.0;
		}
	}
	/* initialize stress field */
	local_loop{
		jph = phase_f[pIDX];
		chg_basis(aux6,aux33,C_gr[pIDX],cc,3);
		T2_loop{
			Sig[pIDX][mi][mj] = 0.0;
			if(!Type_phases[jph-1]){
				T2p_loop{
					Sig[pIDX][mi][mj] += cc[mi][mj][mip][mjp]*(epsavg0[mip][mjp]+
							Eps0_last[pIDX][mip][mjp]);
				}
			}
		}
	}
	T2_loop{
		local_sigavg = 0.0;
		local_sigavg1 = 0.0;
		local_loop{
			local_sigavg += Sig[pIDX][mi][mj];
			if(phase_f[pIDX]==1)
				local_sigavg1 += Sig[pIDX][mi][mj];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&local_sigavg, &sigavg[mi][mj], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&local_sigavg1, &sigavg1[mi][mj], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		sigavg[mi][mj] *= WGT;
		sigavg1[mi][mj] *= WGT/Wgt_ph1;
	}

	while((iter<itermax)&&(err2mod>Err)){
		iter++;
		if(mpirank==0){
			printf("\nITER = %d\n",iter);
			printf("Direct FFT of polarization and Lagrange multiplier fields\n");
		}

		/* stress polarization field, stored in Tau, in real space */
		local_loop{
			T2_loop{
				Tau[pIDX][mi][mj] = Sig[pIDX][mi][mj];	// Sig is updated every while step
				T2p_loop{
					Tau[pIDX][mi][mj] -= C0[mi][mj][mip][mjp]*Eps0_last[pIDX][mip][mjp];
				}
			}
		}

		/* forward FFT to get stress polarization field in k-space */
		T2_loop{
			local_loop{
				fft_data[pIDX].re = Tau[pIDX][mi][mj];
				fft_data[pIDX].im = 0.0;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
			local_loop{
				kSig_r[pIDX][mi][mj] = fft_data[pIDX].re;
				kSig_i[pIDX][mi][mj] = fft_data[pIDX].im;
			}
		}

		/* multiply Green operator to get strain polarization field in k-space */
		T2_loop{
			local_loop{
				fft_data[pIDX].re = 0.0;
				fft_data[pIDX].im = 0.0;
				T2p_loop{
					fft_data[pIDX].re += GAMMA[pIDX][mi][mj][mip][mjp]*kSig_r[pIDX][mip][mjp];
					fft_data[pIDX].im += GAMMA[pIDX][mi][mj][mip][mjp]*kSig_i[pIDX][mip][mjp];
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
			local_loop{
				Eps0_last[pIDX][mi][mj] = fft_data[pIDX].re/Nxyz;	// dis_grad fluctuation field
			}
		}
		local_loop{
			T2_loop{
				aux33[mi][mj] = (Eps0_last[pIDX][mi][mj] + Eps0_last[pIDX][mj][mi])/2.0;
			}
			T2_loop{
				Eps0_last[pIDX][mi][mj] = aux33[mi][mj];	// strain fluctuation field
			}
		}

		/* update Lagrange multiplier */
		errd_local = 0.0;
		errs_local = 0.0;
		local_loop{
			jph = phase_f[pIDX];
			if(!Type_phases[jph-1]){
				T2_loop{
					epstot[mi][mj] = epsavg[mi][mj]+Eps0_last[pIDX][mi][mj];	// total strain field
					DisGrad[pIDX][mi][mj] = epstot[mi][mj];
				}
				T2_loop{
					Xloc[mi][mj] = Sig[pIDX][mi][mj];
					T2p_loop{
						Xloc[mi][mj] += C0[mi][mj][mip][mjp]*epstot[mip][mjp];
					}
				}
				chg_basis(aux6,aux33,FS_gr[pIDX],FSloc,3);
				T2_loop{
					e2nd[mi][mj] = eps[pIDX][mi][mj];
					T2p_loop{
						// The 2nd term on the r.h.s of Eq. 16
						e2nd[mi][mj] += FSloc[mi][mj][mip][mjp]*Xloc[mip][mjp];
					}
				}
				T2_loop{
					epsdiff[mi][mj] = epstot[mi][mj]-e2nd[mi][mj];
					sigdiff[mi][mj] = 0.0;
					T2p_loop{
						sigdiff[mi][mj] += C0[mi][mj][mip][mjp]*(epstot[mip][mjp]-e2nd[mip][mjp]);
					}
					Sig[pIDX][mi][mj] += sigdiff[mi][mj];	// update full stress field
				}
				errd_local += Norm_t2nd(epsdiff);
				errs_local += Norm_t2nd(sigdiff);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&errd_local, &errd, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&errs_local, &errs, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		errd *= WGT;
		errs *= WGT;
		

		/* update stress avg */
		T2_loop{
			local_sigavg = 0.0;
			local_sigavg1 = 0.0;
			local_loop{
				local_sigavg += Sig[pIDX][mi][mj];
				if(phase_f[pIDX]==1)
					local_sigavg1 += Sig[pIDX][mi][mj];
			}
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allreduce(&local_sigavg, &sigavg[mi][mj], 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&local_sigavg1, &sigavg1[mi][mj], 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
			sigavg[mi][mj] *= WGT;
			sigavg1[mi][mj] *= WGT/Wgt_ph1;
		}
		/* update strain avg */
		T2_loop{
			local_epsavg = 0.0;
			local_loop{
				local_epsavg += DisGrad[pIDX][mi][mj];
			}
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allreduce(&local_epsavg, &epsavg[mi][mj], 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
			epsavg[mi][mj] *= WGT;
		}

		epsnorm = Norm_t2nd(epsavg);
		signorm = Norm_t2nd(sigavg);
		errs /= signorm;
		errd /= epsnorm;
		err2mod = (errs+errd)/2.0;

		if(mpirank==0){
			printf("Strain norm = %e, Strain field error = %f\n", epsnorm, errd);
			printf("Stress norm = %e, Stress field error = %e\n", signorm, errs);
		}
	}

	return;
}/*end ElasticEVP_2()*/
#endif
