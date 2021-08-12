#include "evp.h"

#define TINY (1.0E-20)

void VoigtToFull(real *v, ten2nd m)
{
	m[0][0] = v[0];	m[0][1] = v[5];	m[0][2] = v[4];
	m[1][0] = v[5];	m[1][1] = v[1];	m[1][2] = v[3];
	m[2][0] = v[4];	m[2][1] = v[3];	m[2][2] = v[2];
	
	return;
}/*end VoigtToFull()*/

void chg_basis(voigt V6, ten2nd T2, voigt66 C2, ten4th T4, int opt)
{
	/* The V6 notation here adopts the notation of Lebensohn and Tome, 93
	   opt = 0: Define the global ten2nd B[6]
	   opt = 1: V6 -> T2
	   opt = 2: T2 -> V6
	   opt = 3: C2 -> T4 
	   opt = 4: T4 -> C2 (Cijkl -> Cij) */
	real RSQ2 = 0.70710678118654744;
    real RSQ3 = 0.57735026918962584;
    real RSQ6 = 0.40824829046386304;

	int i, j;

	if(opt==0){
		for(i=0;i<6;i++){
			T2_loop{
				B(mi+1,mj+1,i+1) = 0.0;
			}
		}

        B(1,1,2)=-1.*RSQ6;
        B(2,2,2)=-1.*RSQ6;
        B(3,3,2)= 2.*RSQ6;

        B(1,1,1)=-1.*RSQ2;
        B(2,2,1)= RSQ2;

        B(2,3,3)=RSQ2;
        B(3,2,3)=RSQ2;

        B(1,3,4)=RSQ2;
        B(3,1,4)=RSQ2;

        B(1,2,5)=RSQ2;
        B(2,1,5)=RSQ2;

        B(1,1,6)=RSQ3;
        B(2,2,6)=RSQ3;
        B(3,3,6)=RSQ3;
	}
	else if(opt==1){
		T2_loop{
			T2[mi][mj] = 0.0;
			for(i=0;i<6;i++)
				T2[mi][mj] += V6[i]*B[i][mi][mj];
		}
	}
	else if(opt==2){
		for(i=0;i<6;i++){
			V6[i] = 0.0;
			T2_loop{
				V6[i] += T2[mi][mj]*B[i][mi][mj];
			}
		}
	}
	else if(opt==3){
		C4_loop{
			T4[mi][mj][mk][ml] = 0.0;
			for(i=0;i<6;i++)
				for(j=0;j<6;j++)
					T4[mi][mj][mk][ml] += C2[i][j]*B[i][mi][mj]*B[j][mk][ml];
		}
	}
	else if(opt==4){
		for(i=0;i<6;i++){
			for(j=0;j<6;j++){
				C2[i][j] = 0.0;
				C4_loop{
					C2[i][j] += T4[mi][mj][mk][ml]*B[i][mi][mj]*B[j][mk][ml];
				}
			}
		}
	}
	else{
		PError("chg_basis: Wrong option!!",4);
	}

	return;
}/*end chg_basis()*/

void chg_basis5(voigt5 V6, ten2nd T2, voigt66 C2, ten4th T4, int opt)
{
	/* The V6 notation here adopts the notation of Lebensohn and Tome, 93
	   opt = 0: Define the global ten2nd B[6]
	   opt = 1: V6 -> T2
	   opt = 2: T2 -> V6
	   opt = 3: C2 -> T4 
	   opt = 4: T4 -> C2 (Cijkl -> Cij) */
	real RSQ2 = 0.70710678118654744;
    real RSQ3 = 0.57735026918962584;
    real RSQ6 = 0.40824829046386304;

	int i, j;

	if(opt==0){
		for(i=0;i<6;i++){
			T2_loop{
				B(mi+1,mj+1,i+1) = 0.;
			}
		}

        B(1,1,2)=-RSQ6;
        B(2,2,2)=-RSQ6;
        B(3,3,2)= 2.*RSQ6;

        B(1,1,1)=-RSQ2;
        B(2,2,1)= RSQ2;

        B(2,3,3)=RSQ2;
        B(3,2,3)=RSQ2;

        B(1,3,4)=RSQ2;
        B(3,1,4)=RSQ2;

        B(1,2,5)=RSQ2;
        B(2,1,5)=RSQ2;

        B(1,1,6)=RSQ3;
        B(2,2,6)=RSQ3;
        B(3,3,6)=RSQ3;
	}
	else if(opt==1){
		T2_loop{
			T2[mi][mj] = 0.0;
			for(i=0;i<5;i++)
				T2[mi][mj] += V6[i]*B[i][mi][mj];
		}
	}
	else if(opt==2){
		for(i=0;i<5;i++){
			V6[i] = 0.0;
			T2_loop{
				V6[i] += T2[mi][mj]*B[i][mi][mj];
			}
		}
	}
	else if(opt==3){
		C4_loop{
			T4[mi][mj][mk][ml] = 0.0;
			for(i=0;i<5;i++)
				for(j=0;j<5;j++)
					T4[mi][mj][mk][ml] += C2[i][j]*B[i][mi][mj]*B[j][mk][ml];
		}
	}
	else if(opt==4){
		for(i=0;i<5;i++){
			for(j=0;j<5;j++){
				C2[i][j] = 0.0;
				C4_loop{
					C2[i][j] += T4[mi][mj][mk][ml]*B[i][mi][mj]*B[j][mk][ml];
				}
			}
		}
	}
	else{
		PError("chg_basis: Wrong option!!",4);
	}

	return;
}/*end chg_basis5()*/

void SymAntDecompose(ten2nd t, ten2nd s, ten2nd a)
{
	T2_loop{
		s[mi][mj] = (t[mi][mj]+t[mj][mi])/2.;
		a[mi][mj] = (t[mi][mj]-t[mj][mi])/2.;
	}
	
	return;
}/*end SymAntDecompose()*/


void EulerToTransMatrix(real *p_t1, real *p_ph, real *p_t2, ten2nd a, int opt)
{
	/* Using Bunge angles (rotate from sample to xtal frames)
	   opt = 1: matrix -> angles
	   opt = 2:	angles -> matrix
	   *angles are in radian */
	real ct1,st1,cph,sph,ct2,st2;
	real t1, ph, t2;

	t1 = *p_t1; ph = *p_ph; t2 = *p_t2;

	if(opt==1){
		ph = acos(a[2][2]);
		if(fabs(a[2][2])>0.9999){
			t2 = 0.;	// For this special case (ph=0), the rotation is in-plane
						// and one need to put a constraint, e.g., t1=t2. Here we
						// set t2 = 0.
			t1 = atan2(a[0][1],a[0][0]);
		}
		else{
			sph = sin(ph);
			t1 = atan2(a[2][0]/sph, -1.0*a[2][1]/sph);
			t2 = atan2(a[0][2]/sph, a[1][2]/sph);
		}
		ph *= 180./PI;
		t1 *= 180./PI;
		t2 *= 180./PI;
		*p_ph = ph;
		*p_t1 = t1;
		*p_t2 = t2;
	}
	else if(opt==2){
		ct1 = cos(t1);
		st1 = sin(t1);
		cph = cos(ph);
		sph = sin(ph);
		ct2 = cos(t2);
		st2 = sin(t2);
		a[0][0] = ct1*ct2-st1*st2*cph;
		a[0][1] = st1*ct2+ct1*st2*cph;
		a[0][2] = st2*sph;
		a[1][0] = -1.0*ct1*st2-st1*ct2*cph;
		a[1][1] = -1.0*st1*st2+ct1*ct2*cph;
		a[1][2] = ct2*sph;
		a[2][0] = st1*sph;
		a[2][1] = -1.0*ct1*sph;
		a[2][2] = cph;
	}
	else{
		PError("Wrong option for EulerToTransMatrix()!!",604);
	}

	return;
}/*end EulerToTransMatrix()*/


void Ten4thTransform(ten4th c1, ten2nd a, ten4th c2, int opt)
{	
	/* assume a is for 1->2 */
	int i,j,k,l;
	int m,n,p,q;

	if(opt==1){		// Given c1, determine c2
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
					for(l=0;l<3;l++){
						c2[i][j][k][l] = 0.0;
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								for(p=0;p<3;p++)
									for(q=0;q<3;q++)
										c2[i][j][k][l] += a[i][m]*a[j][n]*a[k][p]*a[l][q]*c1[m][n][p][q];
					}
	}
	else if(opt==2){	// Given c2, deterine c1
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
					for(l=0;l<3;l++){
						c1[i][j][k][l] = 0.0;
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								for(p=0;p<3;p++)
									for(q=0;q<3;q++)
										c1[i][j][k][l] += a[m][i]*a[n][j]*a[p][k]*a[q][l]*c2[m][n][p][q];
					}
	}
	else{
		PError("Wrong option for Ten4thTransfrom()!!",445);
	}

	return;
}/*end Ten4thTransform()*/

static void LU_dcmp(real **a, int n, int *indx, real *d)
{
	/* Adopted from "Numerical Recipeis in C" by W Press et. al.

	   Given a matrix a[1..n][1..n], this function REPLACE it by
	   the LU decomposition of a rowwise permutation of itself.
INPUT:
	a -- original matrix
	n -- dimension
OUTPUT:
	a -- matrix containing L and U
	indx -- row permuation effected by the partial pivoting
	d -- +1/-1, depending on the number of row interchanges was even or odd, respectively. */

	int i,imax,j,k;
	real big,dum,sum,temp;
	real *vv; //vv stores the implicit scaling of each row.

	vv=(real*)malloc(n*sizeof(real));
	*d=1.0;		//No row interchanges yet.
	for (i=1;i<=n;i++) { //Loop over rows to get the implicit scaling information
		big=0.0; 
		for (j=1;j<=n;j++)
		if ((temp=fabs(a[i-1][j-1])) > big) big=temp;
		if (fabs(big)<EPSILON_FFT){
			printf("Singular matrix in function LU_dcmp()!!\n");
			for(j=1;j<=n;j++){
					printf("pyz:rank#%d:LU_dcmp: a[%d][%d]=%lf\n",mpirank,i,j,a[i-1][j-1]);
					fflush(stdout);
			}
			exit(35);
		}
		//No nonzero largest element.
		vv[i-1]=1.0/big; //Save the scaling.
	}
	for (j=1;j<=n;j++) {	//This is the loop over columns of Crout's method.
		for (i=1;i<j;i++) {		//This is equation (2.3.12) except for i = j.
			sum=a[i-1][j-1];
			for (k=1;k<i;k++) sum -= a[i-1][k-1]*a[k-1][j-1];
			a[i-1][j-1]=sum;
		}
		big=0.0;	//Initialize for the search for largest pivot element.
		for (i=j;i<=n;i++) {	//This is i = j of equation (2.3.12) and i = j+1...N
			sum=a[i-1][j-1];		//of equation (2.3.13).
			for (k=1;k<j;k++) sum -= a[i-1][k-1]*a[k-1][j-1];
			a[i-1][j-1]=sum;
			if ( (dum=vv[i-1]*fabs(sum)) >= big) {
				// Is the figure of merit for the pivot better than the best so far?
				big=dum;
				imax=i;
			}
		}
		//Do we need to interchange rows?
		if (j != imax) {	//Yes, do so...
			for (k=1;k<=n;k++) {	
				dum=a[imax-1][k-1];
				a[imax-1][k-1]=a[j-1][k-1];
				a[j-1][k-1]=dum;
			}
			*d = -(*d);		//...and change the parity of d.
			vv[imax-1]=vv[j-1];	//Also interchange the scale factor.
		}
		indx[j-1]=imax;
		if (fabs(a[j-1][j-1]) < 1E-5) a[j-1][j-1]=TINY;
		/* If the pivot element is zero the matrix is singular (at least to the precision of the
		algorithm). For some applications on singular matrices, it is desirable to substitute
		TINY for zero.*/
		if (j != n) {	//Now, finally, divide by the pivot element.
			dum=1.0/(a[j-1][j-1]);
			for (i=j+1;i<=n;i++) a[i-1][j-1] *= dum;
		}
	}	//Go back for the next column in the reduction.
	
	free(vv);

	return;
}/*end LU_dcmp()*/

static void LU_bksb(real **a, int n, int *indx, real b[])
{
	/* Adopted from "Numerical Recipeis in C" by W Press et. al.

	   Solve the set of n linear equations A*X = B.
INPUT:
	a -- matrix A in the LU decomposition form obtained through LU_dcmp()
	n -- dimension
	indx -- output of LU_dcmp()
	b -- vextor B
OUTPUT:
	b -- solution X.							*/

	int i,ii=0,ip,j;
	real sum;

	/* When ii is set to a positive value, it will become the
		index of the first nonvanishing element of b. We now
		do the forward substitution, equation (2.3.6). The
		only new wrinkle is to unscramble the permutation
		as we go.*/
	for (i=1;i<=n;i++) {
		ip=indx[i-1];
		sum=b[ip-1];
		b[ip-1]=b[i-1];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i-1][j-1]*b[j-1];
		else if (fabs(sum)>1E-5) ii=i;		//A nonzero element was encountered, so from now on we
								//will have to b[i]=sum; do the sums in the loop above.
		b[i-1] = sum;
	}
	for (i=n;i>=1;i--) {	//Now we do the backsubstitution, equation (2.3.7).
		sum=b[i-1];
		for (j=i+1;j<=n;j++) sum -= a[i-1][j-1]*b[j-1];
		b[i-1]=sum/a[i-1][i-1];	//Store a component of the solution vector X.
	}	// All done!
}/*end LU_bksb()*/

void LU_inv_66(voigt66 c)
{
	/* Inverse the matrix using LU decomposition.
	   The original matrix will be replaced with
	   its inverse.
		*A special case for 6x6 matrix with voigt66 type,
		modified from LU_inverse(real**, int) */

	int indx[6];
	int n;
	real d;
	int i, j;
	voigt col;
	voigt66 y;
	real *a[6] = {c[0],c[1],c[2],c[3],c[4],c[5]};

	n = 6;
	LU_dcmp(a,n,indx,&d);
	for(j=0;j<n;j++){
		for(i=0;i<n;i++) col[i] = 0.0;
		col[j] = 1.0;
		LU_bksb(a,n,indx,col);
		for(i=0;i<n;i++) y[i][j] = col[i];
	}

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			c[i][j] = y[i][j];

	return;
}/*end LU_inv_66()*/

void update_schmid(void)
{
	/* rotate Schmid tensor of each grain from crystal to sample axes */
	voigt5 aux5;
	ten2nd aux33,aux33r;
	voigt66 aux66;
	ten4th aux3333;
	int jph,i,j;

	local_loop{
		jph = phase_f[pIDX];

		if(!Type_phases[jph-1]){	// NOT gas!!
			for(i=0;i<nSYS[jph-1];i++){
				for(j=0;j<5;j++) aux5[j] = Schm_xt[jph-1][i][j];
				chg_basis5(aux5,aux33,aux66,aux3333,1);
				// transfrom from xtal to sample(grain)
				T2_loop{
					aux33r[mi][mj] = 0.0;
					T2p_loop{
						aux33r[mi][mj] += TranMat_xt2sa[pIDX][mi][mip]*TranMat_xt2sa[pIDX][mj][mjp]*aux33[mip][mjp];
					}
				}
				chg_basis5(aux5,aux33r,aux66,aux3333,2);
				// update Schmid tensor at the grain point
				for(j=0;j<5;j++){
					Schm_gr[pIDX][i][j] = aux5[j];
				}
			}
		}
	}

	return;
}/*end update_schmid()*/

real VonMises(ten2nd t)
{
	/* Calculate the von Mises equivalent of a non-symmetric,
	   non-traceless (disp. gradent or velocity gradent) tenor */
	real trace;
	ten2nd dt;
	real vm;

	trace = t[0][0]+t[1][1]+t[2][2];
	T2_loop{
		dt[mi][mj] = (t[mi][mj]+t[mj][mi])/2. - (real)(mi==mj)*trace/3.0;
	}
	vm = 0.0;
	T2_loop{
		vm += dt[mi][mj]*dt[mi][mj];
	}

	return (sqrt(2./3.*vm));	// NOTE: for strain, it's 2/3, NOT 3/2(for stress)

}/* end VonMises()*/
