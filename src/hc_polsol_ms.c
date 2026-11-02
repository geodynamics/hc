#include "hc.h"
/*

MULTIPLE SHOOTING SOLVER FOR THE POLOIDAL PROBLEM

This file is an alternative to the propagator matrix solution in
hc_polsol.c and is selected with -solver 1 (the default since
11/2026; -solver 0 gives the original method).

Written by Claude (Anthropic) in collaboration with Thorsten Becker,
November 2026.


1. The problem

For each spherical harmonic degree l, the poloidal flow problem in a
radially layered, incompressible (or Panasyuk-type compressible)
viscous shell is a linear two point boundary value problem for the
six component state vector

   s(r) = [u1 u2 u3 u4 p1 p2]

with u1, u2 the radial and poloidal velocity, u3, u4 the radial and
tangential stress (u3 in the "u3new" convention of hc_polsol, i.e.
including the rho*phi term), and p1, p2 the gravitational potential
and its radial derivative term (see hc_evppot). Between levels where
viscosity, density or the output requirement change, the exact
solution is given by the analytical propagators P (4 x 4, hc_evalpa,
Richards & Hager 1984) and Q (2 x 2, hc_evppot):

   s(r_{k+1}) = D_k s(r_k),   D_k = blockdiag(P_k, Q_k)

At the levels, the state jumps by the density anomaly sources
(density anomaly times layer thickness, dfact) and, for the
compressible case, by the reference density jumps. The boundary
conditions are: at the CMB, no radial velocity, free slip, and a
regular potential below the CMB with the gravity jump due to CMB
topography; at the surface, no radial velocity, free slip or a
prescribed poloidal velocity, and a regular potential above the
surface with the gravity jump due to surface topography.


2. Why the propagator method fails at high degree

hc_polsol integrates from the CMB upwards: it starts three
independent homogeneous solutions (one for each free parameter at the
CMB) and one inhomogeneous solution with all sources, propagates
them level by level, and then solves a 3 x 3 system for the
coefficients that satisfy the three surface conditions. The
propagators contain terms growing like r^(l+1) and decaying like
r^(-l-2). Across the whole mantle, their ratio is
(R/R_cmb)^(2l+3) ~ 1e17 at l = 31 and ~1e33 at l = 63. Any
component along the growing modes gets amplified by this factor,
so after propagation the three homogeneous solutions are dominated
by the same growing mode and are numerically parallel. The 3 x 3
surface system is then singular to within double precision round
off, and the result is noise (in practice, degrees above ~50-80 are
wrong, with the geoid power increasing instead of decreasing with
l). Higher precision arithmetic only moves this limit (the old
"quad precision" compile option), it does not remove it. The CMB
kludge (-cbckl) changes the problem for high l to make it less
sensitive and is a workaround, not a solution.


3. The multiple shooting approach

Instead of integrating across the whole shell, the state vectors at
all N levels are kept as unknowns,

   X = [s_0, s_1, ..., s_{N-1}]      (6 N unknowns)

and the propagation relations are imposed as linear equations

   s_{k+1} - J_{k+1} D_k s_k = f_{k+1},    k = 0 ... N-2

where J_{k+1} is the (linear) reference density jump operator at
level k+1 and f_{k+1} holds the density anomaly and phase boundary
sources at that level (exactly the operations that hc_polsol applies
to the solution vectors after each propagation step). Together with
three boundary conditions on s_0 and three on s_{N-1} this gives a
square, block bidiagonal system A X = F.

Within one interval, the entries of D_k still contain the growth
factor (r_{k+1}/r_k)^(2l+3), but this is a property of a single
block, and Gaussian elimination with partial pivoting (hc_matrix.c,
hc_band_lu_factor) is backward stable: the computed solution is the
exact solution of a problem perturbed by round off times the local
block norms. Since no information has to be carried across the
shell through products of exponentially growing matrices, the
solution is accurate as long as the boundary value problem itself is
well posed, which it is. To keep the local growth moderate, the
intervals are subdivided (hc_ms_init_levels) such that
(r_{k+1}/r_k)^(2l+3) <= max_growth (-msg, 1e4 by default). Tests show
the results are insensitive to this parameter even without any
subdivision; it is a safety margin. The subdivided levels carry no
sources, no density jumps and no output.

The system matrix A depends on l only (through the propagators and
the boundary conditions), the right hand side F on (l, m, A/B)
through the density coefficients and the plate velocities. A is
therefore assembled and factored once per l (hc_ms_factor) and
reused for all m and A/B coefficients (hc_ms_solve). The bandwidth
of A is 9 for the row ordering used here, so the factorization
costs O(N) per l and the solves O(N) per right hand side, which is
comparable to the propagator method for the degrees it can handle.

The results agree with hc_polsol to ~1e-10 for l <= 31 (this is the
round off level of the propagator method there), do not change when
viscosity layers are split (which hc_polsol's results do at the
1e-9 level already at l = 31), and have been tested up to l = 255
(see test/run_test.sh and test/propagator_stability.sh).


4. Row and column ordering of the banded system

column 6 k + c       component c (0...5) of the state at level k
rows 0 ... ncmb-1    CMB conditions (ncmb = 3, or 4 for l = 1)
rows ncmb + 6 k + c  propagation from level k to k+1, component c
rows ncmb + 6 (N-1) ... 6 N - 1   surface conditions (6 - ncmb rows)

With this ordering, all non zero entries lie within 9 diagonals of
the main diagonal, see HC_MS_KL and HC_MS_KU.

For l = 1, hc_polsol does not impose the surface gravity condition
(the degree one potential is special, the center of mass defines the
reference frame) but sets the amplitude of the core potential
solution to zero. The same is achieved here by imposing p1 = 0 at
the CMB as a fourth CMB condition and dropping the surface gravity
condition, which keeps the system square.


5. Usage from hc_polsol

hc_polsol sets up the coarse level structure (rprops, pvisc, den,
qwrite, rho_zero) as before. Then, per call,

   hc_ms_init_levels  builds the refined level structure when the
                      levels, the viscosity or lmax changed
   hc_ms_factor       per degree l: propagators, matrix, LU
   hc_ms_solve        per (l, m, A/B): right hand side, solve,
                      extract the state at the output layers

and hc_polsol writes the output layers to the pol_sol expansions in
the same format as for the propagator method, so that everything
downstream (geoid, velocities, tractions, output) is unchanged.

*/

#define HC_MS_NV 6		/* number of unknowns per level */
#define HC_MS_KL 9		/* lower bandwidth of the system */
#define HC_MS_KU 9		/* upper bandwidth */

/*

set up the refined level structure from the coarse hc->rprops
levels, has to be called after the propagator level section of
hc_polsol has been executed

coarse interval i runs from rprops[i] to rprops[i+1] with viscosity
pvisc[i], the density source den[i] and the output flag qwrite[i]
belong to level i+1

*/
void hc_ms_init_levels(struct hcs *hc, int lmax, int npb, HC_PREC *rpb,
		       hc_boolean verbose)
{
  struct hc_ms *ms;
  int i,j,nsub,n,k,ninho,jpb;
  HC_PREC growth,ratio,lratio,frac;

  ms = &(hc->ms);
  growth = ms->max_growth;
  if(growth < 10.0){
    HC_ERROR("hc_ms_init_levels","max growth factor has to be at least 10");
  }
  /*
     count the refined levels
  */
  n = 1;
  for(i=0;i < hc->nprops;i++){
    nsub = hc_ms_nsub(hc->rprops[i],hc->rprops[i+1],lmax,growth);
    n += nsub;
  }
  if(n > ms->nalloc){
    ms->nalloc = n;
    hc_vecrealloc(&ms->r,n,"hc_ms_init_levels: r");
    hc_vecrealloc(&ms->visc,n,"hc_ms_init_levels: visc");
    hc_vecrealloc(&ms->rho,n,"hc_ms_init_levels: rho");
    hc_vecrealloc(&ms->den,n,"hc_ms_init_levels: den");
    ms->isrc = (int *)realloc(ms->isrc,sizeof(int)*n);
    ms->ipb  = (int *)realloc(ms->ipb,sizeof(int)*n);
    ms->iout = (int *)realloc(ms->iout,sizeof(int)*n);
    if((!ms->isrc) || (!ms->ipb) || (!ms->iout)){
      HC_MEMERROR("hc_ms_init_levels");
    }
  }
  ms->n = n;
  /*
     fill in
  */
  k = 0;
  ninho = 0;
  jpb = 0;
  ms->r[0] = hc->rprops[0];
  ms->rho[0] = hc->rho_zero[0];
  ms->den[0] = 0.0;
  ms->isrc[0] = -1;
  ms->ipb[0] = -1;
  ms->iout[0] = 0;			/* the CMB is output layer 0 */
  ms->nout = 1;
  for(i=0;i < hc->nprops;i++){
    nsub = hc_ms_nsub(hc->rprops[i],hc->rprops[i+1],lmax,growth);
    ratio = hc->rprops[i+1]/hc->rprops[i];
    lratio = log(ratio);
    for(j=1;j <= nsub;j++){
      /*
	 interval k goes from level k to level k+1
      */
      ms->visc[k] = hc->pvisc[i];
      k++;
      if(j < nsub){
	/*
	   intermediate level, equally spaced in log(r), no jumps,
	   no sources
	*/
	frac = (HC_PREC)j/(HC_PREC)nsub;
	ms->r[k] = hc->rprops[i] * exp(frac*lratio);
	ms->rho[k] = hc->rho_zero[i];
	ms->den[k] = 0.0;
	ms->isrc[k] = -1;
	ms->ipb[k] = -1;
	ms->iout[k] = -1;
      }else{
	/*
	   original level i+1
	*/
	ms->r[k] = hc->rprops[i+1];
	ms->rho[k] = hc->rho_zero[i+1];
	ms->den[k] = hc->den[i];
	ms->isrc[k] = ninho;
	/*
	   phase boundaries, same sequential logic as in hc_polsol
	*/
	if((jpb < npb) && (hc->rprops[i+1] > rpb[jpb] - 0.0001)){
	  ms->ipb[k] = jpb;
	  jpb++;
	}else{
	  ms->ipb[k] = -1;
	}
	if(fabs(hc->den[i]) > HC_EPS_PREC){
	  ninho++;
	}
	if(hc->qwrite[i]){
	  ms->iout[k] = ms->nout;
	  ms->nout++;
	}else{
	  ms->iout[k] = -1;
	}
      }
    }
  }
  if(k != n-1){
    HC_ERROR("hc_ms_init_levels","level count mismatch");
  }
  /* viscosity above the last level, for the jump operator */
  ms->visc[n-1] = ms->visc[n-2];
  /*
     band matrix storage and work arrays
  */
  ms->nsys = HC_MS_NV * n;
  ms->ldab = 2*HC_MS_KL + HC_MS_KU + 1;
  hc_vecrealloc(&ms->ab,ms->ldab * ms->nsys,"hc_ms_init_levels: ab");
  hc_vecrealloc(&ms->rhs,ms->nsys,"hc_ms_init_levels: rhs");
  hc_vecrealloc(&ms->prop,16*(n-1),"hc_ms_init_levels: prop");
  hc_vecrealloc(&ms->ppot,4*(n-1),"hc_ms_init_levels: ppot");
  ms->ipiv = (int *)realloc(ms->ipiv,sizeof(int)*ms->nsys);
  if(!ms->ipiv){
    HC_MEMERROR("hc_ms_init_levels: ipiv");
  }
  ms->init = TRUE;
  ms->lfact = -1;
  if(verbose){
    fprintf(stderr,"hc_ms_init_levels: %i coarse levels refined to %i for lmax %i, max growth %g, system size %i\n",
	    hc->nprops+1,n,lmax,(double)growth,ms->nsys);
  }
}
/*
   number of subintervals for an interval r1 ... r2 such that
   (r2/r1)^(2l+3) <= growth
*/
int hc_ms_nsub(HC_PREC r1, HC_PREC r2, int lmax, HC_PREC growth)
{
  HC_PREC x;
  int nsub;
  x = (2.0*(HC_PREC)lmax + 3.0) * log(r2/r1) / log(growth);
  nsub = (int)ceil(x - 1e-10);
  if(nsub < 1){
    nsub = 1;
  }
  return nsub;
}
/*
   free the level structure
*/
void hc_ms_free(struct hc_ms *ms)
{
  if(ms->nalloc){
    free(ms->r);
    free(ms->visc);
    free(ms->rho);
    free(ms->den);
    free(ms->isrc);
    free(ms->ipb);
    free(ms->iout);
    free(ms->ab);
    free(ms->rhs);
    free(ms->prop);
    free(ms->ppot);
    free(ms->ipiv);
  }
  ms->nalloc = 0;
  ms->init = FALSE;
}
/*

assemble and factor the system matrix for degree l

row ordering: ncmb CMB rows, then 6 rows per interval k (k=0...n-2)
for the levels k+1, then 6-ncmb surface rows. column 6*k+c is
component c of the state at level k. ncmb is 3, except for l = 1
where the potential is handled with an extra condition at the CMB
and one less at the surface, see below.

*/
void hc_ms_factor(struct hcs *hc, int l, hc_boolean free_slip)
{
  struct hc_ms *ms;
  int n,k,c,d,row,col,nsys,kl,ku,ldab,ierr,ncmb;
  HC_PREC el,beta,a[6][6],drho,afac,vsum,*p,*q;

  ms = &(hc->ms);
  n = ms->n;
  nsys = ms->nsys;
  kl = HC_MS_KL;
  ku = HC_MS_KU;
  ldab = ms->ldab;
  el = (HC_PREC)l;
  beta = hc->psp.beta;
  /*
     propagators for all intervals, only depend on l
  */
  for(k=0;k < n-1;k++){
    hc_evalpa(l,ms->r[k],ms->r[k+1],ms->visc[k],(ms->prop+16*k));
    hc_evppot(l,(ms->r[k]/ms->r[k+1]),(ms->ppot+4*k));
  }
  /*
     clear the band storage, including the fill in rows
  */
  for(k=0;k < ldab*nsys;k++){
    ms->ab[k] = 0.0;
  }
  /*
     CMB conditions, level 0

     row 0: u1 = 0
     row 1: u4 = 0 (free slip core)
     row 2: p2 - l p1 - beta r0 (u3 - (rho_mantle - rho_core) p1) = 0

     (the potential below the CMB is regular, p2 = l p1, and jumps
     across the CMB by the gravity effect of the CMB topography, as in
     hc_polsol)

     row 3, only for l = 1: p1 = 0. hc_polsol does not impose the
     gravity condition at the surface for l = 1 and sets the core
     potential to zero instead, which is equivalent to this
  */
  if(l == 1){
    ncmb = 4;
    ms->ab[HC_BAND_INDEX(3,4,kl,ku,ldab)] = 1.0;
  }else{
    ncmb = 3;
  }
  ms->ab[HC_BAND_INDEX(0,0,kl,ku,ldab)] = 1.0;
  ms->ab[HC_BAND_INDEX(1,3,kl,ku,ldab)] = 1.0;
  ms->ab[HC_BAND_INDEX(2,5,kl,ku,ldab)] = 1.0;
  ms->ab[HC_BAND_INDEX(2,4,kl,ku,ldab)] = -el + beta * ms->r[0] *
    (ms->rho[0] - hc->rho_zero[-1]);
  ms->ab[HC_BAND_INDEX(2,2,kl,ku,ldab)] = -beta * ms->r[0];
  /*
     propagation rows
  */
  for(k=0;k < n-1;k++){
    p = ms->prop + 16*k;
    q = ms->ppot + 4*k;
    /*
       A = J D, with D = blockdiag(P,Q)
    */
    for(c=0;c < 6;c++){
      for(d=0;d < 6;d++){
	a[c][d] = 0.0;
      }
    }
    for(c=0;c < 4;c++){
      for(d=0;d < 4;d++){
	a[c][d] = p[c*4+d];
      }
    }
    a[4][4] = q[0];
    a[4][5] = q[1];
    a[5][4] = q[2];
    a[5][5] = q[3];
    /*
       density jump at level k+1, as in hc_polsol:
       u1 += u1 drho/rho, u3 -= 2 (v_k + v_{k+1}) u1 drho/rho + drho p1,
       u4 += (v_k + v_{k+1}) u1 drho/rho
    */
    drho = ms->rho[k] - ms->rho[k+1];
    if(fabs(drho) > 0.0){
      afac = drho/ms->rho[k+1];
      vsum = ms->visc[k] + ms->visc[k+1];
      for(d=0;d < 6;d++){
	a[2][d] -= 2.0 * vsum * afac * a[0][d] + drho * a[4][d];
	a[3][d] +=       vsum * afac * a[0][d];
	a[0][d] +=              afac * a[0][d];
      }
    }
    /*
       s_{k+1} - A s_k = f
    */
    for(c=0;c < 6;c++){
      row = ncmb + 6*k + c;
      col = 6*(k+1) + c;
      ms->ab[HC_BAND_INDEX(row,col,kl,ku,ldab)] = 1.0;
      for(d=0;d < 6;d++){
	col = 6*k + d;
	ms->ab[HC_BAND_INDEX(row,col,kl,ku,ldab)] = -a[c][d];
      }
    }
  }
  /*
     surface conditions, level n-1

     row A: u1 = 0
     row B: u4 = 0 (free slip) or u2 = plate velocity (no slip/plates)
     row C: (l+1) p1 + p2 - beta r_s (u3 - rho_s p1) = 0, i.e. the
            potential is regular above the surface after the gravity
            jump due to surface topography (not for l = 1, see above)
  */
  row = ncmb + 6*(n-1);
  col = 6*(n-1);
  ms->ab[HC_BAND_INDEX(row,col+0,kl,ku,ldab)] = 1.0;
  if(free_slip){
    ms->ab[HC_BAND_INDEX(row+1,col+3,kl,ku,ldab)] = 1.0;
  }else{
    ms->ab[HC_BAND_INDEX(row+1,col+1,kl,ku,ldab)] = 1.0;
  }
  if(l != 1){
    ms->ab[HC_BAND_INDEX(row+2,col+4,kl,ku,ldab)] = (el + 1.0) + beta * ms->r[n-1] * ms->rho[n-1];
    ms->ab[HC_BAND_INDEX(row+2,col+5,kl,ku,ldab)] = 1.0;
    ms->ab[HC_BAND_INDEX(row+2,col+2,kl,ku,ldab)] = -beta * ms->r[n-1];
  }
  ms->ncmb = ncmb;
  /*
     factor
  */
  ierr = hc_band_lu_factor(nsys,kl,ku,ms->ab,ms->ipiv);
  if(ierr){
    fprintf(stderr,"hc_ms_factor: singular system for l %i, column %i\n",l,ierr-1);
    exit(-1);
  }
  ms->lfact = l;
}
/*

assemble the right hand side for the density coefficients b[] (one
per density layer, as in hc_polsol) and the plate velocity
coefficient clm, solve, and return the state vectors at the output
layers in sol[nout*6] (layer-major), with the surface potential
derivative term replaced by the value above the surface, as in
hc_polsol

*/
void hc_ms_solve(struct hcs *hc, int l, HC_PREC *b, HC_PREC clm,
		 hc_boolean free_slip, int npb, HC_PREC *rpb, HC_PREC *fpb,
		 HC_PREC *sol)
{
  struct hc_ms *ms;
  int n,k,c,row,os,nsys;
  HC_PREC beta,dadd,*x;

  ms = &(hc->ms);
  if(ms->lfact != l){
    HC_ERROR("hc_ms_solve","system not factored for this l");
  }
  n = ms->n;
  nsys = ms->nsys;
  beta = hc->psp.beta;
  x = ms->rhs;
  for(k=0;k < nsys;k++){
    x[k] = 0.0;
  }
  /*
     sources at levels k+1, rows 3+6k+2 (u3) and 3+6k+5 (p2)
  */
  for(k=0;k < n-1;k++){
    row = ms->ncmb + 6*k;
    if(ms->isrc[k+1] >= 0){
      dadd = ms->den[k+1] * b[ms->isrc[k+1]];
      if(dadd != 0.0){
	x[row+2] += dadd;
	x[row+5] += beta * dadd * hc->rden[ms->isrc[k+1]];
      }
      if(ms->ipb[k+1] >= 0){
	x[row+2] -= fpb[ms->ipb[k+1]] * b[ms->isrc[k+1]];
	x[row+5] -= beta * rpb[ms->ipb[k+1]] * fpb[ms->ipb[k+1]] *
	  b[ms->isrc[k+1]] * hc->psp.rho_scale;
      }
    }
  }
  /*
     surface velocity
  */
  if(!free_slip){
    row = ms->ncmb + 6*(n-1) + 1;
    x[row] = clm;
  }
  hc_band_lu_solve(nsys,HC_MS_KL,HC_MS_KU,ms->ab,ms->ipiv,x);
  /*
     extract output layers
  */
  for(k=0;k < n;k++){
    if(ms->iout[k] >= 0){
      os = ms->iout[k] * 6;
      for(c=0;c < 6;c++){
	sol[os+c] = x[6*k+c];
      }
    }
  }
  /*
     potential derivative term above the surface
  */
  k = n-1;
  os = ms->iout[k] * 6;
  sol[os+5] = x[6*k+5] - beta * ms->r[k] * (x[6*k+2] - ms->rho[k] * x[6*k+4]);
}
