#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_An_Xn_ww(void *vp)
{
  int verbose=1;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  /* struct M_handle *M_St = (struct M_handle *)(vpra[ip++]); */ struct M_handle *M_St = NULL;
  struct M_handle *M_Xn = (struct M_handle *)(vpra[ip++]); 
  struct M_handle *M_Xt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_An_Xn_ww = (struct L_handle *)(vpra[ip++]);
  int trm_flag = *(int *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_s = *(int *)(vpra[ip++]);
  int output_spacing_w = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols);
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  double a0dh[A_pcols],a1dh[A_pcols];
  int ns_j=0,ns_a=0,ns_b=0,tab_s=0,tab_s_stride=0,ma_j=0,ma_a=0,ma_b=0,tab_a=0,tab_a_stride=0,nz_j=0,nz_a=0,nz_b=0,tab_z=0,nw_j=0,nw_a=0,nw_b=0,tab_w=0,tab_w_stride=0;
  double output_tmp=0;
  int vA=0;
  unsigned char *An_tag=NULL;
  __m128i *wAn_tag=NULL,*wXn_tag=NULL,*mcan_tag=NULL,*mcan_end=NULL;
  long long int n2=0;
  int M_St__rpop_j = /* (M_St!=NULL ? M_St->rpop_j : 1) */ 1;
  int M_St__rpop_b = /* (M_St!=NULL ? M_St->rpop_b : 1) */ 1;
  int M_St__nrows = /* (M_St!=NULL ? M_St->nrows : 1) */ 1;
  double output_An_Xn_base=0,output_an_Xn_base=0,output____Xn_base=0,output____Xn_base_[M_St__rpop_j],*dinp=NULL,dtmp=0,output_An_Xn_tmp=0,output_an_Xn_tmp=0,output____Xn_tmp_[M_Xt->rpop_j*M_St__rpop_j];
  int tab_x=0,mr=0,mx=0,mx_j=0;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering get_An_Xn_ww] tidx %d\n",tidx);}
  if (verbose>1){ printf(" %% Calculating output_An_Xn_ww\n");}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St__rpop_j; break; case SPACING_b: tab_s_stride = M_St__rpop_b; break; case SPACING_a: tab_s_stride = M_St__nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_An->rpop_j; break; case SPACING_b: tab_a_stride = M_An->rpop_b; break; case SPACING_a: tab_a_stride = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_w){ case SPACING_j: tab_w_stride = M_Xt->rpop_j; break; case SPACING_b: tab_w_stride = M_Xt->rpop_b; break; case SPACING_a: tab_w_stride = M_Xt->nrows; break; default: break; /* switch (output_spacing_w){ } */}
  output_An_Xn_ww->spacing_row = output_spacing_a; output_An_Xn_ww->row_stride = tab_a_stride;
  output_An_Xn_ww->spacing_col = output_spacing_w; output_An_Xn_ww->col_stride = tab_w_stride;
  output_An_Xn_ww->spacing_lyr = output_spacing_s; output_An_Xn_ww->lyr_stride = tab_s_stride;
  if (strstr(GLOBAL_skip,"An_Xn_ww")){ goto skip_An_Xn_ww;}
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    for (nz_a=0;nz_a<A_pcols;nz_a++){ a0dh[nz_a] = sqrt(D_An[nz_a]); a1dh[nz_a] = a_An[nz_a]*sqrt(D_An[nz_a]);}
    ns_j=0;
    while (ns_j<M_St__rpop_j){
      ns_a = ns_j; ns_b = ns_j; if (M_St!=NULL){ ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];} 
      switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
      wAn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
      wXn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
      mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
      dinp = &(a1dh[0]);
      dtmp = popcount_pm0_lf(&wAn_tag,&wXn_tag,&mcan_tag,&mcan_end,&dinp);
      output____Xn_base = dtmp*M_Xn->min_d*M_Xn->mlt_d;
      memset(output____Xn_tmp_,0,M_Xt->rpop_j*sizeof(double));
      nw_j=0;
      while (nw_j<M_Xt->rpop_j){
	nw_a = M_Xt->m_a_[nw_j]; nw_b = M_Xt->m_b_[nw_j];
	switch (output_spacing_w){ case SPACING_j: tab_w=nw_j; break; case SPACING_b: tab_w=nw_b; break; case SPACING_a: tab_w=nw_a; break; default: break; /* switch (output_spacing_w){ } */}
	output_an_Xn_tmp = 0; 
	mr=M_Xn->ncols_per_z*nw_j/* spacing_j */; n2=1;
	for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){
	  wAn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
	  wXn_tag = (__m128i*)&(M_Xn->wX[(mr+mx)*M_Xn->mc_length]);
	  mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	  dinp = &(a1dh[0]);
	  dtmp = popcount_pm0_lf(&wAn_tag,&wXn_tag,&mcan_tag,&mcan_end,&dinp);
	  output_an_Xn_tmp += n2*dtmp;
	  n2*=2;
	  /* for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){ } */}
	output____Xn_tmp_[nw_j] = output_an_Xn_tmp;
	nw_j++; /* while (nw_j<M_Xt->rpop_j){ } */}
      ma_j=0;
      while (ma_j<M_An->rpop_j){
	ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=ma_j; break; case SPACING_b: tab_a=ma_b; break; case SPACING_a: tab_a=ma_a; break; default: break; /* switch (output_spacing_a){ } */}
	wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	wXn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
	mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	dinp = &(a0dh[0]);
	dtmp = popcount_pm0_lf(&wAn_tag,&wXn_tag,&mcan_tag,&mcan_end,&dinp);
	output_An_Xn_base = dtmp*M_Xn->min_d*M_Xn->mlt_d;
	nw_j=0;
	while (nw_j<M_Xt->rpop_j){
	  nw_a = M_Xt->m_a_[nw_j]; nw_b = M_Xt->m_b_[nw_j];
	  switch (output_spacing_w){ case SPACING_j: tab_w=nw_j; break; case SPACING_b: tab_w=nw_b; break; case SPACING_a: tab_w=nw_a; break; default: break; /* switch (output_spacing_w){ } */}
	  tab_x = tab_a + tab_w*tab_a_stride + tab_s*tab_a_stride*tab_w_stride;
	  output_An_Xn_tmp = 0;
	  mr=M_Xn->ncols_per_z*nw_j/* spacing_j */; n2=1;
	  for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){
	    wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	    wXn_tag = (__m128i*)&(M_Xn->wX[(mr+mx)*M_Xn->mc_length]);
	    mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	    dinp = &(a0dh[0]);
	    dtmp = popcount_pm0_lf(&wAn_tag,&wXn_tag,&mcan_tag,&mcan_end,&dinp);
	    output_An_Xn_tmp += n2*dtmp;
	    n2*=2;
	    /* for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){ } */}
	  output_an_Xn_base = output____Xn_base; output_an_Xn_tmp = output____Xn_tmp_[nw_j];
	  output_An_Xn_ww->lf[tab_x] = (double)(output_An_Xn_base + output_An_Xn_tmp)/M_Xn->mlt_d - (double)(output_an_Xn_base + output_an_Xn_tmp)/M_Xn->mlt_d;
	  nw_j++; /* while (nw_j<M_Xt->rpop_j){ } */}
	ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
      GLOBAL_ops_count_one(tidx,0,1*(unsigned long long int)M_An->rpop_j*(unsigned long long int)M_Xt->rpop_j*(unsigned long long int)M_Xn->ncols_per_z*(unsigned long long int)M_Xn->mc_length*BIT8);
      ns_j++; /* while (ns_j<M_St__rpop_j){ } */}
    if (verbose>1){ raprintf(output_An_Xn_ww->lf,"double",tab_a_stride,tab_s_stride*tab_w_stride," %% output_An_Xn_ww->lf: ");}
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  if (GLOBAL_omp_type==GLOBAL_omp__on){
    for (nz_a=0;nz_a<A_pcols;nz_a++){ a0dh[nz_a] = sqrt(D_An[nz_a]); a1dh[nz_a] = a_An[nz_a]*sqrt(D_An[nz_a]);}
    memset(output____Xn_base_,0,M_St__rpop_j*sizeof(double));
    memset(output____Xn_tmp_,0,M_Xt->rpop_j*M_St__rpop_j*sizeof(double));
    ns_j=0;
    while (ns_j<M_St__rpop_j){
      ns_a = ns_j; ns_b = ns_j; if (M_St!=NULL){ ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];} 
      switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
      wAn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
      wXn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
      mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
      dinp = &(a1dh[0]);
      dtmp = popcount_pm0_lf(&wAn_tag,&wXn_tag,&mcan_tag,&mcan_end,&dinp);
      output____Xn_base_[ns_j] = dtmp*M_Xn->min_d*M_Xn->mlt_d;
      nw_j=0;
      while (nw_j<M_Xt->rpop_j){
	nw_a = M_Xt->m_a_[nw_j]; nw_b = M_Xt->m_b_[nw_j];
	switch (output_spacing_w){ case SPACING_j: tab_w=nw_j; break; case SPACING_b: tab_w=nw_b; break; case SPACING_a: tab_w=nw_a; break; default: break; /* switch (output_spacing_w){ } */}
	output_an_Xn_tmp = 0; 
	mr=M_Xn->ncols_per_z*nw_j/* spacing_j */; n2=1;
	for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){
	  wAn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
	  wXn_tag = (__m128i*)&(M_Xn->wX[(mr+mx)*M_Xn->mc_length]);
	  mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	  dinp = &(a1dh[0]);
	  dtmp = popcount_pm0_lf(&wAn_tag,&wXn_tag,&mcan_tag,&mcan_end,&dinp);
	  output_an_Xn_tmp += n2*dtmp;
	  n2*=2;
	  /* for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){ } */}
	output____Xn_tmp_[nw_j+ns_j*M_Xt->rpop_j] = output_an_Xn_tmp;
	nw_j++; /* while (nw_j<M_Xt->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_St__rpop_j){ } */}
    GLOBAL_ops_count_one(tidx,0,1*(unsigned long long int)M_St__rpop_j*(unsigned long long int)M_Xt->rpop_j*(unsigned long long int)M_Xn->ncols_per_z*(unsigned long long int)M_Xn->mc_length*BIT8);
#pragma omp parallel private(mx_j,ns_j,ns_a,ns_b,ma_j,ma_b,ma_a,nw_j,nw_b,nw_a,mr,mx,tab_s,tab_a,tab_x,tab_w,wAn_tag,wXn_tag,mcan_tag,mcan_end,n2,dtmp,dinp,output_An_Xn_base,output_an_Xn_base,output_An_Xn_tmp,output_an_Xn_tmp)
    { /* begin omp parallel */
      mx_j=0;
#pragma omp for schedule(dynamic)
      for (mx_j=0;mx_j<M_An->rpop_j*M_St__rpop_j;mx_j++){
	ns_j = mx_j / M_An->rpop_j; 
	ma_j = mx_j % M_An->rpop_j; 
	ns_a = ns_j; ns_b = ns_j; if (M_St!=NULL){ ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];} 
	switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
	ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=ma_j; break; case SPACING_b: tab_a=ma_b; break; case SPACING_a: tab_a=ma_a; break; default: break; /* switch (output_spacing_a){ } */}
	wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	wXn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
	mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	dinp = &(a0dh[0]);
	dtmp = popcount_pm0_lf(&wAn_tag,&wXn_tag,&mcan_tag,&mcan_end,&dinp);
	output_An_Xn_base = dtmp*M_Xn->min_d*M_Xn->mlt_d;
	nw_j=0;
	while (nw_j<M_Xt->rpop_j){
	  nw_a = M_Xt->m_a_[nw_j]; nw_b = M_Xt->m_b_[nw_j];
	  switch (output_spacing_w){ case SPACING_j: tab_w=nw_j; break; case SPACING_b: tab_w=nw_b; break; case SPACING_a: tab_w=nw_a; break; default: break; /* switch (output_spacing_w){ } */}
	  tab_x = tab_a + tab_w*tab_a_stride + tab_s*tab_a_stride*tab_w_stride;
	  output_An_Xn_tmp = 0;
	  mr=M_Xn->ncols_per_z*nw_j/* spacing_j */; n2=1;
	  for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){
	    wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	    wXn_tag = (__m128i*)&(M_Xn->wX[(mr+mx)*M_Xn->mc_length]);
	    mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	    dinp = &(a0dh[0]);
	    dtmp = popcount_pm0_lf(&wAn_tag,&wXn_tag,&mcan_tag,&mcan_end,&dinp);
	    output_An_Xn_tmp += n2*dtmp;
	    n2*=2;
	    /* for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){ } */}
	  output_an_Xn_base = output____Xn_base_[ns_j]; output_an_Xn_tmp = output____Xn_tmp_[nw_j+ns_j*M_Xt->rpop_j];
	  output_An_Xn_ww->lf[tab_x] = (double)(output_An_Xn_base + output_An_Xn_tmp)/M_Xn->mlt_d - (double)(output_an_Xn_base + output_an_Xn_tmp)/M_Xn->mlt_d;
	  nw_j++; /* while (nw_j<M_Xt->rpop_j){ } */}
	/* for (mx_j=0;mx_j<M_An->rpop_j*M_St__rpop_j;mx_j++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,0,(unsigned long long int)M_St__rpop_j*(unsigned long long int)M_An->rpop_j*(unsigned long long int)M_Xt->rpop_j*(unsigned long long int)M_Xn->ncols_per_z*(unsigned long long int)M_Xn->mc_length*BIT8);
    if (verbose>1){ raprintf(output_An_Xn_ww->lf,"double",tab_a_stride,tab_s_stride*tab_w_stride," %% output_An_Xn_ww->lf: ");}
    /* if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
 skip_An_Xn_ww:
  if (verbose>1){ printf(" %% [finished get_An_Xn_ww] tidx %d\n",tidx);}
  return NULL;
}

int wrap_An_Xn_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,int spacing_w,struct M_handle *M_An,/* struct M_handle *M_St, */struct M_handle *M_Xn,struct M_handle *M_Xt,double *A_ajdk,int trm_flag,struct L_handle **output_An_Xn_ww_p)
{
  /* calls get_An_Xn_ww;
     variable space in **vpra (should be at least size 10)
   */
  int verbose=0;
  unsigned long long int length_a=0,length_s=0,length_w=0,length=0,ip=0;
  int M_St__rpop_j = /* (M_St!=NULL ? M_St->rpop_j : 1) */ 1;
  int M_St__rpop_b = /* (M_St!=NULL ? M_St->rpop_b : 1) */ 1;
  int M_St__nrows = /* (M_St!=NULL ? M_St->nrows : 1) */ 1;
  if (verbose){ printf(" %% [entering wrap_An_Xn_ww__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  /* if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");} */
  if (verbose){ M_handle_printf(M_Xt,verbose," %% M_Xt: ");}
  switch (spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: length_s = M_St__rpop_j; break; case SPACING_b: length_s = M_St__rpop_b; break; case SPACING_a: length_s = M_St__nrows; break; default: break; /* switch (spacing_s){ } */}
  switch (spacing_w){ case SPACING_j: length_w = M_Xt->rpop_j; break; case SPACING_b: length_w = M_Xt->rpop_b; break; case SPACING_a: length_w = M_Xt->nrows; break; default: break; /* switch (spacing_w){ } */}
  length = length_a*length_s*length_w; if (verbose){ printf(" %% length %llu*%llu*%llu=%llu\n",length_a,length_s,length_w,length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Xt->mr_b,M_Xt->bitj,1,M_Xt->nrows," %% M_Xt->mr_b: ");}
  if (verbose>2){ bprintf(M_Xt->mr_j,M_Xt->bitj,1,M_Xt->nrows," %% M_Xt->mr_j: ");}
  if (verbose>2){ bprintf(M_Xt->mc_b,M_Xt->bitj,1,M_Xt->ncols," %% M_Xt->mc_b: ");}
  if (verbose>2){ bprintf(M_Xt->mc_j,M_Xt->bitj,1,M_Xt->ncols," %% M_Xt->mc_j: ");}
  if (*output_An_Xn_ww_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_Xn_ww_p = L_handle_make(length);}
  if ((*output_An_Xn_ww_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_An_Xn_ww__run\n",(*output_An_Xn_ww_p)->length,length);}
  memset((*output_An_Xn_ww_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; /* vpra[ip++] = M_St; */ vpra[ip++] = M_Xn; vpra[ip++] = M_Xt; vpra[ip++] = A_ajdk; vpra[ip++] = *output_An_Xn_ww_p;
  switch (trm_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break; default: break; /* switch (trm_flag){ } */}
  switch (spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_s){ } */}
  switch (spacing_w){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_w){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_Xn_ww,vpra)){ printf("Warning! cannot create thread %d in wrap_An_Xn_ww__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_Xn_ww(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_An_Xn_ww__run] tidx %d\n",*tidx);}
  return length;
}

void *get_An_Xn_uu(void *vp)
{
  int verbose=1;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  /* struct M_handle *M_St = (struct M_handle *)(vpra[ip++]); */ struct M_handle *M_St = NULL;
  struct L_handle *lf_Xn = (struct L_handle *)(vpra[ip++]);
  struct M_handle *M_Xt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_An_Xn_uu = (struct L_handle *)(vpra[ip++]);
  int trm_flag = *(int *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_s = *(int *)(vpra[ip++]);
  int output_spacing_w = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols);
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int ns_j=0,ns_a=0,ns_b=0,tab_s=0,tab_s_stride=0,ma_j=0,ma_a=0,ma_b=0,tab_a=0,tab_a_stride=0,nz_j=0,nz_a=0,nz_b=0,tab_z=0,nw_j=0,nw_a=0,nw_b=0,tab_w=0,tab_w_stride=0;
  double output_tmp=0;
  int vA=0;
  unsigned char *An_tag=NULL;
  __m128i *wAn_tag=NULL,*wXn_tag=NULL,*mcan_tag=NULL,*mcan_end=NULL;
  long long int n2=0;
  int M_St__rpop_j = /* (M_St!=NULL ? M_St->rpop_j : 1) */ 1;
  int M_St__rpop_b = /* (M_St!=NULL ? M_St->rpop_b : 1) */ 1;
  int M_St__nrows = /* (M_St!=NULL ? M_St->nrows : 1) */ 1;
  double output_An_Xn_base=0,output_an_Xn_base=0,output____Xn_base=0,*dinp=NULL,dtmp=0,output_An_Xn_tmp=0,output_an_Xn_tmp=0,output____Xn_tmp_[M_Xt->nrows];
  int tab_x=0,mr=0,mx=0;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering get_An_Xn_uu] tidx %d\n",tidx);}
  if (verbose>1){ printf(" %% Calculating output_An_Xn_uu\n");}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St__rpop_j; break; case SPACING_b: tab_s_stride = M_St__rpop_b; break; case SPACING_a: tab_s_stride = M_St__nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_An->rpop_j; break; case SPACING_b: tab_a_stride = M_An->rpop_b; break; case SPACING_a: tab_a_stride = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_w){ case SPACING_j: tab_w_stride = M_Xt->rpop_j; break; case SPACING_b: tab_w_stride = M_Xt->rpop_b; break; case SPACING_a: tab_w_stride = M_Xt->nrows; break; default: break; /* switch (output_spacing_w){ } */}
  if (verbose>1){ printf(" %% tab_s_stride %d tab_a_stride %d tab_w_stride %d\n",tab_s_stride,tab_a_stride,tab_w_stride);}
  output_An_Xn_uu->spacing_row = output_spacing_a; output_An_Xn_uu->row_stride = tab_a_stride;
  output_An_Xn_uu->spacing_col = output_spacing_w; output_An_Xn_uu->col_stride = tab_w_stride;
  output_An_Xn_uu->spacing_lyr = output_spacing_s; output_An_Xn_uu->lyr_stride = tab_s_stride;
  if (strstr(GLOBAL_skip,"An_Xn_uu")){ goto skip_An_Xn_uu;}
  if (verbose>2){ M_handle_printf(M_An,1," %% M_An: ");}
  if (verbose>2){ M_handle_printf(M_St,1," %% M_St: ");}
  if (verbose>1){ raprintf(lf_Xn->lf,"double",lf_Xn->row_stride,lf_Xn->col_stride*lf_Xn->lyr_stride," %% lf_Xn->lf: ");}
  if (verbose>1){ raprintf(D_An,"double",1,A_pcols," %% D_An: "); raprintf(a_An,"double",1,A_pcols," %% a_An: ");}
  if (verbose>1){ printf(" %% POPLENGTH %d M_An->mc_length %d\n",POPLENGTH,M_An->mc_length);}
  if (verbose>1){ if (M_St!=NULL){ printf(" %% POPLENGTH %d M_St->mc_length %d\n",POPLENGTH,M_St->mc_length);}}
  ns_j=0;
  while (ns_j<M_St__rpop_j){
    ns_a = ns_j; ns_b = ns_j; if (M_St!=NULL){ ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];}
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
    ma_j=0;
    while (ma_j<M_An->rpop_j){
      ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=ma_j; break; case SPACING_b: tab_a=ma_b; break; case SPACING_a: tab_a=ma_a; break; default: break; /* switch (output_spacing_a){ } */}
      An_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length]));
      nw_j=0;
      while (nw_j<M_Xt->rpop_j){
	nw_a = M_Xt->m_a_[nw_j]; nw_b = M_Xt->m_b_[nw_j];
	switch (output_spacing_w){ case SPACING_j: tab_w=nw_j; break; case SPACING_b: tab_w=nw_b; break; case SPACING_a: tab_w=nw_a; break; default: break; /* switch (output_spacing_w){ } */}
	output_tmp=0;
	nz_j=0;
	while (nz_j<M_An->cpop_j){
	  nz_a = M_An->n_a_[nz_j]; nz_b = M_An->n_b_[nz_j];
	  vA = bget____(An_tag,nz_a);
	  output_tmp += (vA - a_An[nz_a/POPLENGTH])*sqrt(D_An[nz_a/POPLENGTH]) * (trm_flag ? (*L3_get(lf_Xn,nw_j,nw_b,nw_a,nz_j,nz_b,nz_a,ns_j,ns_b,ns_a)) : (*L3_get(lf_Xn,nz_j,nz_b,nz_a,nw_j,nw_b,nw_a,ns_j,ns_b,ns_a)));
	  nz_j++; /* while (nz_j<M_An->cpop_j){ } */}
	if (verbose>3){ printf("(%d,%d,%d) %lf ",tab_a,tab_w,tab_s,output_tmp);}
	output_An_Xn_uu->lf[tab_a + tab_w*tab_a_stride + tab_s*tab_a_stride*tab_w_stride] = output_tmp;
	nw_j++; /* while (nw_j<M_Xt->rpop_j){ } */}
      if (verbose>3){ printf("\n");}
      ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St__rpop_j){ } */}
  if (verbose>1){ raprintf(output_An_Xn_uu->lf,"double",tab_a_stride,tab_w_stride*tab_s_stride," %% output_An_Xn_uu->lf: ");}
  GLOBAL_ops_count_one(tidx,(unsigned long long int)M_St__rpop_j*(unsigned long long int)M_An->rpop_j*(unsigned long long int)M_Xt->rpop_j,0);
 skip_An_Xn_uu:
  if (verbose>1){ printf(" %% [finished get_An_Xn_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_An_Xn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,int spacing_w,struct M_handle *M_An,/* struct M_handle *M_St, */struct L_handle *lf_Xn,struct M_handle *M_Xt,double *A_ajdk,int trm_flag,struct L_handle **output_An_Xn_uu_p)
{
  /* calls get_An_Xn_uu;
     variable space in **vpra (should be at least size 10)
   */
  int verbose=0;
  unsigned long long int length_a=0,length_s=0,length_w=0,length=0,ip=0;
  int M_St__rpop_j = /* (M_St!=NULL ? M_St->rpop_j : 1) */ 1;
  int M_St__rpop_b = /* (M_St!=NULL ? M_St->rpop_b : 1) */ 1;
  int M_St__nrows = /* (M_St!=NULL ? M_St->nrows : 1) */ 1;
  if (verbose){ printf(" %% [entering wrap_An_Xn_uu__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  /* if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");} */
  if (verbose){ M_handle_printf(M_Xt,verbose," %% M_Xt: ");}
  switch (spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: length_s = M_St__rpop_j; break; case SPACING_b: length_s = M_St__rpop_b; break; case SPACING_a: length_s = M_St__nrows; break; default: break; /* switch (spacing_s){ } */}
  switch (spacing_w){ case SPACING_j: length_w = M_Xt->rpop_j; break; case SPACING_b: length_w = M_Xt->rpop_b; break; case SPACING_a: length_w = M_Xt->nrows; break; default: break; /* switch (spacing_w){ } */}
  length = length_a*length_s*length_w; if (verbose){ printf(" %% length %llu*%llu*%llu=%llu\n",length_a,length_s,length_w,length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Xt->mr_b,M_Xt->bitj,1,M_Xt->nrows," %% M_Xt->mr_b: ");}
  if (verbose>2){ bprintf(M_Xt->mr_j,M_Xt->bitj,1,M_Xt->nrows," %% M_Xt->mr_j: ");}
  if (verbose>2){ bprintf(M_Xt->mc_b,M_Xt->bitj,1,M_Xt->ncols," %% M_Xt->mc_b: ");}
  if (verbose>2){ bprintf(M_Xt->mc_j,M_Xt->bitj,1,M_Xt->ncols," %% M_Xt->mc_j: ");}
  if (*output_An_Xn_uu_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_Xn_uu_p = L_handle_make(length);}
  if ((*output_An_Xn_uu_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_An_Xn_uu__run\n",(*output_An_Xn_uu_p)->length,length);}
  memset((*output_An_Xn_uu_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; /* vpra[ip++] = M_St; */ vpra[ip++] = lf_Xn; vpra[ip++] = M_Xt; vpra[ip++] = A_ajdk; vpra[ip++] = *output_An_Xn_uu_p; 
  switch (trm_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break; default: break; /* switch (trm_flag){ } */}
  switch (spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_s){ } */}
  switch (spacing_w){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_w){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_Xn_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_An_Xn_uu__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_Xn_uu(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_An_Xn_uu__run] tidx %d\n",*tidx);}
  return length;
}

void wrap_An_Xn_ww_test()
{
  /* test for errors with input file: An_Xn_ww_error.in ;
  */
  int St_flag = 0;
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  int X_n_cols=0,nr_j=0,nr_b=0,nr_a=0,nc_j=0,nc_b=0,nc_a=0;
  struct L_handle **lf_Xn=NULL;
  unsigned char *bmc1_p=NULL;
  int ns_j=0,ns_b=0,ns_a=0;
  struct M_handle **M_Xn=NULL,**M_Xt=NULL;
  struct L_handle **lf_An_Xn_ww=NULL; int *length_An_Xn_ww=NULL;
  struct L_handle **lf_An_Xn_uu=NULL; int *length_An_Xn_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_An_Xn_ww_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  nb=0;
  int M_St_nb___rpop_j = (St_flag ? M_St[nb]->rpop_j : 1);
  int M_St_nb___rpop_b = (St_flag ? M_St[nb]->rpop_b : 1);
  int M_St_nb___nrows = (St_flag ? M_St[nb]->nrows : 1);
  X_n_cols = 8;
  bmc1_p = wkspace_all0c(bsize(X_n_cols)); 
  bset_off(bmc1_p,0); bset__on(bmc1_p,1); bset_off(bmc1_p,2); bset__on(bmc1_p,3); bset_off(bmc1_p,4); bset_off(bmc1_p,5); bset__on(bmc1_p,6); bset__on(bmc1_p,7); 
  M_Xt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_Xt[nb] = M_handle_v_make(BITJ,X_n_cols,A_n_cols,NULL,NULL,bmc1_p,M_An[nb]->mc_b);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_Xn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_Xn[nb] = L_handle_make((unsigned long long int)X_n_cols*(unsigned long long int)A_n_cols*(unsigned long long int)M_St_nb___nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  M_Xn = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_Xn[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,A_n_cols,X_n_cols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_An_Xn_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_Xn_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_An_Xn_ww[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_St_nb___nrows*(unsigned long long int)X_n_cols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_An_Xn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_Xn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_An_Xn_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_St_nb___nrows*(unsigned long long int)X_n_cols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
  	if (verbose){ printf(" %% %s; %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_An[nb],A_ajdk,&(lf_An_ajdk[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	for (nb=0;nb<nbins;nb++){
	  M_mxset(M_Xt[nb],bmc1_p,M_An[nb]->mc_j);
	  if (verbose>2){ M_handle_printf(M_Xt[nb],1," %% M_Xt[nb]: ");}
	  lf_Xn[nb]->spacing_row = n_spacing_A; 
	  switch(lf_Xn[nb]->spacing_row){ case SPACING_j: lf_Xn[nb]->row_stride = M_Xt[nb]->cpop_j; break; case SPACING_b: lf_Xn[nb]->row_stride = M_Xt[nb]->cpop_b; break; case SPACING_a: lf_Xn[nb]->row_stride = M_Xt[nb]->ncols; break; default: break; /* switch(lf_Xn[nb]->spacing_row){ } */}
	  lf_Xn[nb]->spacing_col = n_spacing_A;
	  switch(lf_Xn[nb]->spacing_col){ case SPACING_j: lf_Xn[nb]->col_stride = M_Xt[nb]->rpop_j; break; case SPACING_b: lf_Xn[nb]->col_stride = M_Xt[nb]->rpop_b; break; case SPACING_a: lf_Xn[nb]->col_stride = M_Xt[nb]->nrows; break; default: break; /* switch(lf_Xn[nb]->spacing_col){ } */}
	  lf_Xn[nb]->spacing_lyr = n_spacing_A;
	  switch(lf_Xn[nb]->spacing_lyr){ case SPACING_j: lf_Xn[nb]->lyr_stride = M_St_nb___rpop_j; break; case SPACING_b: lf_Xn[nb]->lyr_stride = M_St_nb___rpop_b; break; case SPACING_a: lf_Xn[nb]->lyr_stride = M_St_nb___nrows; break; default: break; /* switch(lf_Xn[nb]->spacing_lyr){ } */}
	  for (ns_j=0;ns_j<M_St_nb___rpop_j;ns_j++){
	    ns_a = ns_j; ns_b = ns_j; if (St_flag){ ns_a = M_St[nb]->m_a_[ns_j]; ns_b = M_St[nb]->m_b_[ns_j];}
	    for (nc_j=0;nc_j<M_Xt[nb]->rpop_j;nc_j++){ 
	      nc_a = M_Xt[nb]->m_a_[nc_j]; nc_b = M_Xt[nb]->m_b_[nc_j];
	      for (nr_j=0;nr_j<M_Xt[nb]->cpop_j;nr_j++){
		nr_a = M_Xt[nb]->n_a_[nr_j]; nr_b = M_Xt[nb]->n_b_[nr_j];
		if (verbose>2){ printf(" %% nb: %d --> nr (%d,%d,%d)/%d , nc (%d,%d,%d)/%d\n",nb,nr_j,nr_b,nr_a,M_Xt[nb]->cpop_j,nc_j,nc_b,nc_a,M_Xt[nb]->rpop_j);}
		L3_set(lf_Xn[nb],nr_j,nr_b,nr_a,nc_j,nc_b,nc_a,ns_j,ns_b,ns_a,rand01);
		/* for (nr_j=0;nr_j<M_Xt[nb]->cpop_j;nr_j++){ } */}
	      /* for (nc_j=0;nc_j<M_Xt[nb]->rpop_j;nc_j++){ } */}
	    /* for (ns_j=0;ns_j<M_St_nb___rpop_j;ns_j++){ } */}
	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
	  for (ns_j=0;ns_j<M_St_nb___rpop_j;ns_j++){
	    ns_a = ns_j; ns_b = ns_j; if (St_flag){ ns_a = M_St[nb]->m_a_[ns_j]; ns_b = M_St[nb]->m_b_[ns_j];}
	    GLOBAL_pthread_tic(); 
	    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Xt[nb]->mc_j,M_Xt[nb]->mc_b,M_Xt[nb]->mr_j,M_Xt[nb]->mr_b,lf_Xn[nb],L3_lf_get(lf_Xn[nb],lf_Xn[nb]->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_Xt[nb]->ncols),&(M_Xt[nb]->nrows),&(M_Xn[nb]),&(GLOBAL_B_MLT),(addressable_0));
	    GLOBAL_pthread_toc();
	    /* for (ns_j=0;ns_j<M_St_nb___rpop_j;ns_j++){ } */}
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% xcalc: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
	  length_An_Xn_ww[nb] = wrap_An_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_An[nb],M_Xn[nb],M_Xt[nb],A_ajdk,(addressable_0),&(lf_An_Xn_ww[nb]));
	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_Xn_ww : ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic(); 
	    length_An_Xn_uu[nb] = wrap_An_Xn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_An[nb],lf_Xn[nb],M_Xt[nb],A_ajdk,(addressable_0),&(lf_An_Xn_uu[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc(); 
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_Xn_uu : ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_An_Xn_ww[%d] error %0.16f\n",nb,dra_diff(lf_An_Xn_ww[nb]->lf,lf_An_Xn_uu[nb]->lf,length_An_Xn_ww[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
  	/* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_An_Xn_ww_test]\n");}
}
