function A_ajdk_ = calc_A_ajdk(An_pcols,A_p_);
AJDK_AXDX_define;
A_ajdk_ = zeros(AJDK_TOT*An_pcols,1);
%%%%%%%%%%%%%%%% ;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_0_0*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^0 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^0 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_1_0*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^1 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^0 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_2_0*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^2 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^0 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_3_0*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^3 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^0 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_4_0*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^4 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^0 ; end;%for nc=0:An_pcols-1;
%%%%%%%%%%%%%%%% ;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_0_1*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^0 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^1 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_1_1*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^1 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^1 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_2_1*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^2 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^1 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_3_1*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^3 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^1 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_4_1*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^4 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^1 ; end;%for nc=0:An_pcols-1;
%%%%%%%%%%%%%%%% ;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_0_2*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^0 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^2 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_1_2*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^1 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^2 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_2_2*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^2 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^2 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_3_2*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^3 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^2 ; end;%for nc=0:An_pcols-1;
for nc=0:An_pcols-1; A_ajdk_(1+nc+AJDK_4_2*An_pcols) = (A_p_(1+nc)-(1-A_p_(1+nc))).^4 * (1.0./max(0.01,(4.0*A_p_(1+nc)*(1-A_p_(1+nc))))).^2 ; end;%for nc=0:An_pcols-1;