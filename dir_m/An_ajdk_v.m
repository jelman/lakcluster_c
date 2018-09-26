function output_ = An_ajdk_v(uAn_,A_ajdk_,umc_,POPLENGTH);

if (nargin<1);

disp(sprintf(' '));
disp(' testing An_ajdk_v: ');
AJDK_AXDX_define;
popcount_ = make_popcount_();
BIT8 = 8; bitj = 16; POPLENGTH = 1920; POPSUB = POPLENGTH/BIT8;
At_nrows = round(2*POPLENGTH + 13);
At_ncols = 12; Zt_ncols = 17;
At_ = randn(At_nrows,At_ncols)>0;
Zt_ = randn(At_nrows,Zt_ncols)>0;
uAn_ = tutorial_integer_to_binary(bitj,At_nrows,At_ncols,At_);
uZn_ = tutorial_integer_to_binary(bitj,At_nrows,Zt_ncols,Zt_);
iAn_ = tutorial_binary_to_integer(bitj,At_nrows,At_ncols,uAn_);
iZn_ = tutorial_binary_to_integer(bitj,At_nrows,Zt_ncols,uZn_);
disp(sprintf('error |At_-iAn_| = %f',norm(At_-iAn_)));
disp(sprintf('error |Zt_-iZn_| = %f',norm(Zt_-iZn_)));
At_mr_ = 1.0*(randn(At_nrows,1)>0);
uAn_mr_ = tutorial_integer_to_binary(bitj,At_nrows,1,At_mr_);
At_prows = rup(At_nrows,POPLENGTH)/POPLENGTH;
A_p_ = 0.25 + 0.5*rand(At_prows,1);
A_ajdk_ = calc_A_ajdk(At_prows,A_p_);
An_ajdk_ = An_ajdk_v(uAn_,A_ajdk_,uAn_mr_,POPLENGTH);
A1_ = zeros(At_ncols*AJDK_TOT,1);
for nc1=0:At_ncols-1;
for nr1=0:At_nrows-1;
tmp_At = 2*iAn_(1+nr1,1+nc1)-1;
tmp_a = A_ajdk_(1 + floor(nr1/POPLENGTH) + AJDK_1_0*At_prows);
tmp_D = A_ajdk_(1 + floor(nr1/POPLENGTH) + AJDK_0_1*At_prows);
ajdk_x_x=0;
for nd=0:2; for na=0:4;
tmp_A1 = tmp_At * At_mr_(1+nr1) * tmp_a.^na * tmp_D.^nd;
A1_(1+nc1 + ajdk_x_x*At_ncols) = A1_(1+nc1 + ajdk_x_x*At_ncols) + tmp_A1;
ajdk_x_x = ajdk_x_x+1;
end;end;%for nd=0:2; for na=0:4;
end;%for nr1=0:At_nrows-1;
end;%for nc1=0:At_ncols-1;
A1_ = reshape(A1_,At_ncols,AJDK_TOT);
disp(sprintf('error |An_ajdk_-A1_| = %f',norm(An_ajdk_-A1_)));
Zn_ajdk_ = An_ajdk_v(uZn_,A_ajdk_,uAn_mr_,POPLENGTH);
Z1_ = zeros(Zt_ncols*AJDK_TOT,1);
for nc1=0:Zt_ncols-1;
for nr1=0:At_nrows-1;
tmp_Zt = 2*iZn_(1+nr1,1+nc1)-1;
tmp_a = A_ajdk_(1 + floor(nr1/POPLENGTH) + AJDK_1_0*At_prows);
tmp_D = A_ajdk_(1 + floor(nr1/POPLENGTH) + AJDK_0_1*At_prows);
ajdk_x_x=0;
for nd=0:2; for na=0:4;
tmp_Z1 = tmp_Zt * At_mr_(1+nr1) * tmp_a.^na * tmp_D.^nd;
Z1_(1+nc1 + ajdk_x_x*Zt_ncols) = Z1_(1+nc1 + ajdk_x_x*Zt_ncols) + tmp_Z1;
ajdk_x_x = ajdk_x_x+1;
end;end;%for nd=0:2; for na=0:4;
end;%for nr1=0:At_nrows-1;
end;%for nc1=0:Zt_ncols-1;
Z1_ = reshape(Z1_,Zt_ncols,AJDK_TOT);
disp(sprintf('error |Zn_ajdk_-Z1_| = %f',norm(Zn_ajdk_-Z1_)));

return;
end;%if (nargin<1);

if (nargin<4); POPLENGTH = 1920; end;
AJDK_AXDX_define;
BIT8=8; POPSUB = POPLENGTH/BIT8;
[An_bcols,An_nrows] = size(uAn_);
An_pcols = rup(An_bcols,POPSUB)/POPSUB;
popcount_ = make_popcount_();
for ajdk_x_x=0:AJDK_TOT-1;
dtmp_base = 0;
for nc=0:An_bcols-1;
dtmp_base = dtmp_base + A_ajdk_(1+floor(nc/POPSUB) + ajdk_x_x*An_pcols)*double(intlut(umc_(1+nc),popcount_));
end;%for nc=0:An_bcols-1;
for nr=0:An_nrows-1;
dtmp = dtmp_base;
output_(1+nr + ajdk_x_x*An_nrows) = -dtmp;
dtmp = popcount_lf(uAn_(:,1+nr),umc_,An_bcols,A_ajdk_(1 + (0:An_pcols-1) + ajdk_x_x*An_pcols),popcount_,POPLENGTH);
output_(1+nr + ajdk_x_x*An_nrows) = output_(1+nr + ajdk_x_x*An_nrows) + 2*dtmp;
end;%for nr=0:An_nrows-1;
end;%for ajdk_x_x=0:AJDK_TOT-1;

output_ = reshape(output_,An_nrows,AJDK_TOT);

