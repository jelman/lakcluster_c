function hmap_autumn_1(M_n_,A_n_rind_,A_n_cind,Z_n_rind_,Y_n_cind,rev_flag,pgap,ggap);
% hmap for AZWY;
if nargin<7; pgap = 16; end;
if nargin<8; ggap = 128; end;

nbins=length(M_n_);

M_A___n = []; M_Z___n = []; 
M_Y___n = []; M_W___n = []; 
for nb1=0:nbins-1;
if (rev_flag==0); 
M_A___n = [M_A___n ; M_n_{1+nb1}(A_n_rind_{1+nb1},A_n_cind) ]; 
M_Z___n = [M_Z___n ; M_n_{1+nb1}(Z_n_rind_{1+nb1},A_n_cind) ]; 
M_Y___n = [M_Y___n ; M_n_{1+nb1}(A_n_rind_{1+nb1},Y_n_cind) ]; 
M_W___n = [M_W___n ; M_n_{1+nb1}(Z_n_rind_{1+nb1},Y_n_cind) ]; 
end;%if (rev_flag==0); 
if (rev_flag==1); 
M_A___n = [M_A___n ; M_n_{1+nb1}(Z_n_rind_{1+nb1},A_n_cind) ]; 
M_Z___n = [M_Z___n ; M_n_{1+nb1}(A_n_rind_{1+nb1},A_n_cind) ]; 
M_Y___n = [M_Y___n ; M_n_{1+nb1}(Z_n_rind_{1+nb1},Y_n_cind) ]; 
M_W___n = [M_W___n ; M_n_{1+nb1}(A_n_rind_{1+nb1},Y_n_cind) ]; 
end;%if (rev_flag==1); 
end;%for nb1=0:nbins-1;

cmap_autumn = colormap('autumn'); 
cmap_autumn(1,:) = [1,1,1]; 
cmap_autumn(end-1,:) = 0.85*[1,1,1]; 
cmap_autumn(end-0,:) = 0.00*[1,1,1]; 
clen = size(cmap_autumn,1); 
colormap(cmap_autumn);
cbot=-1.5; ctop=1.5; cvals = linspace(cbot,ctop,clen); 
cbotl=cvals(2); ctopl = cvals(end-2); cdiff = cvals(end)-cvals(end-1); 
%subplot(1,1,1);
A13 = min(ctopl,max(cbotl,1*M_A___n));
A23 = cbot*ones(pgap,size(M_A___n,2));
A33 = min(ctopl,max(cbotl,1*M_Z___n));
A14 = cbot*ones(size(M_A___n,1),ggap);
A24 = cbot*ones(pgap,ggap);
A34 = cbot*ones(size(M_Z___n,1),ggap);
A15 = min(ctopl,max(cbotl,1*M_Y___n));
A25 = cbot*ones(pgap,size(M_Y___n,2));
A35 = min(ctopl,max(cbotl,1*M_W___n));
if 0;
elseif (size(M_Z___n,1)>0 & size(M_Y___n,2)>0);
imagesc( [ A13 , A14 , A15 ; A23 , A24 , A25 ; A33 , A34 , A35 ] , [cbot,ctop] );
elseif (size(M_Z___n,1)>0 & size(M_Y___n,2)==0);
imagesc( [ A13 ; A23 ; A33 ] , [cbot,ctop] );
elseif (size(M_Z___n,1)==0 & size(M_Y___n,2)>0);
imagesc( [ A13 , A14 , A15 ] , [cbot,ctop] );
elseif (size(M_Z___n,1)==0 & size(M_Y___n,2)==0);
imagesc( [ A13 ] , [cbot,ctop] );
end;%if
set(gca,'Xtick',[],'Ytick',[]);axis off; 
set(gcf,'Position',[100,100,1124,868/1.5]);