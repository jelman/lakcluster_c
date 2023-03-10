%%%%%%%%;
% This script is technically entitled Up99_vs_Ap05, ;
% because I wanted to compare the Up99 vs the Ap05 data. ;
% However, because the Up99 data is large (and takes a while to load), ;
% I used the smaller Up05 dataset instead. ;
% Thus, if you really want to use the Up99 data, ;
% replace all instances of Up05 below with Up99. ;
%%%%%%%%;

clear;
run('/home/jelman/Github/lakcluster_c/dir_m/setup_0'); %<-- set up the paths. ;

flag_verbose = 1;
flag_disp = 1+flag_verbose; nf=0;
flag_replot = 0;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Comparing Up05 with Ap05 data. ;')); end;
if (flag_verbose); disp(sprintf(' %% Assume path is set using dir_lakcluster_c/dir_m/setup_0.m. ;')); end;

if (flag_verbose); disp(sprintf(' %% Comparing Up99 with Up05 data. ;')); end;
memory_GB = 128; %<-- maybe this should be increased? ;
if (flag_verbose); disp(sprintf(' %% trying with memory_GB %d ;',memory_GB)); end;

dir_code = '/home/jelman/Github/lakcluster_c';
dir_trunk = '/home/jelman/Projects/AD_Biclustering/data/ADNI/ADNI_vs_UKB';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

n_dataset = 2;
str_dataset_ = {'Up05','Ap05'};
dataset_ = cell(n_dataset,1);
for ndataset=0:n_dataset-1;
dataset = struct('type','dataset');
dataset.str_dataset = str_dataset_{1+ndataset};
dataset.dir_trunk = sprintf('%s/dir_%s',dir_trunk,dataset.str_dataset);
dataset.str_prefix = 'test2mds_maf01';
dataset.dir_0in = sprintf('%s/dir_%s',dataset.dir_trunk,dataset.str_prefix);
%%%%%%%%;
dataset.fname_bimext = sprintf('%s/%s_bim.ext',dataset.dir_0in,dataset.str_prefix);
if ~exist(dataset.fname_bimext,'file');
dataset.fname_bimext = sprintf('%s/%s_bim.ext',dataset.dir_0in,dataset.str_prefix);
end;%if ~exist(dataset.fname_bimext,'file');
if ~exist(dataset.fname_bimext,'file');
dataset.fname_bimext = sprintf('%s/test0mds_maf01_bim.ext',dataset.dir_0in);
end;%if ~exist(dataset.fname_bimext,'file');
if ~exist(dataset.fname_bimext,'file');
dataset.fname_bimext = sprintf('%s/%s_bim.ext',dataset.dir_0in,dataset.str_dataset);
end;%if ~exist(dataset.fname_bimext,'file');
if (flag_verbose); disp(sprintf(' %% attempting to load %s',dataset.fname_bimext)); end;
tmp_t = tic();
[ ...
,dataset.n_snp ...
,dataset.bim_khr_ ...
,dataset.bim_vid_ ...
,dataset.bim_gdi_ ...
,dataset.bim_pdi_ ...
,dataset.bim_al1_ ...
,dataset.bim_al2_ ...
,dataset.bim_alt_ ...
,dataset.bim_ent_ ...
,dataset.bim_frq_ ...
,dataset.bim_mss_ ...
,dataset.bim_maf_ ...
,dataset.bim_name_ ...
,~ ...
] = ...
load_bimext_ver1( ...
dataset.fname_bimext ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% load_bimext_ver1: %0.6fs',tmp_t)); end;
%%%%%%%%;
dataset.fname_famext = sprintf('%s/%s_fam.ext',dataset.dir_0in,dataset.str_prefix);
if ~exist(dataset.fname_famext,'file');
dataset.fname_famext = sprintf('%s/%s_fam.ext',dataset.dir_0in,dataset.str_prefix);
end;%if ~exist(dataset.fname_famext,'file');
if ~exist(dataset.fname_famext,'file');
dataset.fname_famext = sprintf('%s/test0mds_maf01_fam.ext',dataset.dir_0in);
end;%if ~exist(dataset.fname_famext,'file');
if ~exist(dataset.fname_famext,'file');
dataset.fname_famext = sprintf('%s/%s_fam.ext',dataset.dir_0in,dataset.str_dataset);
end;%if ~exist(dataset.fname_famext,'file');
[ ...
,dataset.n_patient ...
,dataset.fam_fid_ ...
,dataset.fam_iid_ ...
,dataset.fam_yid_ ...
,dataset.fam_xid_ ...
,dataset.fam_sex_ ...
,dataset.fam_dvx_ ...
,dataset.fam_dir_ ...
,dataset.fam_fidandiid_ ...
,~ ...
  ] = ...
load_famext_ver1( ...
 dataset.fname_famext ...
);
%%%%%%%%;
dataset_{1+ndataset} = dataset;
clear dataset;
end;%for ndataset=0:n_dataset-1;

for nd=0:1;
if (flag_verbose); disp(sprintf(' %% nd %d <-- %s',nd,dataset_{1+nd}.dir_0in)); end;
[~,n_row,n_col] = binary_getsize(sprintf('%s/%s_A_full_n.b16',dataset_{1+nd}.dir_0in,dataset_{1+nd}.str_prefix));
if (flag_verbose); disp(sprintf(' %% %% A_full_n_.b16 [%d,%d]',n_row,n_col)); end;
if (flag_verbose); disp(sprintf(' %% %% n_patient [%d]',dataset_{1+nd}.n_patient)); end;
if (flag_verbose); disp(sprintf(' %% %% n_snp [%d]',dataset_{1+nd}.n_snp)); end;
end;%for nd=0:1;

for ndataset=0:n_dataset-1;
parameter = struct('type','parameter');
parameter.dir_0in = dataset_{1+ndataset}.dir_0in;
parameter.str_prefix = dataset_{1+ndataset}.str_prefix;
if (flag_verbose); disp(sprintf(' %% attempting to load mx__ from %s',parameter.dir_0in)); end;
tmp_t = tic();
mx__ = load_mx__from_parameter_ver0(parameter);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% load_mx__from_parameter_ver0: %0.6fs',tmp_t)); end;
dataset_{1+ndataset}.parameter = parameter;
dataset_{1+ndataset}.mx__ = mx__;
end;%for ndataset=0:n_dataset-1;

ndataset_Up05 = 0;
ndataset_Ap05 = 1;
n_snp_Up05 = dataset_{1+ndataset_Up05}.n_snp;
n_snp_Ap05 = dataset_{1+ndataset_Ap05}.n_snp;
if (flag_verbose); disp(sprintf(' %% intersecting allele subsets')); end;
tmp_t = tic();
[allele_cap_,ij_Up05_from_cap_,ij_Ap05_from_cap_] = intersect(dataset_{1+ndataset_Up05}.bim_name_,dataset_{1+ndataset_Ap05}.bim_name_,'stable');
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% intersect: %0.6fs',tmp_t)); end;
index_Up05_from_cap_ = ij_Up05_from_cap_ - 1;
index_Ap05_from_cap_ = ij_Ap05_from_cap_ - 1;
ij_Up05_from_Ap05__ = sparse(ij_Up05_from_cap_,ij_Ap05_from_cap_,1,n_snp_Up05,n_snp_Ap05);
ij_Ap05_from_Up05__ = sparse(ij_Ap05_from_cap_,ij_Up05_from_cap_,1,n_snp_Ap05,n_snp_Up05);
index_Up05_setminus_cap_ = setdiff(efind(dataset_{1+ndataset_Up05}.mx__.mc_A_),index_Up05_from_cap_);
index_Ap05_setminus_cap_ = setdiff(efind(dataset_{1+ndataset_Ap05}.mx__.mc_A_),index_Ap05_from_cap_);
if (flag_verbose);
disp(sprintf(' %% number of alleles in common: %d',numel(allele_cap_)));
disp(sprintf(' %% number of alleles from Up05: %d/%d',numel(ij_Up05_from_cap_),numel(dataset_{1+ndataset_Up05}.bim_name_)));
disp(sprintf(' %% number of alleles from Ap05: %d/%d',numel(ij_Ap05_from_cap_),numel(dataset_{1+ndataset_Ap05}.bim_name_)));
end;%if (flag_verbose);

pca_rank_use = 2; %<-- number of principal-components to compute. ;

%%%%%%%%;
% First calculate pca of D_Up05_p01. ;
%%%%%%%%;
ndataset=ndataset_Up05;
tmp_parameter = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter.dir_0in: %s',tmp_parameter.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter.str_prefix: %s',tmp_parameter.str_prefix)); end;
tmp_parameter.dir_code = dir_code;
tmp_parameter.maf_lo_threshold = 0.01;
tmp_parameter.maf_hi_threshold = 0.50;
tmp_parameter.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter.flag_verbose = 0;
tmp_parameter.flag_force_create = 0;
pca_rank = pca_rank_use;
pca_mr_A_={dataset_{1+ndataset}.mx__.mr_A_full_};
pca_mr_Z_={dataset_{1+ndataset}.mx__.mr_Z_full_};
tmp_fname_bim = sprintf('%s/%s_bim.ext',tmp_parameter.dir_0in,tmp_parameter.str_prefix);
if ~exist(tmp_fname_bim,'file'); disp(sprintf(' %% Warning, %s not found',tmp_fname_bim)); end;
if  exist(tmp_fname_bim,'file');
if (flag_verbose); disp(sprintf(' %% calculating pca_mc_A')); end;
tmp_t = tic();
pca_mc_A = mc_from_bim_ext_ver5(tmp_fname_bim,tmp_parameter.maf_lo_threshold,tmp_parameter.maf_hi_threshold);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% mc_from_bim_ext_ver5: %0.6fs',tmp_t)); end;
end;%if  exist(tmp_fname_bim,'file');
pca_str_infix='D_Up05_p01';
parameter_D_Up05_p01 = tmp_parameter;
clear tmp_parameter;
%%%%;
parameter_Up05_A_p_p01 = parameter_D_Up05_p01;
parameter_Up05_A_p_p01.str_name_s0000 = 'A_p_p01';
[ ...
 parameter_Up05_A_p_p01 ...
,str_Up05_A_p_p01 ...
] = ...
xxxcluster_fromdisk_uADZSZDA_A_p_from_mx_ver16( ...
 parameter_Up05_A_p_p01 ...
,pca_mr_A_ ...
,pca_mr_Z_ ...
,pca_mc_A ...
,pca_str_infix ...
);
A_p_p01_ = mda_read_r8(str_Up05_A_p_p01); 
if (flag_disp);
figure(1+nf);nf=nf+1;clf;figsml;
plot(A_p_p01_);
xlabel('pcol');ylabel('A_p','Interpreter','none');
title(str_Up05_A_p_p01,'Interpreter','none');
end;%if (flag_disp);
%%%%;
parameter_D_Up05_p01.str_A_p = str_Up05_A_p_p01;
parameter_D_Up05_p01.str_name_s0000 = 'pca_D_Up05_p01';
[ ...
 parameter_D_Up05_p01 ...
,AZnV_D_Up05_p01_ ...
,AnV_D_Up05_p01_ ...
,ZnV_D_Up05_p01_ ...
,V_D_Up05_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_D_Up05_p01 ...
,pca_rank ...
,pca_mr_A_ ...
,pca_mr_Z_ ...
,pca_mc_A ...
,pca_str_infix ...
);
%%%%;
tmp_index_A_ = efind(dataset_{1+ndataset}.mx__.mr_A_full_);
tmp_index_Z_ = efind(dataset_{1+ndataset}.mx__.mr_Z_full_);
mr_dvx_Up05_ = zeros(dataset_{1+ndataset}.n_patient,1);
mr_dvx_Up05_(1+tmp_index_A_) = 1;
mr_dvx_Up05_(1+tmp_index_Z_) = 0;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_use = 12;
scatter(AZnV_D_Up05_p01_(:,1+0),AZnV_D_Up05_p01_(:,1+1),markersize_use,mr_dvx_Up05_,'filled','MarkerEdgeColor','k');
axis equal; grid on;
title('AZnV_D_Up05_p01_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
%%%%;


%%%%%%%%;
% Separate the Up05 patients into 3 continents based on previous classification with isosplit5. ;
%%%%%%%%;
fname_continents = sprintf('%s/dir_Up05/mr_Up05_p01_continent_.txt',dir_trunk);
fid = fopen(fname_continents, 'r');
mr_Up05_p01_continent_ = fscanf(fid, '%d');
fclose(fid);
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_use = 12;
scatter(AZnV_D_Up05_p01_(:,1+0),AZnV_D_Up05_p01_(:,1+1),markersize_use,mr_dvx_Up05_ + 2*mr_Up05_p01_continent_,'filled','MarkerEdgeColor','k');
axis equal; grid on;
title('Up05_p01_continent_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;
%%%%;


%%%%%%%%;
% Now calculate pca of D_Up05_cap_Ap05_p01. ;
%%%%%%%%;
fname_fig = sprintf('%s/D_Up05_cap_Ap05_p01_scatter_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
markersize_use = 12;
p_row = 3; p_col= 5; np=0;
rseed = 0;
rng(rseed);
f_ = linspace(0,1,p_row*p_col); n_f = numel(f_);
%%%%%%%%;
for nf=0:n_f-1;
%%%%%%%%;
f = f_(1+nf);
tmp_parameter = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter.dir_0in: %s',tmp_parameter.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter.str_prefix: %s',tmp_parameter.str_prefix)); end;
tmp_parameter.dir_code = dir_code;
tmp_parameter.maf_lo_threshold = 0.01;
tmp_parameter.maf_hi_threshold = 0.50;
tmp_parameter.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter.flag_verbose = 0;
tmp_parameter.flag_force_create = 1; %<-- remake each time. ;
pca_rank = pca_rank_use;
pca_mr_A_={dataset_{1+ndataset}.mx__.mr_A_full_};
pca_mr_Z_={dataset_{1+ndataset}.mx__.mr_Z_full_};
tmp_fname_bim = sprintf('%s/%s_bim.ext',tmp_parameter.dir_0in,tmp_parameter.str_prefix);
pca_mc_A = mx__.mc_A_;
if ~exist(tmp_fname_bim,'file'); disp(sprintf(' %% %s not found, continuing',tmp_fname_bim)); end;
if  exist(tmp_fname_bim,'file');
if (flag_verbose); disp(sprintf(' %% calculating pca_mc_A')); end;
tmp_t = tic();
pca_mc_A = mc_from_bim_ext_ver5(tmp_fname_bim,tmp_parameter.maf_lo_threshold,tmp_parameter.maf_hi_threshold);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% mc_from_bim_ext_ver5: %0.6fs',tmp_t)); end;
end;%if  exist(tmp_fname_bim,'file');
tmp_index_ = randperm(numel(index_Up05_setminus_cap_),floor(f*numel(index_Up05_setminus_cap_))) - 1;
pca_mc_A(1+index_Up05_setminus_cap_(1+tmp_index_)) = 0;
pca_str_infix='D_Up05_cap_Ap05_p01';
parameter_D_Up05_cap_Ap05_p01 = tmp_parameter;
clear tmp_parameter;
%%%%;
parameter_D_Up05_cap_Ap05_p01.str_A_p = str_Up05_A_p_p01;
parameter_D_Up05_cap_Ap05_p01.str_name_s0000 = 'pca_D_Up05_cap_Ap05_p01';
[ ...
 parameter_D_Up05_cap_Ap05_p01 ...
,AZnV_D_Up05_cap_Ap05_p01_ ...
,AnV_D_Up05_cap_Ap05_p01_ ...
,ZnV_D_Up05_cap_Ap05_p01_ ...
,V_D_Up05_cap_Ap05_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_D_Up05_cap_Ap05_p01 ...
,pca_rank ...
,pca_mr_A_ ...
,pca_mr_Z_ ...
,pca_mc_A ...
,pca_str_infix ...
);
%%%%;
tmp_index_A_ = efind(dataset_{1+ndataset}.mx__.mr_A_full_);
tmp_index_Z_ = efind(dataset_{1+ndataset}.mx__.mr_Z_full_);
mr_dvx_Up05_ = zeros(dataset_{1+ndataset}.n_patient,1);
mr_dvx_Up05_(1+tmp_index_A_) = 1;
mr_dvx_Up05_(1+tmp_index_Z_) = 0;
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_D_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Up05_cap_Ap05_p01_(:,1+1),markersize_use,mr_Up05_p01_continent_,'filled','MarkerEdgeColor','k');
axis equal; grid on;
title(sprintf('%0.3f (%d)',f,sum(pca_mc_A))); xlabel('PC1'); ylabel('PC2');
drawnow();
%%%%%%%%;
end;%for nf=0:n_f-1;
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now rerun the final projection of D_Up05_cap_Ap05_p01. ;
% Make sure that this is restricted to the intersection. ;
%%%%%%%%;
tmp_parameter = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter.dir_0in: %s',tmp_parameter.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter.str_prefix: %s',tmp_parameter.str_prefix)); end;
tmp_parameter.dir_code = dir_code;
tmp_parameter.maf_lo_threshold = 0.01;
tmp_parameter.maf_hi_threshold = 0.50;
tmp_parameter.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter.flag_verbose = 0;
tmp_parameter.flag_force_create = 0;
pca_rank = pca_rank_use;
pca_mr_A_={dataset_{1+ndataset}.mx__.mr_A_full_};
pca_mr_Z_={dataset_{1+ndataset}.mx__.mr_Z_full_};
tmp_fname_bim = sprintf('%s/%s_bim.ext',tmp_parameter.dir_0in,tmp_parameter.str_prefix);
pca_mc_A = mx__.mc_A_;
if ~exist(tmp_fname_bim,'file'); disp(sprintf(' %% %s not found, continuing',tmp_fname_bim)); end;
if  exist(tmp_fname_bim,'file');
if (flag_verbose); disp(sprintf(' %% calculating pca_mc_A')); end;
tmp_t = tic();
pca_mc_A = mc_from_bim_ext_ver5(tmp_fname_bim,tmp_parameter.maf_lo_threshold,tmp_parameter.maf_hi_threshold);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% mc_from_bim_ext_ver5: %0.6fs',tmp_t)); end;
end;%if  exist(tmp_fname_bim,'file');
pca_mc_A(1+index_Up05_setminus_cap_) = 0;
pca_str_infix='D_Up05_cap_Ap05_p01';
parameter_D_Up05_cap_Ap05_p01 = tmp_parameter;
clear tmp_parameter;
%%%%;
parameter_D_Up05_cap_Ap05_p01.str_A_p = str_Up05_A_p_p01;
parameter_D_Up05_cap_Ap05_p01.str_name_s0000 = 'pca_D_Up05_cap_Ap05_p01';
[ ...
 parameter_D_Up05_cap_Ap05_p01 ...
,AZnV_D_Up05_cap_Ap05_p01_ ...
,AnV_D_Up05_cap_Ap05_p01_ ...
,ZnV_D_Up05_cap_Ap05_p01_ ...
,V_D_Up05_cap_Ap05_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_D_Up05_cap_Ap05_p01 ...
,pca_rank ...
,pca_mr_A_ ...
,pca_mr_Z_ ...
,pca_mc_A ...
,pca_str_infix ...
);
%%%%;
tmp_index_A_ = efind(dataset_{1+ndataset}.mx__.mr_A_full_);
tmp_index_Z_ = efind(dataset_{1+ndataset}.mx__.mr_Z_full_);
mr_dvx_Up05_ = zeros(dataset_{1+ndataset}.n_patient,1);
mr_dvx_Up05_(1+tmp_index_A_) = 1;
mr_dvx_Up05_(1+tmp_index_Z_) = 0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
scatter(AZnV_D_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Up05_cap_Ap05_p01_(:,1+1),markersize_use,mr_Up05_p01_continent_,'filled','MarkerEdgeColor','k');
axis equal; grid on;
title(sprintf('D_Up05_cap_Ap05_p01'),'Interpreter','none'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;
%%%%%%%%;
parameter_D_Up05_cap_Ap05_p01.pca_str_infix = pca_str_infix;

%%%%%%%%;
% Now project Ap05 data onto the same principal-components. ;
%%%%%%%%;
parameter_D_Up05_cap_Ap05_p01.str_V_ = sprintf('%s/dir_pca/dir_pca_mda/pca_D_%s_k%d_B44_V_.mda',parameter_D_Up05_cap_Ap05_p01.dir_out_s0000,parameter_D_Up05_cap_Ap05_p01.pca_str_infix,pca_rank);
if ~exist(parameter_D_Up05_cap_Ap05_p01.str_V_,'file');
disp(sprintf(' %% Warning, %s not found',parameter_D_Up05_cap_Ap05_p01.str_V_));
end;%if ~exist(parameter_D_Up05_cap_Ap05_p01.str_V_,'file');
D_Up05_cap_Ap05_p01_V_ = mda_read_r8(parameter_D_Up05_cap_Ap05_p01.str_V_);
if (flag_verbose); disp(sprintf(' %% size of intersection: %d',numel(index_Up05_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% Support of V_: %d',sum(abs(D_Up05_cap_Ap05_p01_V_(:,1))>0))); end;
if (flag_verbose); disp(sprintf(' %% intersection of supp(V_) with cap: %d',numel(intersect(index_Up05_from_cap_,efind(abs(D_Up05_cap_Ap05_p01_V_(:,1))>0))))); end;
D_Ap05_from_Up05_cap_Ap05_p01_V_ = ij_Ap05_from_Up05__*D_Up05_cap_Ap05_p01_V_;
if (flag_verbose); disp(sprintf(' %% fnorm(D_Ap05_from_Up05_cap_Ap05_p01_V_(ij_Ap05_from_cap_,:) - D_Up05_cap_Ap05_p01_V_(ij_Up05_from_cap_,:)): %0.16f',fnorm(D_Ap05_from_Up05_cap_Ap05_p01_V_(ij_Ap05_from_cap_,:) - D_Up05_cap_Ap05_p01_V_(ij_Up05_from_cap_,:)))); end;

%%%%%%%%;
% Now calculate pca of D_Ap05_p01. ;
%%%%%%%%;
ndataset=ndataset_Ap05;
tmp_parameter = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter.dir_0in: %s',tmp_parameter.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter.str_prefix: %s',tmp_parameter.str_prefix)); end;
tmp_parameter.dir_code = dir_code;
tmp_parameter.maf_lo_threshold = 0.01;
tmp_parameter.maf_hi_threshold = 0.50;
tmp_parameter.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter.flag_verbose = 0;
tmp_parameter.flag_force_create = 0;
pca_rank = pca_rank_use;
pca_mr_A_={dataset_{1+ndataset}.mx__.mr_A_full_};
pca_mr_Z_={dataset_{1+ndataset}.mx__.mr_Z_full_};
tmp_fname_bim = sprintf('%s/%s_bim.ext',tmp_parameter.dir_0in,tmp_parameter.str_prefix);
if ~exist(tmp_fname_bim,'file'); disp(sprintf(' %% Warning, %s not found',tmp_fname_bim)); end;
if  exist(tmp_fname_bim,'file');
if (flag_verbose); disp(sprintf(' %% calculating pca_mc_A')); end;
tmp_t = tic();
pca_mc_A = mc_from_bim_ext_ver5(tmp_fname_bim,tmp_parameter.maf_lo_threshold,tmp_parameter.maf_hi_threshold);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% mc_from_bim_ext_ver5: %0.6fs',tmp_t)); end;
end;%if  exist(tmp_fname_bim,'file');
pca_str_infix='D_Ap05_p01';
parameter_D_Ap05_p01 = tmp_parameter;
clear tmp_parameter;
%%%%;
parameter_Ap05_A_p_p01 = parameter_D_Ap05_p01;
parameter_Ap05_A_p_p01.str_name_s0000 = 'A_p_p01';
[ ...
 parameter_Ap05_A_p_p01 ...
,str_Ap05_A_p_p01 ...
] = ...
xxxcluster_fromdisk_uADZSZDA_A_p_from_mx_ver16( ...
 parameter_Ap05_A_p_p01 ...
,pca_mr_A_ ...
,pca_mr_Z_ ...
,pca_mc_A ...
,pca_str_infix ...
);
A_p_p01_ = mda_read_r8(str_Ap05_A_p_p01); 
if (flag_disp);
figure(1+nf);nf=nf+1;clf;figsml;
plot(A_p_p01_);
xlabel('pcol');ylabel('A_p','Interpreter','none');
title(str_Ap05_A_p_p01,'Interpreter','none');
end;%if (flag_disp);
%%%%;
parameter_D_Ap05_p01.str_A_p = str_Ap05_A_p_p01;
parameter_D_Ap05_p01.str_name_s0000 = 'pca_D_Ap05_p01';
[ ...
 parameter_D_Ap05_p01 ...
,AZnV_D_Ap05_p01_ ...
,AnV_D_Ap05_p01_ ...
,ZnV_D_Ap05_p01_ ...
,V_D_Ap05_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_D_Ap05_p01 ...
,pca_rank ...
,pca_mr_A_ ...
,pca_mr_Z_ ...
,pca_mc_A ...
,pca_str_infix ...
);
%%%%;
tmp_index_A_ = efind(dataset_{1+ndataset}.mx__.mr_A_full_);
tmp_index_Z_ = efind(dataset_{1+ndataset}.mx__.mr_Z_full_);
mr_dvx_Ap05_ = zeros(dataset_{1+ndataset}.n_patient,1);
mr_dvx_Ap05_(1+tmp_index_A_) = 1;
mr_dvx_Ap05_(1+tmp_index_Z_) = 0;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_use = 12;
scatter(AZnV_D_Ap05_p01_(:,1+0),AZnV_D_Ap05_p01_(:,1+1),markersize_use,mr_dvx_Ap05_,'filled','MarkerEdgeColor','k');
axis equal; grid on;
title('AZnV_D_Ap05_p01_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
%%%%;

%%%%%%%%;
% Now project Ap05 data onto D_Ap05_from_Up05_cap_Ap05_p01_V_. ;
%%%%%%%%;
ndataset=ndataset_Ap05;
tmp_parameter = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter.dir_0in: %s',tmp_parameter.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter.str_prefix: %s',tmp_parameter.str_prefix)); end;
tmp_parameter.dir_code = dir_code;
tmp_parameter.maf_lo_threshold = 0.01;
tmp_parameter.maf_hi_threshold = 0.50;
tmp_parameter.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter.flag_verbose = 0;
tmp_parameter.flag_force_create = 1;
pca_rank = pca_rank_use;
pca_mr_A_={dataset_{1+ndataset}.mx__.mr_A_full_};
pca_mr_Z_={dataset_{1+ndataset}.mx__.mr_Z_full_};
tmp_fname_bim = sprintf('%s/%s_bim.ext',tmp_parameter.dir_0in,tmp_parameter.str_prefix);
if ~exist(tmp_fname_bim,'file'); disp(sprintf(' %% Warning, %s not found',tmp_fname_bim)); end;
if  exist(tmp_fname_bim,'file');
if (flag_verbose); disp(sprintf(' %% calculating pca_mc_A')); end;
tmp_t = tic();
pca_mc_A = mc_from_bim_ext_ver5(tmp_fname_bim,tmp_parameter.maf_lo_threshold,tmp_parameter.maf_hi_threshold);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% mc_from_bim_ext_ver5: %0.6fs',tmp_t)); end;
end;%if  exist(tmp_fname_bim,'file');
pca_str_infix='D_Ap05_from_Up05_cap_Ap05_p01';
parameter_D_Ap05_from_Up05_cap_Ap05_p01 = tmp_parameter;
clear tmp_parameter;
%%%%;
parameter_D_Ap05_from_Up05_cap_Ap05_p01.str_A_p = str_Ap05_A_p_p01;
parameter_D_Ap05_from_Up05_cap_Ap05_p01.str_name_s0000 = 'pca_D_Ap05_from_Up05_cap_Ap05_p01';
[ ...
 parameter_D_Ap05_from_Up05_cap_Ap05_p01 ...
,AZnV_D_Ap05_from_Up05_cap_Ap05_p01_ ...
,AnV_D_Ap05_from_Up05_cap_Ap05_p01_ ...
,ZnV_D_Ap05_from_Up05_cap_Ap05_p01_ ...
,tmp_V_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_proj_from_mx_ver16( ...
 parameter_D_Ap05_from_Up05_cap_Ap05_p01 ...
,pca_rank ...
,D_Ap05_from_Up05_cap_Ap05_p01_V_ ...
,pca_mr_A_ ...
,pca_mr_Z_ ...
,pca_mc_A ...
,pca_str_infix ...
);
%%%%;
tmp_index_A_ = efind(dataset_{1+ndataset}.mx__.mr_A_full_);
tmp_index_Z_ = efind(dataset_{1+ndataset}.mx__.mr_Z_full_);
mr_dvx_Ap05_ = zeros(dataset_{1+ndataset}.n_patient,1);
mr_dvx_Ap05_(1+tmp_index_A_) = 1;
mr_dvx_Ap05_(1+tmp_index_Z_) = 0;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_use = 12;
scatter(AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+1),markersize_use,mr_dvx_Ap05_,'filled','MarkerEdgeColor','k');
axis equal; grid on;
title('AZnV_D_Ap05_from_Up05_cap_Ap05_p01_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
%%%%;

%%%%%%%%;
fname_fig = sprintf('%s/D_Ap05_from_Up05_cap_Ap05_p01_scatter_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
markersize_use = 12;
np=0;
%%%%;
subplot(2,2,1+np);np=np+1;
scatter(AZnV_D_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Up05_cap_Ap05_p01_(:,1+1),markersize_use,mr_Up05_p01_continent_,'filled','MarkerEdgeColor',0.15*[1,1,1],'LineWidth',0.5);
axis equal; grid on;
title('AZnV_D_Up05_cap_Ap05_p01_ (by continent)','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
%%%%;
subplot(2,2,1+np);np=np+1;
scatter(AZnV_D_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Up05_cap_Ap05_p01_(:,1+1),markersize_use,mr_dvx_Up05_,'filled','MarkerEdgeColor',0.15*[1,1,1],'LineWidth',0.5);
axis equal; grid on;
title('AZnV_D_Up05_cap_Ap05_p01_ (case-vs-ctrl)','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
%%%%;
subplot(2,2,1+np);np=np+1;
scatter(AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+1),markersize_use,mr_dvx_Ap05_,'filled','MarkerEdgeColor',0.85*[0,1,0],'LineWidth',0.5);
axis equal; grid on;
title('AZnV_D_Ap05_from_Up05_cap_Ap05_p01_ (case-vs-ctrl)','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
%%%%;
subplot(2,2,1+np);np=np+1;
hold on;
scatter(AZnV_D_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Up05_cap_Ap05_p01_(:,1+1),markersize_use,mr_dvx_Up05_,'filled','MarkerEdgeColor',0.15*[1,1,1],'LineWidth',0.5);
scatter(AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+1),markersize_use,mr_dvx_Ap05_,'filled','MarkerEdgeColor',0.85*[0,1,0],'LineWidth',0.5);
hold off;
axis equal; grid on;
title('both together (case-vs-ctrl)','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now save the loadings. ;
%%%%%%%%;
flag_recalc = 0;
fname_ssv = sprintf('%s/dir_Up05/AZnV_D_Up05_p01_.txt',dir_trunk);
if (flag_recalc | ~exist(fname_ssv,'file')); save(fname_ssv,'AZnV_D_Up05_p01_','-ascii'); end;
fname_ssv = sprintf('%s/dir_Up05/AZnV_D_Up05_cap_Ap05_p01_.txt',dir_trunk);
if (flag_recalc | ~exist(fname_ssv,'file')); save(fname_ssv,'AZnV_D_Up05_cap_Ap05_p01_','-ascii'); end;
fname_ssv = sprintf('%s/dir_Up05/V_D_Up05_p01_.txt',dir_trunk);
if (flag_recalc | ~exist(fname_ssv,'file')); save(fname_ssv,'V_D_Up05_p01_','-ascii'); end;
fname_ssv = sprintf('%s/dir_Up05/V_D_Up05_cap_Ap05_p01_.txt',dir_trunk);
if (flag_recalc | ~exist(fname_ssv,'file')); save(fname_ssv,'V_D_Up05_cap_Ap05_p01_','-ascii'); end;

%%%%%%%%;
% use AZnV_D_Ap05_from_Up05_cap_Ap05_p01_ to separate the Ap05 patients into 3 continents; 
% Get mean of PC1 for each Up05 continent. Then assign Ap05 subjects based
% on closest value of PC1.
%%%%%%%%;
PC1_mean_Up05_cap_Ap05_p01 = grpstats(AZnV_D_Up05_cap_Ap05_p01_(:,1), mr_Up05_p01_continent_, 'mean');

V = AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1);
N = PC1_mean_Up05_cap_Ap05_p01;
A = repmat(N,[1 length(V)]);
[minValue,closestIndex] = min(abs(A-V'));
mr_Ap05_from_Up05_cap_Ap05_p01_continent_ = (closestIndex-1)';
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_use = 12;
scatter(AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+1),markersize_use,mr_dvx_Ap05_ + 2*mr_Ap05_from_Up05_cap_Ap05_p01_continent_,'filled','MarkerEdgeColor','k');
axis equal; grid on;
title('Ap05_from_Up05_cap_Ap05_p01_continent_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;
%%%%;
flag_recalc=0;
fname_tsv = sprintf('%s/dir_Ap05/mr_Ap05_from_Up05_cap_Ap05_p01_continent_.txt',dir_trunk);
if (flag_recalc | ~exist(fname_tsv,'file'));
fp = fopen(fname_tsv,'w');
fprintf(fp,'%d\n',mr_Ap05_from_Up05_cap_Ap05_p01_continent_);
fclose(fp);
end;%if (flag_recalc | ~exist(fname_tsv,'file'));
fname_ssv = sprintf('%s/dir_Ap05/AZnV_D_Ap05_from_Up05_cap_Ap05_p01_.txt',dir_trunk);
if (flag_recalc | ~exist(fname_ssv,'file')); save(fname_ssv,'AZnV_D_Ap05_from_Up05_cap_Ap05_p01_','-ascii'); end;


%%%%%%%%;
% Plot UKB and ADNI data together, do not color by continent;
%%%%%%%%;
fname_fig = sprintf('%s/UKB_vs_ADNI_pca',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
markersize_use = 16;
np=0;
hold on;
scatter(AZnV_D_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Up05_cap_Ap05_p01_(:,1+1),markersize_use,'filled','MarkerFaceColor','#0072B2','MarkerEdgeColor','k');
scatter(AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+1),markersize_use,'filled','MarkerFaceColor','#E69F00','MarkerEdgeColor','k');
legend('UKB','ADNI','Location','Best');
axisnotick; title('UKB vs ADNI'); xlabel('PC1'); ylabel('PC2');
hold off
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;

%%%%%%%%;
% Plot UKB and ADNI data together, color by continent;
%%%%%%%%;
fname_fig = sprintf('%s/UKB_vs_ADNI_pca_continents',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
markersize_use = 20;
np=0;
cmap = {'#D55E00', '#0072B2', '#009E73'};
cmap = validatecolor(cmap, 'multiple');
colormap(cmap);
hold on;
scatter(nan, nan,markersize_use+8, 's','MarkerEdgeColor','k','DisplayName','UKB');
scatter(nan, nan,markersize_use, 'o','MarkerEdgeColor','k','DisplayName','ADNI');
scatter(AZnV_D_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Up05_cap_Ap05_p01_(:,1+1),markersize_use+4,mr_Up05_p01_continent_,'s');
scatter(AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+0),AZnV_D_Ap05_from_Up05_cap_Ap05_p01_(:,1+1),markersize_use,mr_Ap05_from_Up05_cap_Ap05_p01_continent_,'o','filled','MarkerEdgeColor','k');
legend('UKB','ADNI','Location','Best','FontSize',14);
axisnotick; xlabel('PC1','FontSize',14); ylabel('PC2','FontSize',14);
hold off
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;
