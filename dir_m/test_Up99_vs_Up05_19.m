clear;
%%%%%%%%;
% Clustering to assign continent labels requires isosplit. This can be ;
% downloaded and installed from
% https://github.com/flatironinstitute/isosplit5 ;
%%%%%%%%;

%%%%%%%%;
% Here we set the matlab path and define str_home for naming. ;
%%%%%%%%;
% platform = 'rusty';
% if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
% if (strcmp(platform,'access1')); str_home = 'data'; end;
% if (strcmp(platform,'OptiPlex')); str_home = 'home'; end;
% if (strcmp(platform,'eval1')); str_home = 'home'; end;
% if (strcmp(platform,'rusty')); str_home = 'mnt/home'; end;
%%%%%%%%;
run('/home/jelman/Github/lakcluster_c/dir_m/setup_0'); %<-- set up the paths. ;
flag_verbose = 1;
flag_disp = 1+flag_verbose; nf=0;
flag_replot = 0;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Assume path is set using dir_lakcluster_c/dir_m/setup_0.m. ;')); end;

if (flag_verbose); disp(sprintf(' %% Comparing Up99 with Up05 data. ;')); end;
memory_GB = 128;
if (flag_verbose); disp(sprintf(' %% trying with memory set to %d ;', memory_GB));end;

dir_code = '/home/jelman/Github/lakcluster_c';
dir_trunk = '/home/jelman/Projects/AD_Biclustering/data/UKB/ukb_pca_p05-p1';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

n_dataset = 2; %<-- load two different datasets. ;
str_dataset_ = {'Up99','Up05'};
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
dataset.fname_bimext = sprintf('%s/test2mds_maf01_bim.ext',dataset.dir_0in);
end;%if ~exist(dataset.fname_bimext,'file');
if ~exist(dataset.fname_bimext,'file');
dataset.fname_bimext = sprintf('%s/%s_bim.ext',dataset.dir_0in,dataset.str_dataset);
end;%if ~exist(dataset.fname_bimext,'file');
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% load_bimext_ver1 ...')); end;
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
dataset.fname_famext = sprintf('%s/test2mds_maf01_fam.ext',dataset.dir_0in);
end;%if ~exist(dataset.fname_famext,'file');
if ~exist(dataset.fname_famext,'file');
dataset.fname_famext = sprintf('%s/%s_fam.ext',dataset.dir_0in,dataset.str_dataset);
end;%if ~exist(dataset.fname_famext,'file');
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% load_famext_ver1 ...')); end;
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
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% load_famext_ver1: %0.6fs',tmp_t)); end;
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
parameter.slurm_memdecl = memory_GB; %<--- May need to be increased ;
parameter.dir_0in = dataset_{1+ndataset}.dir_0in;
parameter.str_prefix = dataset_{1+ndataset}.str_prefix;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% load_mx__from_parameter_ver0 ...')); end;
mx__ = load_mx__from_parameter_ver0(parameter);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% load_mx__from_parameter_ver0: %0.6fs',tmp_t)); end;
dataset_{1+ndataset}.parameter = parameter;
dataset_{1+ndataset}.mx__ = mx__;
end;%for ndataset=0:n_dataset-1;

ndataset_Up99 = 0; %<-- index used for dataset Up99. ;
ndataset_Up05 = 1; %<-- index used for dataset Up05. ;
dir_trunk_Up99 = dataset_{1+ndataset_Up99}.dir_trunk;
dir_trunk_Up05 = dataset_{1+ndataset_Up05}.dir_trunk;
n_patient_Up99 = dataset_{1+ndataset_Up99}.n_patient;
n_patient_Up05 = dataset_{1+ndataset_Up05}.n_patient;
n_snp_Up99 = dataset_{1+ndataset_Up99}.n_snp;
n_snp_Up05 = dataset_{1+ndataset_Up05}.n_snp;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% intersect ...')); end;
[allele_cap_,ij_Up99_from_cap_,ij_Up05_from_cap_] = intersect(dataset_{1+ndataset_Up99}.bim_name_,dataset_{1+ndataset_Up05}.bim_name_,'stable');
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% intersect: %0.6fs',tmp_t)); end;
index_Up99_from_cap_ = ij_Up99_from_cap_ - 1;
index_Up05_from_cap_ = ij_Up05_from_cap_ - 1;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% setdiff ...')); end;
ij_Up99_from_Up05__ = sparse(ij_Up99_from_cap_,ij_Up05_from_cap_,1,n_snp_Up99,n_snp_Up05);
ij_Up05_from_Up99__ = sparse(ij_Up05_from_cap_,ij_Up99_from_cap_,1,n_snp_Up05,n_snp_Up99);
index_Up99_setminus_cap_ = setdiff(efind(dataset_{1+ndataset_Up99}.mx__.mc_A_),index_Up99_from_cap_);
index_Up05_setminus_cap_ = setdiff(efind(dataset_{1+ndataset_Up05}.mx__.mc_A_),index_Up05_from_cap_);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% setdiff: %0.6fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% numel(index_Up99_from_cap_): %d',numel(index_Up99_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% numel(index_Up05_from_cap_): %d',numel(index_Up05_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% numel(index_Up99_setminus_cap_): %d',numel(index_Up99_setminus_cap_))); end;
if (flag_verbose); disp(sprintf(' %% numel(index_Up05_setminus_cap_): %d',numel(index_Up05_setminus_cap_))); end;

pca_rank_use = 2; %<-- number of principal-components. ;
maf_lo_threshold = 0.01; %<-- minor-allele-frequency lower bound. ;

%%%%%%%%;
% Set Up99_A_p_. ;
%%%%%%%%;
ndataset=ndataset_Up99;
mx_tst__ = dataset_{1+ndataset}.mx__;
tmp_parameter_tst = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.dir_0in: %s',tmp_parameter_tst.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.str_prefix: %s',tmp_parameter_tst.str_prefix)); end;
tmp_parameter_tst.dir_code = dir_code;
tmp_parameter_tst.maf_lo_threshold = maf_lo_threshold;
tmp_parameter_tst.maf_hi_threshold = 0.50;
tmp_parameter_tst.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter_tst.flag_verbose = 0;
tmp_parameter_tst.flag_force_create = 0;
pca_rank = pca_rank_use;
mr_A_ori_tst_ = mx_tst__.mr_A_full_;
mr_Z_ori_tst_ = mx_tst__.mr_Z_full_;
mc_A_ori_tst_ = mx_tst__.mc_A_;
pca_str_infix_tst='D_trnUp99_tstUp05_nix_p01';
%%%%;
parameter_A_p_p01 = tmp_parameter_tst;
pca_A_p_str_infix = 'full';
parameter_A_p_p01.str_name_s0000 = 'A_p_p01';
parameter_A_p_p01.flag_force_create = 0;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_A_p_from_mx_ver16 ...')); end;
[ ...
 parameter_A_p_p01 ...
,str_Up99_A_p_p01 ...
] = ...
xxxcluster_fromdisk_uADZSZDA_A_p_from_mx_ver16( ...
parameter_A_p_p01 ...
,{ mr_A_ori_tst_ } ...
,{ mr_Z_ori_tst_ } ...
,mc_A_ori_tst_ ...
,pca_A_p_str_infix ...
);
tmp_Up99_A_p_p01_ = mda_read_r8(str_Up99_A_p_p01); 
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_A_p_from_mx_ver16: %0.6fs',tmp_t)); end;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot(tmp_Up99_A_p_p01_);
xlabel('pcol');ylabel('A_p','Interpreter','none');
title(str_Up99_A_p_p01,'Interpreter','none');
end;%if (flag_disp);

%%%%%%%%;
% Set Up05_A_p_. ;
%%%%%%%%;
ndataset=ndataset_Up05;
mx_tst__ = dataset_{1+ndataset}.mx__;
tmp_parameter_tst = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.dir_0in: %s',tmp_parameter_tst.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.str_prefix: %s',tmp_parameter_tst.str_prefix)); end;
tmp_parameter_tst.dir_code = dir_code;
tmp_parameter_tst.maf_lo_threshold = maf_lo_threshold;
tmp_parameter_tst.maf_hi_threshold = 0.50;
tmp_parameter_tst.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter_tst.flag_verbose = 0;
tmp_parameter_tst.flag_force_create = 0;
pca_rank = pca_rank_use;
mr_A_ori_tst_ = mx_tst__.mr_A_full_;
mr_Z_ori_tst_ = mx_tst__.mr_Z_full_;
mc_A_ori_tst_ = mx_tst__.mc_A_;
pca_str_infix_tst='D_Up05_from_trnUp99_tstUp05_nix_p01';
%%%%;
parameter_A_p_p01 = tmp_parameter_tst;
pca_A_p_str_infix = 'full';
parameter_A_p_p01.str_name_s0000 = 'A_p_p01';
parameter_A_p_p01.flag_force_create = 0;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_A_p_from_mx_ver16 ...')); end;
[ ...
 parameter_A_p_p01 ...
,str_Up05_A_p_p01 ...
] = ...
xxxcluster_fromdisk_uADZSZDA_A_p_from_mx_ver16( ...
parameter_A_p_p01 ...
,{ mr_A_ori_tst_ } ...
,{ mr_Z_ori_tst_ } ...
,mc_A_ori_tst_ ...
,pca_A_p_str_infix ...
);
tmp_Up05_A_p_p01_ = mda_read_r8(str_Up05_A_p_p01); 
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_A_p_from_mx_ver16: %0.6fs',tmp_t)); end;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot(tmp_Up05_A_p_p01_);
xlabel('pcol');ylabel('A_p','Interpreter','none');
title(str_Up05_A_p_p01,'Interpreter','none');
end;%if (flag_disp);

%%%%%%%%;
% if (flag_verbose); disp(sprintf(' %% ')); end;
% if (flag_verbose); disp(sprintf(' %% We divide Up99 data into three continents. ;')); end;
% % dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman',str_home);
% mr_Up99_p01_continent_ = textread(sprintf('%s/dir_Up99/mr_Up99_p01_continent_.txt',dir_trunk));
% n_continent = max(mr_Up99_p01_continent_)+1; %<-- continent index is zero-based. ;
% mr_Up99_p01_continent_pc__ = zeros(n_patient_Up99,n_continent);
% for ncontinent=0:n_continent-1;
% mr_Up99_p01_continent_pc__(:,1+ncontinent) = (mr_Up99_p01_continent_==ncontinent);
% end;%for ncontinent=0:n_continent-1;
% %%%%%%%%;
% if (flag_verbose); disp(sprintf(' %% ')); end;
% if (flag_verbose); disp(sprintf(' %% We divide Up05 data into three continents. ;')); end;
% % dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman',str_home);
% mr_Up05_p01_continent_ = textread(sprintf('%s/dir_Up05/mr_Up05_p01_continent_.txt',dir_trunk));
% n_continent = max(mr_Up05_p01_continent_)+1; %<-- continent index is zero-based. ;
% mr_Up05_p01_continent_pc__ = zeros(n_patient_Up05,n_continent);
% for ncontinent=0:n_continent-1;
% mr_Up05_p01_continent_pc__(:,1+ncontinent) = (mr_Up05_p01_continent_==ncontinent);
% end;%for ncontinent=0:n_continent-1;
%%%%%%%%;

%%%%%%%%;
% Now we load Kunkle_AD_GWAS_pvals.txt ;
%%%%%%%%;
fname_kunkle = sprintf('/home/jelman/Projects/AD_Biclustering/data/Kunkle_AD_GWAS_pvals.txt');
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% kunkle_ ...')); end;
fp = fopen(fname_kunkle,'r');
kunkle_ = textscan(fp,'%s %f','headerlines',1);
fclose(fp);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% kunkle_: %0.6fs',tmp_t)); end;
kunkle_vid_ = kunkle_{1};
kunkle_ADp_ = kunkle_{2};

%%%%%%%%;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% unique ...')); end;
[ ...
 dataset_{1+ndataset_Up99}.u_bim_vid_ ...
,dataset_{1+ndataset_Up99}.ij_nvid_from_nbim_ ...
,dataset_{1+ndataset_Up99}.ij_nbim_from_nvid_ ...
] = ...
unique(dataset_{1+ndataset_Up99}.bim_vid_,'stable');
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% unique: %0.6fs',tmp_t)); end;
%%%%%%%%;
dataset_{1+ndataset_Up99}.n_u_bim_vid = numel(dataset_{1+ndataset_Up99}.u_bim_vid_);
dataset_{1+ndataset_Up99}.ij_nbim_from_nvid__ = sparse(dataset_{1+ndataset_Up99}.ij_nbim_from_nvid_,1:dataset_{1+ndataset_Up99}.n_snp,1,dataset_{1+ndataset_Up99}.n_u_bim_vid,dataset_{1+ndataset_Up99}.n_snp);
%%%%%%%%;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% unique ...')); end;
[ ...
 dataset_{1+ndataset_Up05}.u_bim_vid_ ...
,dataset_{1+ndataset_Up05}.ij_nvid_from_nbim_ ...
,dataset_{1+ndataset_Up05}.ij_nbim_from_nvid_ ...
] = ...
unique(dataset_{1+ndataset_Up05}.bim_vid_,'stable');
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% unique: %0.6fs',tmp_t)); end;
%%%%%%%%;
dataset_{1+ndataset_Up05}.n_u_bim_vid = numel(dataset_{1+ndataset_Up05}.u_bim_vid_);
dataset_{1+ndataset_Up05}.ij_nbim_from_nvid__ = sparse(dataset_{1+ndataset_Up05}.ij_nbim_from_nvid_,1:dataset_{1+ndataset_Up05}.n_snp,1,dataset_{1+ndataset_Up05}.n_u_bim_vid,dataset_{1+ndataset_Up05}.n_snp);
%%%%%%%%;

tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% intersect ...')); end;
[cap_Up99_u_bim_vid_kunkle_vid_,ij_Up99_u_bim_vid_from_cap_,ij_kunkle_vid_from_cap_] = intersect(dataset_{1+ndataset_Up99}.u_bim_vid_,kunkle_vid_,'stable');
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% intersect: %0.6fs',tmp_t)); end;
cap_Up99_u_bim_vid_kunkle_vid_ADp_ = kunkle_ADp_(ij_kunkle_vid_from_cap_);
Up99_u_bim_vid_ADp_ = ones(dataset_{1+ndataset_Up99}.n_u_bim_vid,1);
Up99_u_bim_vid_ADp_(ij_Up99_u_bim_vid_from_cap_) = cap_Up99_u_bim_vid_kunkle_vid_ADp_;

n_test = 8;
for ntest=0:n_test-1;
n_cap = numel(cap_Up99_u_bim_vid_kunkle_vid_);
ncap = max(0,min(n_cap-1,floor(n_cap*ntest/n_test)));
tmp_kunkle_vid = kunkle_vid_{ij_kunkle_vid_from_cap_(1+ncap)};
tmp_u_bim_vid = dataset_{1+ndataset_Up99}.u_bim_vid_{ij_Up99_u_bim_vid_from_cap_(1+ncap)};
tmp_kunkle_ADp = kunkle_ADp_(ij_kunkle_vid_from_cap_(1+ncap));
tmp_Up99_u_bim_vid_ADp = Up99_u_bim_vid_ADp_(ij_Up99_u_bim_vid_from_cap_(1+ncap));
assert(strcmp(tmp_kunkle_vid,tmp_u_bim_vid));
assert(tmp_kunkle_ADp==tmp_Up99_u_bim_vid_ADp);
if (flag_verbose); disp(sprintf(' %% ntest %d/%d: %s vs %s, %0.4f vs %0.4f',ntest,n_test,tmp_kunkle_vid,tmp_u_bim_vid,tmp_kunkle_ADp,tmp_Up99_u_bim_vid_ADp)); end;
end;%for ntest=0:n_test-1;

Up99_bim_ADp_ = transpose(dataset_{1+ndataset_Up99}.ij_nbim_from_nvid__)*Up99_u_bim_vid_ADp_;

n_test = 8;
for ntest=0:n_test-1;
nsnp = max(0,min(dataset_{1+ndataset_Up99}.n_snp,floor(dataset_{1+ndataset_Up99}.n_snp*ntest/n_test)));
tmp_u_bim_vid = dataset_{1+ndataset_Up99}.bim_vid_{1+nsnp};
tmp_Up99_bim_ADp = Up99_bim_ADp_(1+nsnp);
if (tmp_Up99_bim_ADp< 1);
nk = efind(strcmp(kunkle_vid_,tmp_u_bim_vid));
tmp_kunkle_vid = kunkle_vid_{1+nk};
tmp_kunkle_ADp = kunkle_ADp_(1+nk);
if (flag_verbose); disp(sprintf(' %% ntest %d/%d: %s vs %s, %0.4f vs %0.4f',ntest,n_test,tmp_kunkle_vid,tmp_u_bim_vid,tmp_kunkle_ADp,tmp_Up99_bim_ADp)); end;
assert(strcmp(tmp_kunkle_vid,tmp_u_bim_vid));
assert(tmp_kunkle_ADp==tmp_Up99_bim_ADp);
end;%if (tmp_Up99_bim_ADp< 1);
end;%for ntest=0:n_test-1;

if (flag_verbose); disp(sprintf(' %% numel(intersect(dataset_{1+ndataset_Up99}.bim_vid_(1+efind(tmp_mc_)),dataset_{1+ndataset_Up05}.bim_vid_)): %d/%d',numel(intersect(dataset_{1+ndataset_Up99}.bim_vid_(1+efind(Up99_bim_ADp_<=0.05)),dataset_{1+ndataset_Up05}.bim_vid_)),numel(unique(dataset_{1+ndataset_Up05}.bim_vid_)))); end;

flag_calc = 1;
if flag_calc;
%%%%%%%%;
% Now calculate dominant principal-components across the entire Up99 data-set. ;
%%%%%%%%;
str_lak_vs_dex = 'dex';
str_prefix = 'test2mds_maf01';
gamma = 0.05;
n_mds_0in = 2; n_mds_repl = 1; ij_mds_use_ = [1:2];
parameter_Up99 = struct('type','parameter');
parameter_Up99.slurm_memdecl = memory_GB; %<-- this may need to be inreased
parameter_Up99.dir_code = dir_code;
parameter_Up99.dir_trunk = dir_trunk_Up99;
parameter_Up99.str_lak_vs_dex = str_lak_vs_dex;
parameter_Up99.str_prefix = str_prefix;
parameter_Up99.gamma = gamma;
parameter_Up99.n_mds = n_mds_0in;
parameter_Up99.n_mds_repl = n_mds_repl;
parameter_Up99.ij_mds_use_ = ij_mds_use_;
parameter_Up99.flag_force_create = 0; %<-- reload previous run. ;
parameter_Up99.flag_verbose = max(0,flag_verbose-1);
parameter_Up99.maf_lo_threshold = maf_lo_threshold;
parameter_Up99.n_shuffle = 4;
parameter_Up99.dir_0in = sprintf('%s/dir_%s',parameter_Up99.dir_trunk,parameter_Up99.str_prefix);
parameter_Up99.nshuffle = 0;
parameter_Up99 = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_Up99);
pca_rank = 2;
%%%%;
fname_mr_A_full = sprintf('%s/%s_mr_A_full.b16',parameter_Up99.dir_0in,parameter_Up99.str_prefix);
fname_mr_Z_full = sprintf('%s/%s_mr_Z_full.b16',parameter_Up99.dir_0in,parameter_Up99.str_prefix);
fname_mc_A = sprintf('%s/%s_mc_A.b16',parameter_Up99.dir_0in,parameter_Up99.str_prefix);
mr_A_ori_Up99_ = binary_uncompress(fname_mr_A_full)>0; %<-- thresholded to avoid negative values. ;
mr_Z_ori_Up99_ = binary_uncompress(fname_mr_Z_full)>0; %<-- thresholded to avoid negative values. ;
mc_A_ori_Up99_ = binary_uncompress(fname_mc_A)>0; %<-- thresholded to avoid negative values. ;
%%%%;
end;%if flag_calc;

flag_calc = 1;
%%%%%%%%;
% Now calculate dominant _DandX_ principal-components across the entire Up99 data-set. ;
%%%%%%%%;
if flag_calc;

pca_rank = 2;
p_threshold_ = 0.05:0.05:1.00; n_p_threshold = numel(p_threshold_);
% p_threshold_ = 0.60:0.005:0.65; n_p_threshold = numel(p_threshold_);
AZnV_DandX_Up99_p01_pnt___ = zeros(dataset_{1+ndataset_Up99}.n_patient,pca_rank,n_p_threshold);
%%%%%%%%;
for np_threshold=0:n_p_threshold-1;
%%%%%%%%;
p_threshold = p_threshold_(1+np_threshold);
str_p_threshold = sprintf('%.2d',min(99,floor(100*p_threshold)));
pca_mr_A_Up99_ = { 1*mr_A_ori_Up99_ + 1*mr_Z_ori_Up99_ };
pca_mr_Z_Up99_ = { 0*mr_A_ori_Up99_ + 0*mr_Z_ori_Up99_ };
pca_mc_A_Up99 = mc_A_ori_Up99_;
tmp_mc_ = zeros(dataset_{1+ndataset_Up99}.n_snp,1);
tmp_mc_(1+efind(Up99_bim_ADp_<=p_threshold))=1;
pca_mc_A_Up99 = tmp_mc_;
pca_str_infix_Up99=sprintf('Up99t%s_DandX_p01',str_p_threshold);
parameter_DandX_Up99_p01 = parameter_Up99;
parameter_DandX_Up99_p01.flag_force_create = 0;
parameter_DandX_Up99_p01.str_A_p = str_Up99_A_p_p01;
parameter_DandX_Up99_p01.str_name_s0000 = sprintf('pca_Up99t%s_DandX_p01',str_p_threshold);
parameter_DandX_Up99_p01.slurm_memdecl = memory_GB;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16 ...')); end;
[ ...
 parameter_DandX_Up99_p01 ...
,tmp_AZnV_DandX_Up99_p01_ ...
,tmp_AnV_DandX_Up99_p01_ ...
,tmp_ZnV_DandX_Up99_p01_ ...
,tmp_V_DandX_Up99_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_DandX_Up99_p01 ...
,pca_rank ...
,pca_mr_A_Up99_ ...
,pca_mr_Z_Up99_ ...
,pca_mc_A_Up99 ...
,pca_str_infix_Up99 ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16: %0.6fs',tmp_t)); end;
AZnV_DandX_Up99_p01_pnt___(:,:,1+np_threshold) = tmp_AZnV_DandX_Up99_p01_;
%%%%%%%%;
end;%for np_threshold=0:n_p_threshold-1;
%%%%%%%%;
end;%if flag_calc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot first two PCs and color by continent label at p<0.05 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath '/home/jelman/Github/isosplit5/matlab'
% Set filenames of case/control PCs
pca_cases_fname = sprintf('%s/dir_Up99/dir_test2mds_maf01_analyze/dir_test2mds_maf01_dex_p01_D_m2r1_g050/dir_pca/dir_pca_mda/pca_proj_D_Up99t05_DandX_p01_k2_B44_AnV_.mda',dir_trunk);
pca_ctrls_fname = sprintf('%s/dir_Up99/dir_test2mds_maf01_analyze/dir_test2mds_maf01_dex_p01_D_m2r1_g050/dir_pca/dir_pca_mda/pca_proj_D_Up99t05_DandX_p01_k2_B44_ZnV_.mda',dir_trunk);

% Get PCs of cases and add index
pca_cases = mda_read_r8(pca_cases_fname);
idx_cases = find(pca_cases(:,1));
pca_cases = pca_cases(idx_cases, :);
pca_cases = [pca_cases idx_cases ones(length(idx_cases),1)];
% Get PCs of controls and add index
pca_ctrls = mda_read_r8(pca_ctrls_fname);
idx_ctrls = find(pca_ctrls(:,1));
pca_ctrls = pca_ctrls(idx_ctrls, :);
pca_ctrls = [pca_ctrls idx_ctrls zeros(length(idx_ctrls),1)];
% Concat cases and controls, sort by idx
pca_all = [pca_cases; pca_ctrls];
pca_all = sortrows(pca_all, 3);

% Get cluster labels using isosplit
[labels_out_all] = isosplit5_mex(pca_all(:,1).');
% Name clusters based on size (e.g., cluster 1 has most individuals)
[clustcounts_n, clustcounts_idx] = sort(groupcounts(reshape(labels_out_all, [], 1)), 'descend');
orig_labels_out = labels_out_all;
labels_out_all = clustcounts_idx(orig_labels_out);

% % Plot PCs colored by case/control status
% gscatter(pca_all(:,1), pca_all(:,2), {labels_out_all, pca_all(:,4)})
% Set filename of fam file
fam_fname = dataset_{1}.fname_famext;
fam = readtable(fam_fname, 'Delimiter', '\t', 'FileType', 'delimitedtext');
% Concatenate subject ID, cluster number, case-control status
subjectClusters = [table2cell(fam(:, 2)) num2cell(labels_out_all) num2cell(pca_all(:,4)), num2cell(pca_all(:,1:2))];
% write out
clust_fname = sprintf('%s/subjectClusterAssignments_PC0.txt', dir_trunk);
writecell(subjectClusters,clust_fname,'Delimiter','\t');
% Write out 0-based continent labels to use when searching for second bicluster
writematrix(labels_out_all - 1, sprintf('%s/dir_Up99/mr_Up99_p01_continent_.txt',dir_trunk));
writematrix(labels_out_all - 1, sprintf('%s/dir_Up05/mr_Up05_p01_continent_.txt',dir_trunk));
% nf=0;


% platform = 'rusty';
% if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
% if (strcmp(platform,'access1')); str_home = 'data'; end;
% if (strcmp(platform,'OptiPlex')); str_home = 'home'; end;
% if (strcmp(platform,'eval1')); str_home = 'home'; end;
% if (strcmp(platform,'rusty')); str_home = 'mnt/home'; end;
% dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman',str_home);
% dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
% mr_Up99_p01_continent_ = textread(sprintf('%s/dir_Up99/mr_Up99_p01_continent_.txt',dir_trunk));
mr_Up99_p01_continent_ = labels_out_all;

fname_fig_pre = sprintf('%s/pca_proj_D_Up99txx_DandX_p01_k22_B44_AZnV_',dir_jpg); %<-- fname_fig_pre for coarse resolution of p-value-threshold. ;
%fname_fig_pre = sprintf('%s/pca_proj_D_Up99t6x_DandX_p01_k22_B44_AZnV_',dir_jpg); %<-- fname_fig_pre for finer resolution. ;
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
%%%%%%%%;
dir_pca_tmp = sprintf('%s/dir_Up99/dir_test2mds_maf01_analyze/dir_test2mds_maf01_dex_p01_D_m2r1_g050/dir_pca/dir_pca_mda',dir_trunk);
p_threshold_ = 0.05:0.05:1.00; n_p_threshold = numel(p_threshold_); %<-- coarse resolution of p-value-threshold. ;
%p_threshold_ = 0.60:0.005:0.65; n_p_threshold = numel(p_threshold_); %<-- finer resolution. ;
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
markersize_use = 12;
cmap = {'#D55E00', '#0072B2', '#009E73'};
cmap = validatecolor(cmap, 'multiple');
p_row = 4; p_col = ceil(n_p_threshold/p_row); np=0;
for np_threshold=0:n_p_threshold-1;
p_threshold = p_threshold_(1+np_threshold);
str_p_threshold = sprintf('%.2d',min(99,floor(100*p_threshold)));
n_snps = nnz(Up99_bim_ADp_<=p_threshold)/3;
subplot(p_row,p_col,1+np);np=np+1;
tmp_fname_AnV_txx = sprintf('%s/pca_proj_D_Up99t%s_DandX_p01_k2_B44_AnV_.mda',dir_pca_tmp,str_p_threshold);
tmp_fname_ZnV_txx = sprintf('%s/pca_proj_D_Up99t%s_DandX_p01_k2_B44_ZnV_.mda',dir_pca_tmp,str_p_threshold);
if exist(tmp_fname_AnV_txx,'file') & exist(tmp_fname_AnV_txx,'file');
tmp_AnV_txx__ = mda_read_r8(tmp_fname_AnV_txx);
tmp_ZnV_txx__ = mda_read_r8(tmp_fname_ZnV_txx);
tmp_AZnV_txx__ = tmp_AnV_txx__ + tmp_ZnV_txx__;
colormap(cmap);
scatter(tmp_AZnV_txx__(:,1+0),tmp_AZnV_txx__(:,1+1),markersize_use,mr_Up99_p01_continent_,'filled','MarkerEdgeColor','#D3D3D3','Linewidth',0.5,'MarkerEdgeAlpha',.5);
end;%if exist(tmp_fname_AnV_txx,'file') & exist(tmp_fname_AnV_txx,'file');
xlabel('PC1');ylabel('PC2');title(sprintf('p<%0.2f, n SNPs=%d',p_threshold, n_snps));
end;%for np_threshold=0:n_p_threshold-1;
%%%%%%%%;
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close gcf;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run PCA on random SNPs subset of the same size as p<0.05  
% Plot results colored by continent to determine whether    
% structure resembles p<0.05 or p<1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str_randp_threshold = 'random_p05';
pca_mr_A_Up99_ = { 1*mr_A_ori_Up99_ + 1*mr_Z_ori_Up99_ };
pca_mr_Z_Up99_ = { 0*mr_A_ori_Up99_ + 0*mr_Z_ori_Up99_ };
pca_mc_A_Up99 = mc_A_ori_Up99_;
tmp_mc_ = zeros(dataset_{1+ndataset_Up99}.n_snp,1);
n_snp_p05 = sum(Up99_bim_ADp_<=p_threshold);
idx_randp = randsample(length(tmp_mc_),n_snp_p05);
tmp_mc_(sort(idx_randp))=1;
pca_mc_A_Up99 = tmp_mc_;
pca_str_infix_Up99=sprintf('Up99t%s_DandX_p01',str_randp_threshold);
parameter_DandX_Up99_p01 = parameter_Up99;
parameter_DandX_Up99_p01.flag_force_create = 0;
parameter_DandX_Up99_p01.str_A_p = str_Up99_A_p_p01;
parameter_DandX_Up99_p01.str_name_s0000 = sprintf('pca_Up99t%s_DandX_p01',str_randp_threshold);
parameter_DandX_Up99_p01.slurm_memdecl = memory_GB;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16 ...')); end;
[ ...
parameter_DandX_Up99_p01 ...
,randp_AZnV_DandX_Up99_p01_ ...
,randp_AnV_DandX_Up99_p01_ ...
,randp_ZnV_DandX_Up99_p01_ ...
,randp_V_DandX_Up99_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
parameter_DandX_Up99_p01 ...
,pca_rank ...
,pca_mr_A_Up99_ ...
,pca_mr_Z_Up99_ ...
,pca_mc_A_Up99 ...
,pca_str_infix_Up99 ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16: %0.6fs',tmp_t)); end;

% Plot results of random subset of SNPs
dir_pca_tmp = sprintf('%s/dir_Up99/dir_test2mds_maf01_analyze/dir_test2mds_maf01_dex_p01_D_m2r1_g050/dir_pca/dir_pca_mda',dir_trunk);
markersize_use = 12;
cmap = {'#D55E00', '#0072B2', '#009E73'};
cmap = validatecolor(cmap, 'multiple');
randp_fname_AnV_txx = sprintf('%s/pca_proj_D_Up99t%s_DandX_p01_k2_B44_AnV_.mda',dir_pca_tmp,str_randp_threshold);
randp_fname_ZnV_txx = sprintf('%s/pca_proj_D_Up99t%s_DandX_p01_k2_B44_ZnV_.mda',dir_pca_tmp,str_randp_threshold);
randp_AnV_txx__ = mda_read_r8(randp_fname_AnV_txx);
randp_ZnV_txx__ = mda_read_r8(randp_fname_ZnV_txx);
randp_AZnV_txx__ = randp_AnV_txx__ + randp_ZnV_txx__;
colormap(cmap);
scatter(randp_AZnV_txx__(:,1+0),randp_AZnV_txx__(:,1+1),markersize_use,mr_Up99_p01_continent_,'filled','MarkerEdgeColor','k','Linewidth',.5);
xlabel('PC1');ylabel('PC2');
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% Plot p<0.05 SNPs with scatterplot and heatmap (full dataset and
% individual continents)
%%%%%%%%%%%%%%%%%%
parameter_scatterheat = struct('type','parameter_scatterheat');
parameter_scatterheat.min_prctile = 0; %<-- percentile to use for interval minimum (default 0) ; you can increase this to crop outliers.
parameter_scatterheat.max_prctile = 100; %<-- percentile to use for interval maximum (default 100) ; you can decrease this to crop outliers.
parameter_scatterheat.box_expand = 1.25; %<-- relative size of interior region to box margin (default 1.25) ; you can increase this to 'zoom out'.
parameter_scatterheat.h_DvX_scale = 0.015; %<-- colormap max/min (default 0.025); this isn't intelligently chosen; increase it to apply the heatmap over a larger range of histogram values.
parameter_scatterheat.n_h_0= 64; %<-- number of bins in x_0 direction (e.g., horizontal) (default 32+1); increase to refine.
parameter_scatterheat.n_h_1 = 32; %<-- number of bins in x_1 direction (e.g., vertical) (default 31+1); increase to refine.
parameter_scatterheat.n_tick_0 = 4; %<-- place a tick-mark every n_tick_0 bins in x_0 direction (default 2); increase to have fewer tick-marks.
parameter_scatterheat.n_tick_1 = 2; %<-- place a tick-mark every n_tick_1 bins in x_1 direction (default 2); increase to have fewer tick-marks.

label_Y_y_ = pca_all(:,4) + 1; %<-- 1=ctrls, 2=cases ;
p_threshold = 0.05;
str_p_threshold = sprintf('%.2d',min(99,floor(100*p_threshold)));
tmp_fname_AnV_txx = sprintf('%s/pca_proj_D_Up99t%s_DandX_p01_k2_B44_AnV_.mda',dir_pca_tmp,str_p_threshold);
tmp_fname_ZnV_txx = sprintf('%s/pca_proj_D_Up99t%s_DandX_p01_k2_B44_ZnV_.mda',dir_pca_tmp,str_p_threshold);
tmp_AnV_txx__ = mda_read_r8(tmp_fname_AnV_txx);
tmp_ZnV_txx__ = mda_read_r8(tmp_fname_ZnV_txx);
tmp_AZnV_txx__ = tmp_AnV_txx__ + tmp_ZnV_txx__;

% All continents
[~,AZnV_0_lim_,AZnV_0_tick_,AZnV_1_lim_,AZnV_1_tick_] = ...
test_scatter_and_heatmap_1( ...
 parameter_scatterheat ...
, label_Y_y_ ...
, tmp_AZnV_txx__ ...
);

% Continent 1
tmp_index_Up05_ = efind(tmp_AZnV_txx__(:,1)<-30);
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*3,768]);
[~,AZnV_0_lim_,AZnV_0_tick_,AZnV_1_lim_,AZnV_1_tick_] = ...
test_scatter_and_heatmap_1( ...
parameter_scatterheat ...
, label_Y_y_(1+tmp_index_Up05_) ...
, tmp_AZnV_txx__(1+tmp_index_Up05_,:) ...
, AZnV_0_lim_,AZnV_0_tick_,AZnV_1_lim_,AZnV_1_tick_);

% Continent 2
tmp_index_Up05_ = efind(tmp_AZnV_txx__(:,1)>67);
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,768]);
[~,AZnV_0_lim_,AZnV_0_tick_,AZnV_1_lim_,AZnV_1_tick_] = ...
test_scatter_and_heatmap_1( ...
parameter_scatterheat ...
, label_Y_y_(1+tmp_index_Up05_) ...
, tmp_AZnV_txx__(1+tmp_index_Up05_,:) ...
, AZnV_0_lim_,AZnV_0_tick_,AZnV_1_lim_,AZnV_1_tick_);

% Continent 3
tmp_index_Up05_ = efind(tmp_AZnV_txx__(:,2)>154)
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,768]);
[~,AZnV_0_lim_,AZnV_0_tick_,AZnV_1_lim_,AZnV_1_tick_] = ...
test_scatter_and_heatmap_1( ...
parameter_scatterheat ...
, label_Y_y_(1+tmp_index_Up05_) ...
, tmp_AZnV_txx__(1+tmp_index_Up05_,:)...
, AZnV_0_lim_,AZnV_0_tick_,AZnV_1_lim_,AZnV_1_tick_);




%%%%%%%%%%%%%%%%%%
% Plot PCA loadings of p<0.05 SNPs 
%%%%%%%%%%%%%%%%%%
% First get PCA from p<0.05 threshold;
p_threshold = 0.05
str_p_threshold = sprintf('%.2d',min(99,floor(100*p_threshold)));
pca_mr_A_Up99_ = { 1*mr_A_ori_Up99_ + 1*mr_Z_ori_Up99_ };
pca_mr_Z_Up99_ = { 0*mr_A_ori_Up99_ + 0*mr_Z_ori_Up99_ };
pca_mc_A_Up99 = mc_A_ori_Up99_;
tmp_mc_ = zeros(dataset_{1+ndataset_Up99}.n_snp,1);
tmp_mc_(1+efind(Up99_bim_ADp_<=p_threshold))=1;
pca_mc_A_Up99 = tmp_mc_;
pca_str_infix_Up99=sprintf('Up99t%s_DandX_p01',str_p_threshold);
parameter_DandX_Up99_p01 = parameter_Up99;
parameter_DandX_Up99_p01.flag_force_create = 0;
parameter_DandX_Up99_p01.str_A_p = str_Up99_A_p_p01;
parameter_DandX_Up99_p01.str_name_s0000 = sprintf('pca_Up99t%s_DandX_p01',str_p_threshold);
parameter_DandX_Up99_p01.slurm_memdecl = memory_GB;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16 ...')); end;
[ ...
 parameter_DandX_Up99_p01 ...
,tmp_AZnV_DandX_Up99_p01_ ...
,tmp_AnV_DandX_Up99_p01_ ...
,tmp_ZnV_DandX_Up99_p01_ ...
,tmp_V_DandX_Up99_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_DandX_Up99_p01 ...
,pca_rank ...
,pca_mr_A_Up99_ ...
,pca_mr_Z_Up99_ ...
,pca_mc_A_Up99 ...
,pca_str_infix_Up99 ...
);

% Concatenate with chromosome, position, and maf
col_loadings_p05 = [tmp_V_DandX_Up99_p01_ double(dataset_{1}.bim_khr_) double(dataset_{1}.bim_pdi_) dataset_{1}.bim_maf_];
% Sort by chr and pos
col_loadings_p05 = sortrows(col_loadings_p05, [3 4]);
% Filter for only p<0.05 rows. These have non-zero values in first 2 cols;
nonZeroRows = ~(col_loadings_p05(:,1) == 0 & col_loadings_p05(:,2) == 0);
col_loadings_p05 = col_loadings_p05(nonZeroRows, :);
% Add headers: PC1, PC2, chr, pos, maf. The save to file
fname_loadings_p05 = '/home/jelman/Projects/AD_Biclustering/data/UKB/ukb_pca_p05-p1/ukb_pca_p05_snp-loadings.txt';
header = {'PC1', 'PC2', 'chr', 'pos', 'maf'};
col_loadings_p05 = [header; num2cell(col_loadings_p05)];
writecell(col_loadings_p05, fname_loadings_p05, 'Delimiter', '\t');

% Add row indices as columm
rowNums = double(1:size(col_loadings_p05))';
col_loadings_p05 = [col_loadings_p05 rowNums];
% Plot PC1
binscatter(col_loadings_p05(:,6), col_loadings_p05(:,1))
xlabel('Column #');
ylabel('PC1 loading');

% Plot PC2
binscatter(col_loadings_p05(:,6), col_loadings_p05(:,2))
xlabel('Column #');
ylabel('PC2 loading');

max_loadings = [col_loadings_p05 abs(col_loadings_p05(:,1))];
max_loadings = sortrows(max_loadings, -7);
% Top loadings all fall in a region of chr17. Subset these rows;
chr17_loadings = max_loadings(1:7715,:);
% Get boundaries of LD region;
minpos = min(chr17_loadings(:,4));
maxpos = max(chr17_loadings(:,4));