clear;
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

if (flag_verbose); disp(sprintf(' %% Comparing Up99 with Ap99 data. ;')); end;
if (flag_verbose); disp(sprintf(' %% Note that much of this analysis involves ncontinent = 0. ;')); end;

memory_GB = 128;
dir_code = '/home/jelman/Github/lakcluster_c';
dir_trunk = '/home/jelman/Projects/AD_Biclustering/data/ADNI/ADNI_vs_UKB';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

%%%%%%%%;
% first we load the 2 datasets. ;
%%%%%%%%;
n_dataset = 2;
str_dataset_ = {'Up99','Ap99'};
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

%%%%%%%%;
% Now we measure the dataset sizes. ;
%%%%%%%%;
for nd=0:1;
if (flag_verbose); disp(sprintf(' %% nd %d <-- %s',nd,dataset_{1+nd}.dir_0in)); end;
[~,n_row,n_col] = binary_getsize(sprintf('%s/%s_A_full_n.b16',dataset_{1+nd}.dir_0in,dataset_{1+nd}.str_prefix));
if (flag_verbose); disp(sprintf(' %% %% A_full_n_.b16 [%d,%d]',n_row,n_col)); end;
if (flag_verbose); disp(sprintf(' %% %% n_patient [%d]',dataset_{1+nd}.n_patient)); end;
if (flag_verbose); disp(sprintf(' %% %% n_snp [%d]',dataset_{1+nd}.n_snp)); end;
end;%for nd=0:1;

%%%%%%%%;
% Now we load the row- and col-masks for the datasets. ;
%%%%%%%%;
for ndataset=0:n_dataset-1;
parameter = struct('type','parameter');
parameter.slurm_memdecl = memory_GB; %<--- May need to be increased ;
parameter.dir_0in = dataset_{1+ndataset}.dir_0in;
parameter.str_prefix = dataset_{1+ndataset}.str_prefix;
mx__ = load_mx__from_parameter_ver0(parameter);
dataset_{1+ndataset}.parameter = parameter;
dataset_{1+ndataset}.mx__ = mx__;
end;%for ndataset=0:n_dataset-1;

%%%%%%%%;
% Now we define indices used to refer to the datasets. ;
% These will be used below to restrict the snp-subsets to the intersection of the two datasets. ;
% Here 'allele_cap' refers to the intersection of allele-combinations between the two datasets. ;
% index_Up99_from_cap_ cross-references the allele_cap and the Up99 dataset. ;
% index_Ap99_from_cap_ cross-references the allele_cap and the Ap99 dataset. ;
% index_Up99_setminus_cap_ lists those allele-combinations in Up99 which are not part of the intersection. ;
% index_Ap99_setminus_cap_ lists those allele-combinations in Ap99 which are not part of the intersection. ;
%%%%%%%%;
ndataset_Up99 = 0;
ndataset_Ap99 = 1;
n_patient_Up99 = dataset_{1+ndataset_Up99}.n_patient;
n_patient_Ap99 = dataset_{1+ndataset_Ap99}.n_patient;
n_snp_Up99 = dataset_{1+ndataset_Up99}.n_snp;
n_snp_Ap99 = dataset_{1+ndataset_Ap99}.n_snp;
[allele_cap_,ij_Up99_from_cap_,ij_Ap99_from_cap_] = intersect(dataset_{1+ndataset_Up99}.bim_name_,dataset_{1+ndataset_Ap99}.bim_name_,'stable');
index_Up99_from_cap_ = ij_Up99_from_cap_ - 1;
index_Ap99_from_cap_ = ij_Ap99_from_cap_ - 1;
ij_Up99_from_Ap99__ = sparse(ij_Up99_from_cap_,ij_Ap99_from_cap_,1,n_snp_Up99,n_snp_Ap99);
ij_Ap99_from_Up99__ = sparse(ij_Ap99_from_cap_,ij_Up99_from_cap_,1,n_snp_Ap99,n_snp_Up99);
index_Up99_setminus_cap_ = setdiff(efind(dataset_{1+ndataset_Up99}.mx__.mc_A_),index_Up99_from_cap_);
index_Ap99_setminus_cap_ = setdiff(efind(dataset_{1+ndataset_Ap99}.mx__.mc_A_),index_Ap99_from_cap_);
if (flag_verbose); disp(sprintf(' %% numel(index_Up99_from_cap_): %d',numel(index_Up99_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% numel(index_Ap99_from_cap_): %d',numel(index_Ap99_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% numel(index_Up99_setminus_cap_): %d',numel(index_Up99_setminus_cap_))); end;
if (flag_verbose); disp(sprintf(' %% numel(index_Ap99_setminus_cap_): %d',numel(index_Ap99_setminus_cap_))); end;

%%%%%%%%;
% The number of pcs we will calculate is pca_rank_use. ;
% The minor-allele-frequency we will use is maf_lo_threshold. ;
%%%%%%%%;
pca_rank_use = 2;
maf_lo_threshold = 0.01;

%%%%%%%%;
% Over the course of the next few paragraphs we load the Kunkle_AD_GWAS_pvals. ;
% These will be used to determine the disease-related p-value for each snp. ;
%%%%%%%%;

%%%%%%%%;
% First we load Kunkle_AD_GWAS_pvals.txt ;
% kunkle_vid_ is the variant id, ;
% kunkle_ADp_ is the AD-p-value. ;
%%%%%%%%;
fname_kunkle = '/home/jelman/Projects/AD_Biclustering/data/Kunkle_AD_GWAS_pvals.txt';
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% kunkle_ ...')); end;
fp = fopen(fname_kunkle,'r');
kunkle_ = textscan(fp,'%s %f','headerlines',1);
fclose(fp);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% kunkle_: %0.6fs',tmp_t)); end;
kunkle_vid_ = kunkle_{1};
kunkle_ADp_ = kunkle_{2};

%%%%%%%%;
% Now we cross-reference the snp indices. ;
% (this requires some bookkeeping). ;
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
 dataset_{1+ndataset_Ap99}.u_bim_vid_ ...
,dataset_{1+ndataset_Ap99}.ij_nvid_from_nbim_ ...
,dataset_{1+ndataset_Ap99}.ij_nbim_from_nvid_ ...
] = ...
unique(dataset_{1+ndataset_Ap99}.bim_vid_,'stable');
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% unique: %0.6fs',tmp_t)); end;
%%%%%%%%;
dataset_{1+ndataset_Ap99}.n_u_bim_vid = numel(dataset_{1+ndataset_Ap99}.u_bim_vid_);
dataset_{1+ndataset_Ap99}.ij_nbim_from_nvid__ = sparse(dataset_{1+ndataset_Ap99}.ij_nbim_from_nvid_,1:dataset_{1+ndataset_Ap99}.n_snp,1,dataset_{1+ndataset_Ap99}.n_u_bim_vid,dataset_{1+ndataset_Ap99}.n_snp);
%%%%%%%%;

%%%%%%%%;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% intersect ...')); end;
[cap_Up99_u_bim_vid_kunkle_vid_,ij_Up99_u_bim_vid_from_cap_,ij_Up99_kunkle_vid_from_cap_] = intersect(dataset_{1+ndataset_Up99}.u_bim_vid_,kunkle_vid_,'stable');
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% intersect: %0.6fs',tmp_t)); end;
cap_Up99_u_bim_vid_kunkle_vid_ADp_ = kunkle_ADp_(ij_Up99_kunkle_vid_from_cap_);
Up99_u_bim_vid_ADp_ = ones(dataset_{1+ndataset_Up99}.n_u_bim_vid,1);
Up99_u_bim_vid_ADp_(ij_Up99_u_bim_vid_from_cap_) = cap_Up99_u_bim_vid_kunkle_vid_ADp_;
%%%%%%%%;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% intersect ...')); end;
[cap_Ap99_u_bim_vid_kunkle_vid_,ij_Ap99_u_bim_vid_from_cap_,ij_Ap99_kunkle_vid_from_cap_] = intersect(dataset_{1+ndataset_Ap99}.u_bim_vid_,kunkle_vid_,'stable');
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% intersect: %0.6fs',tmp_t)); end;
cap_Ap99_u_bim_vid_kunkle_vid_ADp_ = kunkle_ADp_(ij_Ap99_kunkle_vid_from_cap_);
Ap99_u_bim_vid_ADp_ = ones(dataset_{1+ndataset_Ap99}.n_u_bim_vid,1);
Ap99_u_bim_vid_ADp_(ij_Ap99_u_bim_vid_from_cap_) = cap_Ap99_u_bim_vid_kunkle_vid_ADp_;
%%%%%%%%;

n_test = 8;
for ntest=0:n_test-1;
n_cap = numel(cap_Up99_u_bim_vid_kunkle_vid_);
ncap = max(0,min(n_cap-1,floor(n_cap*ntest/n_test)));
tmp_kunkle_vid = kunkle_vid_{ij_Up99_kunkle_vid_from_cap_(1+ncap)};
tmp_u_bim_vid = dataset_{1+ndataset_Up99}.u_bim_vid_{ij_Up99_u_bim_vid_from_cap_(1+ncap)};
tmp_kunkle_ADp = kunkle_ADp_(ij_Up99_kunkle_vid_from_cap_(1+ncap));
tmp_Up99_u_bim_vid_ADp = Up99_u_bim_vid_ADp_(ij_Up99_u_bim_vid_from_cap_(1+ncap));
assert(strcmp(tmp_kunkle_vid,tmp_u_bim_vid));
assert(tmp_kunkle_ADp==tmp_Up99_u_bim_vid_ADp);
if (flag_verbose); disp(sprintf(' %% ntest %d/%d: %s vs %s, %0.4f vs %0.4f',ntest,n_test,tmp_kunkle_vid,tmp_u_bim_vid,tmp_kunkle_ADp,tmp_Up99_u_bim_vid_ADp)); end;
end;%for ntest=0:n_test-1;

Up99_bim_ADp_ = transpose(dataset_{1+ndataset_Up99}.ij_nbim_from_nvid__)*Up99_u_bim_vid_ADp_;
Ap99_bim_ADp_ = transpose(dataset_{1+ndataset_Ap99}.ij_nbim_from_nvid__)*Ap99_u_bim_vid_ADp_;

%%%%%%%%;
% And now we finally have a usable array Up99_bim_ADp_, where: ;
% Up99_bim_ADp_(1+nsnp) is the AD-p-value associated with dataset_{1+ndataset_Up99}.bim_vid_{1+nsnp}. ;
%%%%%%%%;
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
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% numel(intersect(dataset_{1+ndataset_Up99}.bim_vid_(1+efind(Up99_bim_ADp_<=0.05)),dataset_{1+ndataset_Up99}.bim_vid_)): %d/%d',numel(intersect(dataset_{1+ndataset_Up99}.bim_vid_(1+efind(Up99_bim_ADp_<=0.05)),dataset_{1+ndataset_Up99}.bim_vid_)),numel(unique(dataset_{1+ndataset_Up99}.bim_vid_)))); end;
%%%%%%%%;

%%%%%%%%;
% Similarly, we have a usable array Ap99_bim_ADp_, where: ;
% Ap99_bim_ADp_(1+nsnp) is the AD-p-value associated with dataset_{1+ndataset_Ap99}.bim_vid_{1+nsnp}. ;
%%%%%%%%;
n_test = 8;
for ntest=0:n_test-1;
nsnp = max(0,min(dataset_{1+ndataset_Ap99}.n_snp,floor(dataset_{1+ndataset_Ap99}.n_snp*ntest/n_test)));
tmp_u_bim_vid = dataset_{1+ndataset_Ap99}.bim_vid_{1+nsnp};
tmp_Ap99_bim_ADp = Ap99_bim_ADp_(1+nsnp);
if (tmp_Ap99_bim_ADp< 1);
nk = efind(strcmp(kunkle_vid_,tmp_u_bim_vid));
tmp_kunkle_vid = kunkle_vid_{1+nk};
tmp_kunkle_ADp = kunkle_ADp_(1+nk);
if (flag_verbose); disp(sprintf(' %% ntest %d/%d: %s vs %s, %0.4f vs %0.4f',ntest,n_test,tmp_kunkle_vid,tmp_u_bim_vid,tmp_kunkle_ADp,tmp_Ap99_bim_ADp)); end;
assert(strcmp(tmp_kunkle_vid,tmp_u_bim_vid));
assert(tmp_kunkle_ADp==tmp_Ap99_bim_ADp);
end;%if (tmp_Ap99_bim_ADp< 1);
end;%for ntest=0:n_test-1;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% numel(intersect(dataset_{1+ndataset_Ap99}.bim_vid_(1+efind(Ap99_bim_ADp_<=0.05)),dataset_{1+ndataset_Ap99}.bim_vid_)): %d/%d',numel(intersect(dataset_{1+ndataset_Ap99}.bim_vid_(1+efind(Ap99_bim_ADp_<=0.05)),dataset_{1+ndataset_Ap99}.bim_vid_)),numel(unique(dataset_{1+ndataset_Ap99}.bim_vid_)))); end;
%%%%%%%%;

%%%%%%%%;
% Over the course of the next two paragraphs we define Up99_A_p_ and Ap99_A_p_. ;
% These are the 'baseline' values expected for each allele-combination in each dataset. ;
%%%%%%%%;

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
tmp_parameter_tst.slurm_memdecl = memory_GB; %<-- this may need to be inreased
tmp_parameter_tst.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter_tst.flag_verbose = 0;
tmp_parameter_tst.flag_force_create = 0;
pca_rank = pca_rank_use;
mr_A_ori_tst_ = mx_tst__.mr_A_full_;
mr_Z_ori_tst_ = mx_tst__.mr_Z_full_;
mc_A_ori_tst_ = mx_tst__.mc_A_;
pca_str_infix_tst='D_trnUp99_tstAp99_nix_p01';
%%%%;
parameter_A_p_p01 = tmp_parameter_tst;
pca_A_p_str_infix = 'full';
parameter_A_p_p01.str_name_s0000 = 'A_p_p01';
parameter_A_p_p01.flag_force_create = 0;
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
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot(tmp_Up99_A_p_p01_);
xlabel('pcol');ylabel('A_p','Interpreter','none');
title(str_Up99_A_p_p01,'Interpreter','none');
end;%if (flag_disp);

%%%%%%%%;
% Set Ap99_A_p_. ;
%%%%%%%%;
ndataset=ndataset_Ap99;
mx_tst__ = dataset_{1+ndataset}.mx__;
tmp_parameter_tst = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.dir_0in: %s',tmp_parameter_tst.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.str_prefix: %s',tmp_parameter_tst.str_prefix)); end;
tmp_parameter_tst.dir_code = dir_code;
tmp_parameter_tst.maf_lo_threshold = maf_lo_threshold;
tmp_parameter_tst.maf_hi_threshold = 0.50;
tmp_parameter_tst.slurm_memdecl = memory_GB; %<-- this may need to be inreased
tmp_parameter_tst.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter_tst.flag_verbose = 0;
tmp_parameter_tst.flag_force_create = 0;
pca_rank = pca_rank_use;
mr_A_ori_tst_ = mx_tst__.mr_A_full_;
mr_Z_ori_tst_ = mx_tst__.mr_Z_full_;
mc_A_ori_tst_ = mx_tst__.mc_A_;
pca_str_infix_tst='D_Ap99_from_trnUp99_tstAp99_nix_p01';
%%%%;
parameter_A_p_p01 = tmp_parameter_tst;
pca_A_p_str_infix = 'full';
parameter_A_p_p01.str_name_s0000 = 'A_p_p01';
parameter_A_p_p01.flag_force_create = 0;
[ ...
 parameter_A_p_p01 ...
,str_Ap99_A_p_p01 ...
] = ...
xxxcluster_fromdisk_uADZSZDA_A_p_from_mx_ver16( ...
parameter_A_p_p01 ...
,{ mr_A_ori_tst_ } ...
,{ mr_Z_ori_tst_ } ...
,mc_A_ori_tst_ ...
,pca_A_p_str_infix ...
);
tmp_Ap99_A_p_p01_ = mda_read_r8(str_Ap99_A_p_p01); 
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot(tmp_Ap99_A_p_p01_);
xlabel('pcol');ylabel('A_p','Interpreter','none');
title(str_Ap99_A_p_p01,'Interpreter','none');
end;%if (flag_disp);

%%%%%%%%;
% Now we calculate the principal-component projections of Up99 ;
% using the snps in allele_cap_. ;
% Note that this is not the projection that we want to show when ;
% analyzing Up99 initially, since it involves the replication ;
% dataset Ap99 (i.e., Ap99 helps determine what snps are in allele_cap_). ;
% So these plots are probably only something we should show after ;
% discussing replication. ;
% To indicate this distinction, we name these particular analyses using ;
% infixes of the form 'trnUp99_tstAp99' and 'Ap99_from_trnUp99_tstAp99'. ;
%%%%%%%%;
% Note that the particular datasets we are using here only allow for ;
% p_threshold to range from 0.00 to 0.05, as this is the ;
% disease-related p-value threshold for Up99 and Ap99. ;
% For a larger data-set we should increase the p_threshold range. ;
%%%%%%%%;
str_lak_vs_dex = 'dex';
str_prefix = 'test2mds_maf01';
gamma = 0.05;
n_mds_0in = 2; n_mds_repl = 1; ij_mds_use_ = [1:2];
%%%%;
dir_trunk_Up99 = sprintf('%s/dir_Up99',dir_trunk);
mx_Up99__ = dataset_{1+ndataset_Up99}.mx__;
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
parameter_Up99.slurm_memdecl = memory_GB; %<-- this may need to be inreased
parameter_Up99.n_shuffle = 128;
parameter_Up99.dir_0in = sprintf('%s/dir_%s',parameter_Up99.dir_trunk,parameter_Up99.str_prefix);
parameter_Up99.nshuffle = 0;
%%%%;
dir_trunk_Ap99 = sprintf('%s/dir_Ap99',dir_trunk);
mx_Ap99__ = dataset_{1+ndataset_Ap99}.mx__;
parameter_Ap99.dir_code = dir_code;
parameter_Ap99.dir_trunk = dir_trunk_Ap99;
parameter_Ap99.str_lak_vs_dex = str_lak_vs_dex;
parameter_Ap99.str_prefix = str_prefix;
parameter_Ap99.gamma = gamma;
parameter_Ap99.n_mds = n_mds_0in;
parameter_Ap99.n_mds_repl = n_mds_repl;
parameter_Ap99.ij_mds_use_ = ij_mds_use_;
parameter_Ap99.flag_force_create = 0; %<-- reload previous run. ;
parameter_Ap99.flag_verbose = max(0,flag_verbose-1);
parameter_Ap99.maf_lo_threshold = maf_lo_threshold;
parameter_Ap99.slurm_memdecl = memory_GB; %<-- this may need to be inreased
parameter_Ap99.n_shuffle = 128;
parameter_Ap99.dir_0in = sprintf('%s/dir_%s',parameter_Ap99.dir_trunk,parameter_Ap99.str_prefix);
parameter_Ap99.nshuffle = 0;
%%%%%%%%;
pca_rank = 2;
p_threshold_ = 0.05:0.05:1.00; n_p_threshold = numel(p_threshold_); %<-- This is the range for Up99 and Ap99. ;
% p_threshold_ = 0.01:0.01:0.05; n_p_threshold = numel(p_threshold_); %<-- This is the range for Up05 and Ap05. ;
AZnV_DandX_trnUp99_tstAp99_txx_p01_pnt___ = zeros(dataset_{1+ndataset_Up99}.n_patient,pca_rank,n_p_threshold);
AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_p01_pnt___ = zeros(dataset_{1+ndataset_Ap99}.n_patient,pca_rank,n_p_threshold);
for np_threshold=0:n_p_threshold-1;
%%%%%%%%;
p_threshold = p_threshold_(1+np_threshold);
str_p_threshold = sprintf('%.2d',min(99,floor(100*p_threshold)));
pca_mr_A_Up99_ = { 1*mx_Up99__.mr_A_full_ + 1*mx_Up99__.mr_Z_full_ };
pca_mr_Z_Up99_ = { 0*mx_Up99__.mr_A_full_ + 0*mx_Up99__.mr_Z_full_ };
pca_mc_A_Up99 = mx_Up99__.mc_A_;
tmp_mc_ = zeros(dataset_{1+ndataset_Up99}.n_snp,1);
tmp_mc_(1+efind(Up99_bim_ADp_<=p_threshold))=1;
pca_mc_A_Up99 = pca_mc_A_Up99.*tmp_mc_;
pca_mc_A_Up99(1+index_Up99_setminus_cap_) = 0; %<-- now turn off the allele-combinations not in allele_cap_. ;
pca_str_infix_Up99=sprintf('trnUp99_tstAp99_t%s_DandX_p01',str_p_threshold);
parameter_DandX_trnUp99_tstAp99_txx_p01 = parameter_Up99;
parameter_DandX_trnUp99_tstAp99_txx_p01.flag_force_create = 0;
parameter_DandX_trnUp99_tstAp99_txx_p01.str_A_p = str_Up99_A_p_p01;
parameter_DandX_trnUp99_tstAp99_txx_p01.str_name_s0000 = sprintf('pca_trnUp99_tstAp99_t%s_DandX_p01',str_p_threshold);
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16 ...')); end;
[ ...
 parameter_DandX_trnUp99_tstAp99_txx_p01 ...
,tmp_AZnV_DandX_trnUp99_tstAp99_txx_p01_ ...
,tmp_AnV_DandX_trnUp99_tstAp99_txx_p01_ ...
,tmp_ZnV_DandX_trnUp99_tstAp99_txx_p01_ ...
,tmp_V_DandX_trnUp99_tstAp99_txx_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_DandX_trnUp99_tstAp99_txx_p01 ...
,pca_rank ...
,pca_mr_A_Up99_ ...
,pca_mr_Z_Up99_ ...
,pca_mc_A_Up99 ...
,pca_str_infix_Up99 ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16: %0.6fs',tmp_t)); end;
AZnV_DandX_trnUp99_tstAp99_txx_p01_pnt___(:,:,1+np_threshold) = tmp_AZnV_DandX_trnUp99_tstAp99_txx_p01_;
%%%%%%%%;
% Now at this stage parameter_DandX_trnUp99_tstAp99_txx_p01 has the appropriate paths. ;
% This allows us to use these paths to project Ap99 data onto the principal-components defined using the Up99 data. ;
% Note that we have to use those principal-components defined using the intersection ;
% of the Up99 data with the allele_cap (see pca_str_infix_Ap99 later on below). ;
%%%%%%%%;
pca_mr_A_Ap99_ = { 1*mx_Ap99__.mr_A_full_ + 1*mx_Ap99__.mr_Z_full_ };
pca_mr_Z_Ap99_ = { 0*mx_Ap99__.mr_A_full_ + 0*mx_Ap99__.mr_Z_full_ };
pca_mc_A_Ap99 = mx_Ap99__.mc_A_;
tmp_mc_ = zeros(dataset_{1+ndataset_Ap99}.n_snp,1);
tmp_mc_(1+efind(Ap99_bim_ADp_<=p_threshold))=1;
pca_mc_A_Ap99 = pca_mc_A_Ap99.*tmp_mc_;
pca_mc_A_Ap99(1+index_Ap99_setminus_cap_) = 0; %<-- now turn off the allele-combinations not in allele_cap_. ;
%pca_str_infix_Ap99=sprintf('Ap99_from_trnUp99_tstAp99_t%s_DandX_p01',str_p_threshold); %<-- unused, see below. ;
parameter_DandX_Ap99_from_trnUp99_tstAp99_txx_p01 = parameter_Ap99;
parameter_DandX_Ap99_from_trnUp99_tstAp99_txx_p01.flag_force_create = 0;
parameter_DandX_Ap99_from_trnUp99_tstAp99_txx_p01.str_A_p = str_Ap99_A_p_p01;
parameter_DandX_Ap99_from_trnUp99_tstAp99_txx_p01.str_name_s0000 = sprintf('pca_Ap99_from_trnUp99_tstAp99_t%s_DandX_p01',str_p_threshold);
parameter_DandX_Ap99_from_trnUp99_tstAp99_txx_p01.slurm_memdcl = memory_GB;
%%%%;
pca_str_infix_Ap99=sprintf('trnUp99_tstAp99_t%s_DandX_p01',str_p_threshold); %<-- this defines the pca file to use. ;
parameter_DandX_trnUp99_tstAp99_txx_p01.str_V_ = sprintf('%s/dir_pca/dir_pca_mda/pca_D_%s_k%d_B44_V_.mda',parameter_DandX_trnUp99_tstAp99_txx_p01.dir_out_s0000,pca_str_infix_Ap99,pca_rank);
if ~exist(parameter_DandX_trnUp99_tstAp99_txx_p01.str_V_,'file');
disp(sprintf(' %% Warning, %s not found',parameter_DandX_trnUp99_tstAp99_txx_p01.str_V_));
end;%if ~exist(parameter_DandX_trnUp99_tstAp99_txx_p01.str_V_,'file');
DandX_trnUp99_tstAp99_txx_p01_V_ = mda_read_r8(parameter_DandX_trnUp99_tstAp99_txx_p01.str_V_);
if (flag_verbose); disp(sprintf(' %% size of intersection: %d',numel(index_Up99_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% Support of V_: %d',sum(abs(DandX_trnUp99_tstAp99_txx_p01_V_(:,1))>0))); end;
if (flag_verbose); disp(sprintf(' %% intersection of supp(V_) with cap: %d',numel(intersect(index_Up99_from_cap_,efind(abs(DandX_trnUp99_tstAp99_txx_p01_V_(:,1))>0))))); end;
DandX_Ap99_from_trnUp99_tstAp99_txx_p01_V_ = ij_Ap99_from_Up99__*DandX_trnUp99_tstAp99_txx_p01_V_;
if (flag_verbose); disp(sprintf(' %% fnorm(DandX_Ap99_from_trnUp99_tstAp99_txx_p01_V_(ij_Ap99_from_cap_,:) - DandX_trnUp99_tstAp99_txx_p01_V_(ij_Up99_from_cap_,:)): %0.16f',fnorm(DandX_Ap99_from_trnUp99_tstAp99_txx_p01_V_(ij_Ap99_from_cap_,:) - DandX_trnUp99_tstAp99_txx_p01_V_(ij_Up99_from_cap_,:)))); end;
%%%%;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_proj_from_mx_ver16 ...')); end;
[ ...
 parameter_DandX_Ap99_from_trnUp99_tstAp99_txx_p01 ...
,tmp_AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_ ...
,tmp_AnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_ ...
,tmp_ZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_ ...
,tmp_V_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_proj_from_mx_ver16( ...
 parameter_DandX_Ap99_from_trnUp99_tstAp99_txx_p01 ...
,pca_rank ...
,DandX_Ap99_from_trnUp99_tstAp99_txx_p01_V_ ...
,pca_mr_A_Ap99_ ...
,pca_mr_Z_Ap99_ ...
,pca_mc_A_Ap99 ...
,pca_str_infix_Ap99 ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_proj_from_mx_ver16: %0.6fs',tmp_t)); end;
AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_pnt___(:,:,1+np_threshold) = tmp_AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_;
%%%%%%%%;
end;%for np_threshold=0:n_p_threshold-1;
%%%%%%%%;

%%%%%%%%;
% Read in continent labels for Up05 and Ap05. ;
%%%%%%%%;
mr_Up05_p01_continent_ = textread(sprintf('%s/dir_Up05/mr_Up05_p01_continent_.txt',dir_trunk));
mr_Up99_p01_continent_ = mr_Up05_p01_continent_ ; %<--Up05 and Up99 labels are identical ;
% One subject present in the Ap05 dataset is filtered out of the Ap99 ;
% dataset due to missingness filter. Create new continent label array that;
% removes this subjects ;
mr_Ap05_p01_continent_ = textread(sprintf('%s/dir_Ap05/mr_Ap05_from_Up05_cap_Ap05_p01_continent_.txt',dir_trunk));
Ap99_fam_iid_ = dataset_{2}.fam_iid_ ;
Ap05_fam_fname = strrep(dataset_{2}.fname_famext,'Ap99','Ap05') ;
[ ...
,Ap05_n_patient ...
,Ap05_fam_fid_ ...
,Ap05_fam_iid_ ...
,Ap05_fam_yid_ ...
,Ap05_.fam_xid_ ...
,Ap05_fam_sex_ ...
,Ap05_fam_dvx_ ...
,Ap05_fam_dir_ ...
,Ap05_fam_fidandiid_ ...
,~ ...
  ] = ...
load_famext_ver1( ...
 Ap05_fam_fname ...
);
[~,idx_Ap99_in_Ap05] = ismember(Ap99_fam_iid_,Ap05_fam_iid_) ; 
mr_Ap99_p01_continent_ = mr_Ap05_p01_continent_(idx_Ap99_in_Ap05) ;


%%%%%%%%;
% Now put together an array of scatterplots. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
markersize_use = 4;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
mr_dvx_trnUp99_tstAp99_ = 0.0*mx_Up99__.mr_A_full_;
mr_dvx_trnUp99_tstAp99_(1+efind(mx_Up99__.mr_Z_full_)) = 1; %<-- color ctrl-patients. ;
mr_dvx_trnUp99_tstAp99_(1+efind(mx_Up99__.mr_A_full_)) = 2; %<-- color case-patients. ;
mr_dvx_Ap99_from_trnUp99_tstAp99_ = 0.0*mx_Ap99__.mr_A_full_;
mr_dvx_Ap99_from_trnUp99_tstAp99_(1+efind(mx_Ap99__.mr_Z_full_)) = 1; %<-- color ctrl-patients. ;
mr_dvx_Ap99_from_trnUp99_tstAp99_(1+efind(mx_Ap99__.mr_A_full_)) = 2; %<-- color case-patients. ;
p_row = 4; p_col = n_p_threshold;
for np_threshold=0:n_p_threshold-1;
pcol = np_threshold;
p_threshold = p_threshold_(1+np_threshold);
str_p_threshold = sprintf('%.2d',min(99,floor(100*p_threshold)));
AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_pn__ = AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_pnt___(:,:,1+np_threshold);
AZnV_DandX_trnUp99_tstAp99_txx_p01_pn__ = AZnV_DandX_trnUp99_tstAp99_txx_p01_pnt___(:,:,1+np_threshold);
AZnV_0_min = min(min(AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_pn__(:,1+0)),min(AZnV_DandX_trnUp99_tstAp99_txx_p01_pn__(:,1+0)));
AZnV_0_max = max(max(AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_pn__(:,1+0)),max(AZnV_DandX_trnUp99_tstAp99_txx_p01_pn__(:,1+0)));
AZnV_0_lim_ = [AZnV_0_min,AZnV_0_max]; AZnV_0_lim_ = mean(AZnV_0_lim_) + 1.25*0.5*diff(AZnV_0_lim_)*[-1,+1];
AZnV_1_min = min(min(AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_pn__(:,1+1)),min(AZnV_DandX_trnUp99_tstAp99_txx_p01_pn__(:,1+1)));
AZnV_1_max = max(max(AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_pn__(:,1+1)),max(AZnV_DandX_trnUp99_tstAp99_txx_p01_pn__(:,1+1)));
AZnV_1_lim_ = [AZnV_1_min,AZnV_1_max]; AZnV_1_lim_ = mean(AZnV_1_lim_) + 1.25*0.5*diff(AZnV_1_lim_)*[-1,+1];
for prow=0:p_row-1;
if prow==0; tmp_pn__ = AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_pn__; tmp_str_0 = 'Ap99'; end;
if prow==1; tmp_pn__ = AZnV_DandX_Ap99_from_trnUp99_tstAp99_txx_p01_pn__; tmp_str_0 = 'Ap99'; end;
if prow==2; tmp_pn__ = AZnV_DandX_trnUp99_tstAp99_txx_p01_pn__; tmp_str_0 = 'Up99'; end;
if prow==3; tmp_pn__ = AZnV_DandX_trnUp99_tstAp99_txx_p01_pn__; tmp_str_0 = 'Up99'; end;
if prow==0; tmp_p_ = mr_dvx_Ap99_from_trnUp99_tstAp99_; tmp_str_1 = 'DvX'; tmp_p_min = 1; tmp_p_max = 2; end;
if prow==1; tmp_p_ = mr_Ap99_p01_continent_; tmp_str_1 = 'continent'; tmp_p_min = 0; tmp_p_max = 2; end;
if prow==2; tmp_p_ = mr_dvx_trnUp99_tstAp99_; tmp_str_1 = 'DvX'; tmp_p_min = 1; tmp_p_max = 2; end;
if prow==3; tmp_p_ = mr_Up99_p01_continent_; tmp_str_1 = 'continent'; tmp_p_min = 0; tmp_p_max = 2; end;
if prow==0; tmp_legend_ = {'case','ctrl'}; end;
if prow==1; tmp_legend_ = {'cont0','cont1','cont2'}; end;
if prow==2; tmp_legend_ = {'case','ctrl'}; end;
if prow==3; tmp_legend_ = {'cont0','cont1','cont2'}; end;
subplot(p_row,p_col,1+pcol+prow*p_col);cla;
hold on;
for tmp_p=tmp_p_min:tmp_p_max;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*(tmp_p-tmp_p_min)/max(1e-12,tmp_p_max-tmp_p_min))));
tmp_index_ = efind(tmp_p_==tmp_p);
plot(tmp_pn__(1+tmp_index_,1+0),tmp_pn__(1+tmp_index_,1+1),'o','MarkerSize',markersize_use,'MarkerEdgeColor','none','MarkerFaceColor',c_80s__(1+nc_80s,:));
end;%for tmp_p=tmp_p_min:tmp_p_max;
hold off;
xlim(AZnV_0_lim_); xlabel('pc0'); ylim(AZnV_1_lim_); ylabel('pc1'); grid on;
legend(tmp_legend_,'Location','NorthWest');
title(sprintf('%s_t%s_%s',tmp_str_0,str_p_threshold,tmp_str_1),'Interpreter','none');
end;%for prow=0:p_row-1;
end;%for np_threshold=0:n_p_threshold-1;
%%%%%%%%;
fname_fig_pre = sprintf('%s/AZnV_xxxx_trnUp99_tstAp99_txx_DandX_p01_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre); fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
%close(gcf);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;
