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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now we continue with our usual business, ;
% Finding biclusters for specific continents, etc. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Here we load the continent-labels for the datasets. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% We divide Up99 data into three continents. ;')); end;
mr_Up99_p01_continent_ = textread(sprintf('%s/dir_Up99/mr_Up99_p01_continent_.txt',dir_trunk));
n_continent = max(mr_Up99_p01_continent_)+1;
mr_Up99_p01_continent_pc__ = zeros(n_patient_Up99,n_continent);
for ncontinent=0:n_continent-1;
mr_Up99_p01_continent_pc__(:,1+ncontinent) = (mr_Up99_p01_continent_==ncontinent);
end;%for ncontinent=0:n_continent-1;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% We divide Ap99 data into three continents. ;')); end;
mr_Ap99_p01_continent_ = textread(sprintf('%s/dir_Ap99/mr_Ap99_from_0UKB_cap_Ap99_p01_continent_.txt',dir_trunk));
n_continent = max(mr_Ap99_p01_continent_)+1;
mr_Ap99_p01_continent_pc__ = zeros(n_patient_Ap99,n_continent);
for ncontinent=0:n_continent-1;
mr_Ap99_p01_continent_pc__(:,1+ncontinent) = (mr_Ap99_p01_continent_==ncontinent);
end;%for ncontinent=0:n_continent-1;
%%%%%%%%;

%%%%%%%%;
% Now we specify the biclustering run to use. ;
% To begin with, we limit ourselves to the training dataset Up99, and ncontinent=0. ;
% This continent will be used to define: ;
% mr_A_alt_trn_, the list of case-patients in that continent (within the training dataset Up99). ;
% mr_Z_alt_trn_, the list of ctrl-patients in that continent (within the training dataset Up99). ;
% Then, within that continent (and training dataset), we pick an iteration. ;
% The iteration we use is defined by the (1-based) 'tmp_ij_nlpR', ;
% which corresponds to the (0-based) tmp_index_nlpR = tmp_ij_nlpR-1. ;
% In this case we define this using a simple 'max' function. ;
% In a more general case one might consider defining tmp_index_nlpR using Z_imax_zerobased. ;
% This particular iteration is then used to define: ;
% index_r_rem_trn_, which lists the rows (patients) to retain in the training dataset. ;
% index_c_rem_trn_, which lists the cols (allele-combinations) to retain in the training dataset. ;
% These will be used below to define row- and col-masks for the pca. ;
%%%%%%%%;
dir_trunk_Up99 = sprintf('%s/dir_Up99',dir_trunk);
ncontinent = 0; %<-- here we set the continent to use below. ;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% maf_lo_threshold %0.2f; ncontinent %d/%d;',maf_lo_threshold,ncontinent,n_continent)); end;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% setting parameters. ;')); end;
str_lak_vs_dex = 'dex';
str_prefix = 'test2mds_maf01';
gamma = 0.05;
n_mds_0in = 2; n_mds_repl = 1; ij_mds_use_ = [1:2];
parameter_trn.dir_code = dir_code;
parameter_trn.dir_trunk = dir_trunk_Up99;
parameter_trn.str_lak_vs_dex = str_lak_vs_dex;
parameter_trn.str_prefix = str_prefix;
parameter_trn.gamma = gamma;
parameter_trn.n_mds = n_mds_0in;
parameter_trn.n_mds_repl = n_mds_repl;
parameter_trn.ij_mds_use_ = ij_mds_use_;
parameter_trn.flag_force_create = 0; %<-- reload previous run. ;
parameter_trn.flag_verbose = max(0,flag_verbose-1);
parameter_trn.maf_lo_threshold = maf_lo_threshold;
parameter_trn.slurm_memdcl = memory_GB;
parameter_trn.n_shuffle = 128;
parameter_trn.dir_0in = sprintf('%s/dir_%s',parameter_trn.dir_trunk,parameter_trn.str_prefix);
fname_mr_A_full = sprintf('%s/%s_mr_A_full.b16',parameter_trn.dir_0in,parameter_trn.str_prefix);
fname_mr_Z_full = sprintf('%s/%s_mr_Z_full.b16',parameter_trn.dir_0in,parameter_trn.str_prefix);
fname_mc_A = sprintf('%s/%s_mc_A.b16',parameter_trn.dir_0in,parameter_trn.str_prefix);
mr_A_ori_trn_ = binary_uncompress(fname_mr_A_full)>0; %<-- thresholded to avoid negative values. ;
mr_Z_ori_trn_ = binary_uncompress(fname_mr_Z_full)>0; %<-- thresholded to avoid negative values. ;
mc_A_ori_trn_ = binary_uncompress(fname_mc_A)>0; %<-- thresholded to avoid negative values. ;
assert(numel(mr_A_ori_trn_)==numel(mr_Up99_p01_continent_));
mr_A_alt_trn_ = mr_A_ori_trn_.*mr_Up99_p01_continent_pc__(:,1+ncontinent);
mr_Z_alt_trn_ = mr_Z_ori_trn_.*mr_Up99_p01_continent_pc__(:,1+ncontinent);
parameter_trn.str_mr_0in = sprintf('continent%d',ncontinent);
parameter_trn.nshuffle = 0;
parameter_trn = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_trn);
%%%%%%%%;
parameter_trn.str_out_xdrop_a_s0000 = sprintf('%s/out_xdrop_a.txt',parameter_trn.dir_out_s0000);
xdrop_trn_ = load_out_xdrop_from_str_ver0(parameter_trn.str_out_xdrop_a_s0000);
index_rdrop_trn_ = xdrop_trn_.index_rdrop_;
index_cdrop_trn_ = xdrop_trn_.index_cdrop_;
parameter_trn.dir_out_s0000_jpg = sprintf('%s/dir_out_jpg',parameter_trn.dir_out_s0000);
if ~exist(parameter_trn.dir_out_s0000_jpg,'dir'); disp(sprintf(' %% mkdir %s',parameter_trn.dir_out_s0000_jpg)); mkdir(parameter_trn.dir_out_s0000_jpg); end;
trace_trn____ = load_trace__from_dir_ver0(parameter_trn.dir_out_trace);
[tmp_nlpR,tmp_ij_nlpR] = max(trace_trn____.nlpR_s0000_);
r_rem_trn = trace_trn____.r_rem_s0000_(tmp_ij_nlpR);
c_rem_trn = trace_trn____.c_rem_s0000_(tmp_ij_nlpR);
index_r_rem_trn_ = xdrop_trn_.index_rkeep_(1:r_rem_trn);
index_c_rem_trn_ = xdrop_trn_.index_ckeep_(1:c_rem_trn);
assert(numel(index_r_rem_trn_)==numel(intersect(efind(mr_A_alt_trn_),index_r_rem_trn_)));
assert(numel(index_c_rem_trn_)==numel(intersect(efind(mc_A_ori_trn_),index_c_rem_trn_)));
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
linewidth_sml = 0.5;
linewidth_big = 3;
markersize_big = 16;
hold on;
plot(trace_trn____.niter_s0000_,trace_trn____.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(trace_trn____.niter_s0000_,trace_trn____.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(trace_trn____.niter_s0000_(tmp_ij_nlpR),trace_trn____.ZR_s0000_(tmp_ij_nlpR),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
hold off;
xlim([min(trace_trn____.niter_s0000_),max(trace_trn____.niter_s0000_)]); xlabel('iteration');
ylabel('negative-log-p');
end;%if flag_disp;
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;
figsml;fig80s;
p_row = 1; p_col = 1; np=0;
%%%%;
[A_p_c_,A_p_0_,AZ_rsum_] = load_A_p_c_from_dir_0(parameter_trn.dir_out_s0000);
alpha_c_ = A_p_c_ - (1-A_p_c_);
D_c_ = sqrt(1./max(0.01,4.0*A_p_c_.*(1-A_p_c_)));
mx_trn__ = load_mx__from_parameter_ver0(parameter_trn);
mr_dvx_trn_ = 0.0*mx_trn__.mr_A_full_;
mr_dvx_trn_(1+efind(mx_trn__.mr_A_full_))=2;
mr_dvx_trn_(1+efind(mx_trn__.mr_Z_full_))=1;
markersize_sml = 12;
markersize_big = 16;
linewidth_use = 1.5;
%%%%;
ni=tmp_ij_nlpR-1;
[ ...
 AZnV_nix_driver_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_ni_ver16( ...
 parameter_trn ...
,ni ...
,xdrop_trn_ ...
,trace_trn____ ...
,mx_trn__ ...
);
subplot(p_row,p_col,1+np);np=np+1;
hold on;
scatter(AZnV_nix_driver_(1+efind(mr_dvx_trn_),1),AZnV_nix_driver_(1+efind(mr_dvx_trn_),2),markersize_sml,mr_dvx_trn_(1+efind(mr_dvx_trn_)),'filled','MarkerEdgeColor','k');
scatter(AZnV_nix_driver_(1+index_r_rem_trn_,1),AZnV_nix_driver_(1+index_r_rem_trn_,2),markersize_big,mr_dvx_trn_(1+index_r_rem_trn_),'filled','MarkerEdgeColor','k','LineWidth',linewidth_use,'MarkerEdgeColor',0.85*[1,1,1]);
hold off;
axisnotick; title('bicluster-informed (driver)'); xlabel('PC1'); ylabel('PC2');
%%%%;
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Now that we have picked out a particular iteration (above), ;
% we can calculate the pca for that iteration directly. ;
% In this case we deliberately set the col-mask (i.e., pca_mc_A_trn) ;
% to exclude those elements of Up99 which are not in the allele_cap ;
% defined above (i.e., the allele-combination-intersection of Up99 and Ap99). ;
%%%%%%%%;
tmp_parameter_trn = parameter_trn;
pca_rank = pca_rank_use;
mr_A_rem_trn_ = 0.0*mr_A_alt_trn_; mr_A_rem_trn_(1+index_r_rem_trn_) = 1;
assert(sum(mr_A_rem_trn_)==sum(mr_A_rem_trn_.*mr_A_alt_trn_)); %<-- ensure that we only turned on rows in ncontinent. ;
mr_Z_rem_trn_ = mr_Z_alt_trn_; %<-- limit the ctrls to those in ncontinent. ;
mc_A_rem_trn_ = 0.0*mc_A_ori_trn_; mc_A_rem_trn_(1+index_c_rem_trn_) = 1;
%%%%%%%%;
% Note that mc_A_rem_trn_ is 1 for the allele-combinations in the bicluster ;
% (defined using tmp_ij_nlpR above). ;
%%%%%%%%;
pca_mr_A_trn_ = { mr_A_rem_trn_ };
pca_mr_Z_trn_ = { mr_Z_rem_trn_ };
pca_mc_A_trn = mc_A_rem_trn_; pca_mc_A_trn(1+index_Up99_setminus_cap_) = 0;
%%%%%%%%;
% Note that pca_mc_A_trn is 1 only for the allele-combinations in the bicluster ;
% that are also in the allele_cap. ;
%%%%%%%%;
pca_str_infix_trn=sprintf('trnUp99_tstAp99_ni%d_p01',tmp_ij_nlpR-1);
tmp_parameter_trn.flag_force_create = 0;
parameter_D_trnUp99_tstAp99_nix_p01 = tmp_parameter_trn;
parameter_D_trnUp99_tstAp99_nix_p01.str_A_p = str_Up99_A_p_p01;
parameter_D_trnUp99_tstAp99_nix_p01.str_name_s0000 = 'pca_D_trnUp99_tstAp99_nix_p01';
[ ...
 parameter_D_trnUp99_tstAp99_nix_p01 ...
,AZnV_D_trnUp99_tstAp99_nix_p01_ ...
,AnV_D_trnUp99_tstAp99_nix_p01_ ...
,ZnV_D_trnUp99_tstAp99_nix_p01_ ...
,V_D_trnUp99_tstAp99_nix_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_D_trnUp99_tstAp99_nix_p01 ...
,pca_rank ...
,pca_mr_A_trn_ ...
,pca_mr_Z_trn_ ...
,pca_mc_A_trn ...
,pca_str_infix_trn ...
);
%%%%;
mr_dvx_trnUp99_tstAp99_nix_ = 0.0*mr_A_ori_trn_;
mr_dvx_trnUp99_tstAp99_nix_(1+efind(mr_Z_alt_trn_)) = 1; %<-- color ctrl-patients. ;
mr_dvx_trnUp99_tstAp99_nix_(1+efind(mr_A_alt_trn_)) = 2; %<-- color case-patients. ;
mr_dvx_trnUp99_tstAp99_nix_(1+efind(mr_A_rem_trn_)) = 3; %<-- color bicluster-patients (a subset of cases). ;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_sml = 12;
markersize_big = 16;
hold on;
tmp_index_ = efind( (mr_dvx_trnUp99_tstAp99_nix_==1) | (mr_dvx_trnUp99_tstAp99_nix_==2) );
scatter(AZnV_D_trnUp99_tstAp99_nix_p01_(1+tmp_index_,1+0),AZnV_D_trnUp99_tstAp99_nix_p01_(1+tmp_index_,1+1),markersize_sml,mr_dvx_trnUp99_tstAp99_nix_(1+tmp_index_),'filled','MarkerEdgeColor','k');
tmp_index_ = efind( (mr_dvx_trnUp99_tstAp99_nix_==3) );
scatter(AZnV_D_trnUp99_tstAp99_nix_p01_(1+tmp_index_,1+0),AZnV_D_trnUp99_tstAp99_nix_p01_(1+tmp_index_,1+1),markersize_big,mr_dvx_trnUp99_tstAp99_nix_(1+tmp_index_),'filled','MarkerEdgeColor',0.85*[1,1,1]);
axis equal; grid on;
title('AZnV_D_trnUp99_tstAp99_nix_p01_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;
%%%%;

%%%%%%%%;
% save down first bicluster. ;
%%%%%%%%;
tmp_fam_iid_ = dataset_{1+ndataset_Up99}.fam_iid_(1+index_r_rem_trn_);
tmp_fname = sprintf('%s/%s_ni%d_fam_iid_.txt',parameter_trn.dir_out_s0000,parameter_trn.str_prefix,tmp_ij_nlpR-1);
if ~exist(tmp_fname,'file'); fp = fopen(tmp_fname,'w'); for nl=0:numel(tmp_fam_iid_)-1; fprintf(fp,'%s\n',tmp_fam_iid_{1+nl}); end; fclose(fp); end;
tmp_bim_vid_ = dataset_{1+ndataset_Up99}.bim_vid_(1+index_c_rem_trn_);
tmp_fname = sprintf('%s/%s_ni%d_bim_vid_.txt',parameter_trn.dir_out_s0000,parameter_trn.str_prefix,tmp_ij_nlpR-1);
if ~exist(tmp_fname,'file'); fp = fopen(tmp_fname,'w'); for nl=0:numel(tmp_bim_vid_)-1; fprintf(fp,'%s\n',tmp_bim_vid_{1+nl}); end; fclose(fp); end;
tmp_bim_name_ = dataset_{1+ndataset_Up99}.bim_name_(1+index_c_rem_trn_);
tmp_fname = sprintf('%s/%s_ni%d_bim_name_.txt',parameter_trn.dir_out_s0000,parameter_trn.str_prefix,tmp_ij_nlpR-1);
if ~exist(tmp_fname,'file'); fp = fopen(tmp_fname,'w'); for nl=0:numel(tmp_bim_name_)-1; fprintf(fp,'%s\n',tmp_bim_name_{1+nl}); end; fclose(fp); end;

%%%%%%%%;
% Now set up search for secondary bicluster. ;
% We will do so by 'scrambling' the first bicluster. ;
% The 'r0' suffix will correspond to scrambling the first bicluster. ;
% Note that the secondary bicluster is not very significant. ;
%%%%%%%%;
tmp_scramble_out_rdrop_trn_ = -1*ones(r_rem_trn,2);
tmp_scramble_out_rdrop_trn_(:,1+0) = index_r_rem_trn_;
tmp_scramble_out_cdrop_trn_ = -1*ones(c_rem_trn,2);
tmp_scramble_out_cdrop_trn_(:,1+1) = index_c_rem_trn_;
tmp_scramble_out_xdrop_trn_ = [ tmp_scramble_out_rdrop_trn_ ; tmp_scramble_out_cdrop_trn_ ];
str_scramble_out_xdrop = sprintf('%s/%s_ni%d_out_xdrop_a.txt',parameter_trn.dir_out_s0000,parameter_trn.str_prefix,tmp_ij_nlpR-1);
if ~exist(str_scramble_out_xdrop,'file');
fp=fopen(str_scramble_out_xdrop,'w');
fprintf(fp,'%d %d\n',transpose(tmp_scramble_out_xdrop_trn_));
fclose(fp);
end;%if ~exist(str_scramble_out_xdrop,'file');
%%%%%%%%;
parameter_trn_r0 = parameter_trn;
parameter_trn_r0.n_scramble = 1;
parameter_trn_r0.scramble_out_xdrop_ = {str_scramble_out_xdrop};
parameter_trn_r0.scramble_rseed_ = [1024*768*1];
parameter_trn_r0.nshuffle = 0;
parameter_trn_r0 = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_trn_r0);
%%%%%%%%;
trace_trn_r0____ = load_trace__from_dir_ver0(parameter_trn_r0.dir_out_trace);
trace_trn_r0____.ZR_s0000_ = (trace_trn_r0____.QR_s0000_ - trace_trn____.QR_avg_i_)./max(1e-12,trace_trn____.QR_std_i_);
%%%%%%%%;
fname_fig_pre = sprintf('%s/%s_trace_r0_FIGA',parameter_trn.dir_out_s0000_jpg,parameter_trn.str_name_s0000);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre); fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
figure(1+nf);nf=nf+1;clf;figsml;
linewidth_sml = 0.5;
linewidth_big = 3;
markersize_big = 16;
hold on;
plot(trace_trn____.niter_s0000_,trace_trn____.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(trace_trn____.niter_s0000_,trace_trn____.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(trace_trn____.niter_s0000_(tmp_ij_nlpR),trace_trn____.ZR_s0000_(tmp_ij_nlpR),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
plot(trace_trn____.niter_s0000_,trace_trn_r0____.ZR_s0000_,'g-','LineWidth',linewidth_big);
hold off;
xlim([min(trace_trn____.niter_s0000_),max(trace_trn____.niter_s0000_)]); xlabel('iteration');
ylabel('negative-log-p');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(parameter_trn.str_name_s0000,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;

%%%%%%%%;
% Now search again for secondary bicluster, ;
% this time explicitly excluding bc0 cases and snps. ;
% Note that the secondary bicluster is not very significant. ;
%%%%%%%%;
dir_trunk_Up99 = sprintf('%s/dir_Up99',dir_trunk);
ncontinent = 0;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% maf_lo_threshold %0.2f; ncontinent %d/%d;',maf_lo_threshold,ncontinent,n_continent)); end;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% setting parameters. ;')); end;
str_lak_vs_dex = 'dex';
str_prefix = 'test2mds_maf01';
gamma = 0.05;
n_mds_0in = 2; n_mds_repl = 1; ij_mds_use_ = [1:2];
parameter_0b.dir_code = dir_code;
parameter_0b.dir_trunk = dir_trunk_Up99;
parameter_0b.str_lak_vs_dex = str_lak_vs_dex;
parameter_0b.str_prefix = str_prefix;
parameter_0b.gamma = gamma;
parameter_0b.n_mds = n_mds_0in;
parameter_0b.n_mds_repl = n_mds_repl;
parameter_0b.ij_mds_use_ = ij_mds_use_;
parameter_0b.flag_force_create = 0; %<-- reload previous run. ;
parameter_0b.flag_verbose = max(0,flag_verbose-1);
parameter_0b.maf_lo_threshold = maf_lo_threshold;
parameter_0b.n_shuffle = 128;
parameter_0b.dir_0in = sprintf('%s/dir_%s',parameter_0b.dir_trunk,parameter_0b.str_prefix);
%%%%;
fname_mr_A_full = sprintf('%s/%s_mr_A_full.b16',parameter_0b.dir_0in,parameter_0b.str_prefix);
fname_mr_Z_full = sprintf('%s/%s_mr_Z_full.b16',parameter_0b.dir_0in,parameter_0b.str_prefix);
fname_mc_A = sprintf('%s/%s_mc_A.b16',parameter_0b.dir_0in,parameter_0b.str_prefix);
mr_A_ori_0b_ = binary_uncompress(fname_mr_A_full)>0; %<-- thresholded to avoid negative values. ;
mr_Z_ori_0b_ = binary_uncompress(fname_mr_Z_full)>0; %<-- thresholded to avoid negative values. ;
mc_A_ori_0b_ = binary_uncompress(fname_mc_A)>0; %<-- thresholded to avoid negative values. ;
assert(numel(mr_A_ori_0b_)==numel(mr_Up99_p01_continent_));
mr_A_alt_0b_ = mr_A_ori_0b_.*mr_Up99_p01_continent_pc__(:,1+ncontinent);
mr_A_alt_0b_(1+index_r_rem_trn_) = 0; %<-- Here we explicitly exclude the case-rows in the first bicluster. ;
mr_Z_alt_0b_ = mr_Z_ori_0b_.*mr_Up99_p01_continent_pc__(:,1+ncontinent);
mr_Z_alt_0b_(1+index_r_rem_trn_) = 0; %<-- Here we explicitly exclude the ctrl-rows in the first bicluster. ;
mc_A_alt_0b_ = mc_A_ori_0b_; 
mc_A_alt_0b_(1+index_c_rem_trn_) = 0; %<-- Here we explicitly exclude the cols in the first bicluster. ;
parameter_0b.str_mr_0in = sprintf('continent%d_mrb',ncontinent);
parameter_0b.str_mc_0in = sprintf('continent%d_mcb',ncontinent);
%%%%;
fname_mr_A_alt_0b_full = sprintf('%s/%s_mr_A_%s_full.b16',parameter_0b.dir_0in,parameter_0b.str_prefix,parameter_0b.str_mr_0in);
fname_mr_Z_alt_0b_full = sprintf('%s/%s_mr_Z_%s_full.b16',parameter_0b.dir_0in,parameter_0b.str_prefix,parameter_0b.str_mr_0in);
fname_mc_A_alt_0b = sprintf('%s/%s_mc_A_%s.b16',parameter_0b.dir_0in,parameter_0b.str_prefix,parameter_0b.str_mc_0in);
bitj=16;
binary_compress(bitj,mr_A_alt_0b_>0,fname_mr_A_alt_0b_full);
binary_compress(bitj,mr_Z_alt_0b_>0,fname_mr_Z_alt_0b_full);
binary_compress(bitj,mc_A_alt_0b_>0,fname_mc_A_alt_0b);
fname_mr_A_alt_0b_01 = sprintf('%s/%s_mr_A_%s_01.b16',parameter_0b.dir_0in,parameter_0b.str_prefix,parameter_0b.str_mr_0in);
fname_mr_Z_alt_0b_01 = sprintf('%s/%s_mr_Z_%s_01.b16',parameter_0b.dir_0in,parameter_0b.str_prefix,parameter_0b.str_mr_0in);
bitj=16;
binary_compress(bitj,mr_A_alt_0b_>0,fname_mr_A_alt_0b_01);
binary_compress(bitj,mr_Z_alt_0b_>0,fname_mr_Z_alt_0b_01);
%%%%;
for nshuffle=0:parameter_0b.n_shuffle-1+1;
parameter_0b.nshuffle = nshuffle;
parameter_0b = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_0b);
end;%for nshuffle=0:parameter_0b.n_shuffle-1+1;
%%%%%%%%;
parameter_0b.str_out_xdrop_a_s0000 = sprintf('%s/out_xdrop_a.txt',parameter_0b.dir_out_s0000);
xdrop_0b_ = load_out_xdrop_from_str_ver0(parameter_0b.str_out_xdrop_a_s0000);
index_rdrop_0b_ = xdrop_0b_.index_rdrop_;
index_cdrop_0b_ = xdrop_0b_.index_cdrop_;
parameter_0b.dir_out_s0000_jpg = sprintf('%s/dir_out_jpg',parameter_0b.dir_out_s0000);
if ~exist(parameter_0b.dir_out_s0000_jpg,'dir'); disp(sprintf(' %% mkdir %s',parameter_0b.dir_out_s0000_jpg)); mkdir(parameter_0b.dir_out_s0000_jpg); end;
trace_0b____ = load_trace__from_dir_ver0(parameter_0b.dir_out_trace);
[tmp_nlpR_0b,tmp_ij_nlpR_0b] = max(trace_0b____.nlpR_s0000_);
r_rem_0b = trace_0b____.r_rem_s0000_(tmp_ij_nlpR_0b);
c_rem_0b = trace_0b____.c_rem_s0000_(tmp_ij_nlpR_0b);
index_r_rem_0b_ = xdrop_0b_.index_rkeep_(1:r_rem_0b);
index_c_rem_0b_ = xdrop_0b_.index_ckeep_(1:c_rem_0b);
assert(numel(index_r_rem_0b_)==numel(intersect(efind(mr_A_alt_0b_),index_r_rem_0b_)));
assert(numel(index_c_rem_0b_)==numel(intersect(efind(mc_A_ori_0b_),index_c_rem_0b_)));
%%%%%%%%;
fname_fig_pre = sprintf('%s/%s_trace_FIGA',parameter_0b.dir_out_s0000_jpg,parameter_0b.str_name_s0000);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre); fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
figure(1+nf);nf=nf+1;clf;figsml;
linewidth_sml = 0.5;
linewidth_big = 3;
markersize_big = 16;
hold on;
plot(trace_0b____.niter_s0000_,trace_0b____.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(trace_0b____.niter_s0000_,trace_0b____.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(trace_0b____.niter_s0000_(tmp_ij_nlpR_0b),trace_0b____.ZR_s0000_(tmp_ij_nlpR_0b),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
hold off;
xlim([min(trace_0b____.niter_s0000_),max(trace_0b____.niter_s0000_)]); xlabel('iteration');
ylabel('negative-log-p');
sgtitle(parameter_0b.str_name_s0000,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;

flag_calc = 1;
if flag_calc;
%%%%%%%%;
% Now we calculate many _D_ principal-components across the entire Up99 data-set. ;
% The '_D_' refers to the situation where we calculate standard pcs using the cases ;
% (i.e., maximizing the variance across all the case-patients). ;
% First we set up the basic parameters. ;
%%%%%%%%;
str_lak_vs_dex = 'dex';
str_prefix = 'test2mds_maf01';
gamma = 0.05;
n_mds_0in = 2; n_mds_repl = 1; ij_mds_use_ = [1:2];
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
parameter_Up99.n_shuffle = 128;
parameter_Up99.dir_0in = sprintf('%s/dir_%s',parameter_Up99.dir_trunk,parameter_Up99.str_prefix);
parameter_Up99.nshuffle = 0;
parameter_Up99 = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_Up99);
%%%%%%%%;
fname_mr_A_full = sprintf('%s/%s_mr_A_full.b16',parameter_Up99.dir_0in,parameter_Up99.str_prefix);
fname_mr_Z_full = sprintf('%s/%s_mr_Z_full.b16',parameter_Up99.dir_0in,parameter_Up99.str_prefix);
fname_mc_A = sprintf('%s/%s_mc_A.b16',parameter_Up99.dir_0in,parameter_Up99.str_prefix);
mr_A_ori_Up99_ = binary_uncompress(fname_mr_A_full)>0; %<-- thresholded to avoid negative values. ;
mr_Z_ori_Up99_ = binary_uncompress(fname_mr_Z_full)>0; %<-- thresholded to avoid negative values. ;
mc_A_ori_Up99_ = binary_uncompress(fname_mc_A)>0; %<-- thresholded to avoid negative values. ;
%%%%%%%%;
end;%if flag_calc;

flag_calc = 1;
if flag_calc;
%%%%%%%%;
% And now we actually perform the pca calculation. ;
%%%%%%%%;
pca_rank = 2;
pca_mr_A_trn_ = { 1*mr_A_ori_trn_ };
pca_mr_Z_trn_ = { 1*mr_Z_ori_trn_ };
pca_mc_A_trn = mc_A_ori_trn_;
pca_str_infix_trn=sprintf('Up99_p01');
parameter_D_Up99_p01 = parameter_Up99;
parameter_D_Up99_p01.flag_force_create = 0;
parameter_D_Up99_p01.str_A_p = str_Up99_A_p_p01;
parameter_D_Up99_p01.str_name_s0000 = 'pca_D_Up99_p01';
[ ...
 parameter_D_Up99_p01 ...
,AZnV_D_Up99_p01_ ...
,AnV_D_Up99_p01_ ...
,ZnV_D_Up99_p01_ ...
,V_D_Up99_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_D_Up99_p01 ...
,pca_rank ...
,pca_mr_A_trn_ ...
,pca_mr_Z_trn_ ...
,pca_mc_A_trn ...
,pca_str_infix_trn ...
);
%%%%%%%%;
mx__ = load_mx__from_parameter_ver0(parameter_D_Up99_p01);
mr_dvx_ = 0.0*mx__.mr_A_full_;
mr_dvx_(1+efind(mx__.mr_A_full_))=2;
mr_dvx_(1+efind(mx__.mr_Z_full_))=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;fig80s;
%%%%;
subplot(1,2,1);
scatter(AZnV_D_Up99_p01_(:,1),AZnV_D_Up99_p01_(:,2),8,mr_dvx_,'filled');
xlabel('pc0'); ylabel('pc1'); title('case magenta, ctrl cyan');
axis equal; axis vis3d;
%%%%;
subplot(1,2,2);
scatter(AZnV_D_Up99_p01_(:,1),AZnV_D_Up99_p01_(:,2),8,mr_Up99_p01_continent_,'filled');
xlabel('pc0'); ylabel('pc1'); title('continent 0,1,2');
axis equal; axis vis3d;
%%%%;
fname_fig = sprintf('/%s/rangan/dir_bcc/dir_jelman/dir_Up99/Up99_D_p01_FIGA',str_home);
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if flag_disp;
%%%%%%%%;
end;%if flag_calc;

flag_calc = 1;
%%%%%%%%;
% Now calculate dominant _DandX_ principal-components across the entire Up99 data-set. ;
% The '_DandX_' refers to the situation where we calculate standard pcs using all the patients ;
% (i.e., maximizing the variance across all the case- and ctrl-patients provided). ;
% Note that we modify the parameters used above. ;
%%%%%%%%;
if flag_calc;
pca_rank = 2;
pca_mr_A_Up99_ = { 1*mr_A_ori_Up99_ + 1*mr_Z_ori_Up99_ };
pca_mr_Z_Up99_ = { 0*mr_A_ori_Up99_ + 0*mr_Z_ori_Up99_ };
pca_mc_A_Up99 = mc_A_ori_Up99_;
pca_str_infix_Up99=sprintf('Up99_DandX_p01');
parameter_DandX_Up99_p01 = parameter_Up99;
parameter_DandX_Up99_p01.flag_force_create = 0;
parameter_DandX_Up99_p01.str_A_p = str_Up99_A_p_p01;
parameter_DandX_Up99_p01.str_name_s0000 = 'pca_DandX_Up99_p01';
parameter_DandX_Up99_p01.slurm_memdecl = 48;
[ ...
 parameter_DandX_Up99_p01 ...
,AZnV_DandX_Up99_p01_ ...
,AnV_DandX_Up99_p01_ ...
,ZnV_DandX_Up99_p01_ ...
,V_DandX_Up99_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_DandX_Up99_p01 ...
,pca_rank ...
,pca_mr_A_Up99_ ...
,pca_mr_Z_Up99_ ...
,pca_mc_A_Up99 ...
,pca_str_infix_Up99 ...
);
%%%%%%%%;
mx__ = load_mx__from_parameter_ver0(parameter_DandX_Up99_p01);
mr_dvx_ = 0.0*mx__.mr_A_full_;
mr_dvx_(1+efind(mx__.mr_A_full_))=2;
mr_dvx_(1+efind(mx__.mr_Z_full_))=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;fig80s;
%%%%;
subplot(1,2,1);
%scatter3(AZnV_DandX_Up99_p01_(:,1),AZnV_DandX_Up99_p01_(:,2),AZnV_DandX_Up99_p01_(:,3),8,mr_dvx_,'filled');
scatter(AZnV_DandX_Up99_p01_(:,1),AZnV_DandX_Up99_p01_(:,2),8,mr_dvx_,'filled');
xlabel('pc0'); ylabel('pc1'); title('case magenta, ctrl cyan');
axis equal; axis vis3d;
%%%%;
subplot(1,2,2);
scatter(AZnV_DandX_Up99_p01_(:,1),AZnV_DandX_Up99_p01_(:,2),8,mr_Up99_p01_continent_,'filled');
xlabel('pc0'); ylabel('pc1'); title('continent 0,1,2');
axis equal; axis vis3d;
%%%%;
fname_fig = sprintf('/%s/rangan/dir_bcc/dir_jelman/dir_Up99/Up99_DandX_p01_FIGA',str_home);
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if flag_disp;
%%%%%%%%;
end;%if flag_calc;

flag_calc = 0;
if flag_calc;
%%%%%%%%;
% Now calculate many _DvX_ principal-components across the entire Up99 data-set. ;
% The '_DvX_' refers to the situation where we calculate 'category-corrected' pcs. ;
% These are designed to maximize the variance between cases and ctrls. ;
% Note that we modify the parameters used above. ;
%%%%%%%%;
pca_rank = 2;
pca_mr_A_trn_ = { 1*mr_A_ori_trn_ };
pca_mr_Z_trn_ = { 1*mr_Z_ori_trn_ };
pca_mc_A_trn = mc_A_ori_trn_;
pca_str_infix_trn=sprintf('Up99_p01');
parameter_DvX_Up99_p01 = parameter_Up99;
parameter_DvX_Up99_p01.flag_force_create = 0;
parameter_DvX_Up99_p01.str_A_p = str_Up99_A_p_p01;
parameter_DvX_Up99_p01.str_name_s0000 = 'pca_DvX_Up99_p01';
[ ...
 parameter_DvX_Up99_p01 ...
,AZnV_DvX_Up99_p01_ ...
,AnV_DvX_Up99_p01_ ...
,ZnV_DvX_Up99_p01_ ...
,V_DvX_Up99_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_DvX_from_mx_ver16( ...
 parameter_DvX_Up99_p01 ...
,pca_rank ...
,pca_mr_A_trn_ ...
,pca_mr_Z_trn_ ...
,pca_mc_A_trn ...
,pca_str_infix_trn ...
);
%%%%%%%%;
mx__ = load_mx__from_parameter_ver0(parameter_DvX_Up99_p01);
mr_dvx_ = 0.0*mx__.mr_A_full_;
mr_dvx_(1+efind(mx__.mr_A_full_))=2;
mr_dvx_(1+efind(mx__.mr_Z_full_))=1;
scatter(AZnV_DvX_Up99_p01_(:,1),AZnV_DvX_Up99_p01_(:,2),8,mr_dvx_);
%%%%%%%%;
end;%if flag_calc;

%%%%%%%%;
% Here we can take a quick peek at the histogram of cases (magenta) as ;
% projected onto the dominant principal-component of the bicluster. ;
% The ctrls are shown in blue. ;
% Note that this projection is limited to the training-set (Up99). ;
%%%%%%%%;
x_ = transpose(linspace(-60,+80,128)); 
xlim_ = [min(x_),max(x_)];
tmp_h_Z_ = hist(AZnV_D_trnUp99_tstAp99_nix_p01_(1+efind(mr_Z_alt_trn_),1+0),x_);
tmp_h_Z_ = tmp_h_Z_/sum(tmp_h_Z_);
tmp_h_A_ = hist(AZnV_D_trnUp99_tstAp99_nix_p01_(1+efind(mr_A_alt_trn_),1+0),x_);
tmp_h_A_ = tmp_h_A_/sum(tmp_h_A_);
figure(1+nf);nf=nf+1;figsml;
hold on;
stairs(x_,tmp_h_Z_,'b-');
stairs(x_,tmp_h_A_,'m-');
hold off;
xlim(xlim_);
%%%%%%%%;

%%%%%%%%;
% Now, for the training dataset Up99, ;
% We define the 'projected value' of each patient as the ;
% position of that patient after projection onto the ;
% dominant pc of the first bicluster ;
% (defined using only the allele-combinations in the bicluster which lie within allele_cap). ;
% Then we calculate the auc between cases and ctrls ;
% for those patients above a certain threshold of this projected value. ;
% We store this in: ;
% D_trnUp99_tstAp99_nix_p01_auc_t_, auc values. ;
% D_trnUp99_tstAp99_nix_p01_logp_auc_t_, log p-values for the auc. ;
% D_trnUp99_tstAp99_nix_p01_thr_t_, threshold values. ;
%%%%%%%%
tmp_Z_ = AZnV_D_trnUp99_tstAp99_nix_p01_(1+efind(mr_Z_alt_trn_),1+0);
tmp_A_ = AZnV_D_trnUp99_tstAp99_nix_p01_(1+efind(mr_A_alt_trn_),1+0);
tmp_threshold_ = sort([tmp_Z_;tmp_A_],'ascend'); tmp_n_threshold = numel(tmp_threshold_);
D_trnUp99_tstAp99_nix_p01_auc_t_ = zeros(tmp_n_threshold,1);
D_trnUp99_tstAp99_nix_p01_logp_auc_t_ = zeros(tmp_n_threshold,1);
for tmp_nthreshold=0:tmp_n_threshold-1;
tmp_threshold = tmp_threshold_(1+tmp_nthreshold);
tmp_Z_thr_ = tmp_Z_(1+efind(tmp_Z_>=tmp_threshold));
tmp_A_thr_ = tmp_A_(1+efind(tmp_A_>=tmp_threshold));
tmp_auc = auc_0(tmp_Z_thr_,tmp_A_thr_);
tmp_n = min(numel(tmp_Z_thr_),numel(tmp_A_thr_));
tmp_logp_auc = logp_auc_0(tmp_auc,tmp_n);
if (flag_verbose>1); disp(sprintf(' %% thr %0.2f: tmp_auc: %0.3f <-- nlp %0.2f <-- p %0.6f',tmp_threshold,tmp_auc,-tmp_logp_auc,exp(tmp_logp_auc))); end;
D_trnUp99_tstAp99_nix_p01_auc_t_(1+tmp_nthreshold) = tmp_auc;
D_trnUp99_tstAp99_nix_p01_logp_auc_t_(1+tmp_nthreshold) = tmp_logp_auc;
end;%for tmp_nthreshold=0:tmp_n_threshold-1;
D_trnUp99_tstAp99_nix_p01_thr_t_ = tmp_threshold_;

%%%%%%%%;
% Now project Ap99 data onto the principal-components defined using the bicluster. ;
% Note that we have to use those principal-components defined using the intersection ;
% of the bicluster with the allele_cap (see pca_str_infix_trn on the next line). ;
%%%%%%%%;
pca_str_infix_trn=sprintf('trnUp99_tstAp99_ni%d_p01',tmp_ij_nlpR-1); %<-- this defines the pca file to use. ;
parameter_D_trnUp99_tstAp99_nix_p01.str_V_ = sprintf('%s/dir_pca/dir_pca_mda/pca_D_%s_k%d_B44_V_.mda',parameter_D_trnUp99_tstAp99_nix_p01.dir_out_s0000,pca_str_infix_trn,pca_rank);
if ~exist(parameter_D_trnUp99_tstAp99_nix_p01.str_V_,'file');
disp(sprintf(' %% Warning, %s not found',parameter_D_trnUp99_tstAp99_nix_p01.str_V_));
end;%if ~exist(parameter_D_trnUp99_tstAp99_nix_p01.str_V_,'file');
D_trnUp99_tstAp99_nix_p01_V_ = mda_read_r8(parameter_D_trnUp99_tstAp99_nix_p01.str_V_);
if (flag_verbose); disp(sprintf(' %% size of intersection: %d',numel(index_Up99_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% Support of V_: %d',sum(abs(D_trnUp99_tstAp99_nix_p01_V_(:,1))>0))); end;
if (flag_verbose); disp(sprintf(' %% intersection of supp(V_) with cap: %d',numel(intersect(index_Up99_from_cap_,efind(abs(D_trnUp99_tstAp99_nix_p01_V_(:,1))>0))))); end;
D_Ap99_from_trnUp99_tstAp99_nix_p01_V_ = ij_Ap99_from_Up99__*D_trnUp99_tstAp99_nix_p01_V_;
if (flag_verbose); disp(sprintf(' %% fnorm(D_Ap99_from_trnUp99_tstAp99_nix_p01_V_(ij_Ap99_from_cap_,:) - D_trnUp99_tstAp99_nix_p01_V_(ij_Up99_from_cap_,:)): %0.16f',fnorm(D_Ap99_from_trnUp99_tstAp99_nix_p01_V_(ij_Ap99_from_cap_,:) - D_trnUp99_tstAp99_nix_p01_V_(ij_Up99_from_cap_,:)))); end;

%%%%%%%%;
% Now project Ap99 data onto D_Ap99_from_trnUp99_tstAp99_nix_p01_V_. ;
%%%%%%%%;
ndataset=ndataset_Ap99;
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
assert(numel(mr_A_ori_tst_)==numel(mr_Ap99_p01_continent_));
mr_A_alt_tst_ = mr_A_ori_tst_.*mr_Ap99_p01_continent_pc__(:,1+ncontinent);
mr_Z_alt_tst_ = mr_Z_ori_tst_.*mr_Ap99_p01_continent_pc__(:,1+ncontinent);
mc_A_alt_tst_ = mc_A_ori_tst_;
pca_mr_A_tst_ = { mr_A_alt_tst_ };
pca_mr_Z_tst_ = { mr_Z_alt_tst_ };
pca_mc_A_tst = mc_A_alt_tst_; pca_mc_A_tst(1+index_Ap99_setminus_cap_) = 0;
pca_str_infix_tst='D_Ap99_from_trnUp99_tstAp99_nix_p01';
parameter_D_Ap99_from_trnUp99_tstAp99_nix_p01 = tmp_parameter_tst;
parameter_D_Ap99_from_trnUp99_tstAp99_nix_p01.str_A_p = str_Ap99_A_p_p01;
parameter_D_Ap99_from_trnUp99_tstAp99_nix_p01.str_name_s0000 = ...
  sprintf('pca_D_Ap99c%d_from_trnUp99c%d_tstAp99c%d_nix_p01',ncontinent,ncontinent,ncontinent);
[ ...
 parameter_D_Ap99_from_trnUp99_tstAp99_nix_p01 ...
,AZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_ ...
,AnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_ ...
,ZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_ ...
,tmp_V_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_proj_from_mx_ver16( ...
 parameter_D_Ap99_from_trnUp99_tstAp99_nix_p01 ...
,pca_rank ...
,D_Ap99_from_trnUp99_tstAp99_nix_p01_V_ ...
,pca_mr_A_tst_ ...
,pca_mr_Z_tst_ ...
,pca_mc_A_tst ...
,pca_str_infix_tst ...
);
%%%%;
mr_dvx_Ap99_from_trnUp99_tstAp99_nix_ = 0.0*mr_A_ori_tst_;
mr_dvx_Ap99_from_trnUp99_tstAp99_nix_(1+efind(mr_Z_alt_tst_)) = 1;
mr_dvx_Ap99_from_trnUp99_tstAp99_nix_(1+efind(mr_A_alt_tst_)) = 2;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_use = 12;
tmp_index_ = efind(mr_dvx_Ap99_from_trnUp99_tstAp99_nix_);
scatter(AZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_(1+tmp_index_,1+0),AZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_(1+tmp_index_,1+1),markersize_use,mr_dvx_Ap99_from_trnUp99_tstAp99_nix_(1+tmp_index_),'filled','MarkerEdgeColor','k');
axis equal; grid on;
title('AZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
fname_fig_pre = sprintf('%s/AZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%;

%%%%%%%%;
% Here we can take a quick peek at the histogram of cases (magenta) as ;
% projected onto the dominant principal-component of the bicluster. ;
% The ctrls are shown in blue. ;
% Note that this projection is limited to the testing-set (Ap99). ;
%%%%%%%%;
x_ = transpose(linspace(-60,+80,128)); 
xlim_ = [min(x_),max(x_)];
tmp_h_Z_ = hist(AZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_(1+efind(mr_Z_alt_tst_),1+0),x_);
tmp_h_Z_ = tmp_h_Z_/sum(tmp_h_Z_);
tmp_h_A_ = hist(AZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_(1+efind(mr_A_alt_tst_),1+0),x_);
tmp_h_A_ = tmp_h_A_/sum(tmp_h_A_);
figure(1+nf);nf=nf+1;figsml;
hold on;
stairs(x_,tmp_h_Z_,'b-');
stairs(x_,tmp_h_A_,'m-');
hold off;
xlim(xlim_);

%%%%%%%%;
% Now, for the testing dataset Ap99, ;
% We define the 'projected value' of each patient as the ;
% position of that patient after projection onto the ;
% dominant pc of the first bicluster ;
% (defined using only the allele-combinations in the bicluster which lie within allele_cap). ;
% Then we calculate the auc between cases and ctrls ;
% for those patients above a certain threshold of this projected value. ;
% We store this in: ;
% D_Ap99_from_trnUp99_tstAp99_nix_p01_auc_t_, auc values. ;
% D_Ap99_from_trnUp99_tstAp99_nix_p01_logp_auc_t_, log p-values for the auc. ;
% D_Ap99_from_trnUp99_tstAp99_nix_p01_thr_t_, threshold values. ;
%%%%%%%%
tmp_Z_ = AZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_(1+efind(mr_Z_alt_tst_),1+0);
tmp_A_ = AZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_(1+efind(mr_A_alt_tst_),1+0);
tmp_threshold_ = sort([tmp_Z_;tmp_A_],'ascend'); tmp_n_threshold = numel(tmp_threshold_);
D_Ap99_from_trnUp99_tstAp99_nix_p01_auc_t_ = zeros(tmp_n_threshold,1);
D_Ap99_from_trnUp99_tstAp99_nix_p01_logp_auc_t_ = zeros(tmp_n_threshold,1);
for tmp_nthreshold=0:tmp_n_threshold-1;
tmp_threshold = tmp_threshold_(1+tmp_nthreshold);
tmp_Z_thr_ = tmp_Z_(1+efind(tmp_Z_>=tmp_threshold));
tmp_A_thr_ = tmp_A_(1+efind(tmp_A_>=tmp_threshold));
tmp_auc = auc_0(tmp_Z_thr_,tmp_A_thr_);
tmp_n = min(numel(tmp_Z_thr_),numel(tmp_A_thr_));
tmp_logp_auc = logp_auc_0(tmp_auc,tmp_n);
if (flag_verbose>1); disp(sprintf(' %% thr %0.2f: tmp_auc: %0.3f <-- nlp %0.2f <-- p %0.6f',tmp_threshold,tmp_auc,-tmp_logp_auc,exp(tmp_logp_auc))); end;
D_Ap99_from_trnUp99_tstAp99_nix_p01_auc_t_(1+tmp_nthreshold) = tmp_auc;
D_Ap99_from_trnUp99_tstAp99_nix_p01_logp_auc_t_(1+tmp_nthreshold) = tmp_logp_auc;
end;%for tmp_nthreshold=0:tmp_n_threshold-1;
D_Ap99_from_trnUp99_tstAp99_nix_p01_thr_t_ = tmp_threshold_;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Now we plot the log-p-values for the auc of cases-vs-ctrls ;
% for each of the thresholds used above.
% This is after we project onto the dominant pc of the bicluster, ;
% defined using only the shared allele-combinations. ;
% For this plot the training-data (Up99) is shown in red, ;
% while the testing-data (Ap99) is shown in green. ;
% Note that there is a large range of thresholds for which ;
% the testing-data is very significant. ;
% (significance value of 0.05 is marked with dashed line). ;
% (everything above dashed line is significant). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
ymax = 10;
tmp_dlim_0_ = [min(D_trnUp99_tstAp99_nix_p01_thr_t_);max(D_trnUp99_tstAp99_nix_p01_thr_t_)];
tmp_dlim_1_ = [min(D_Ap99_from_trnUp99_tstAp99_nix_p01_thr_t_);max(D_Ap99_from_trnUp99_tstAp99_nix_p01_thr_t_)];
tmp_dlim_ = [min([tmp_dlim_0_;tmp_dlim_1_]),max([tmp_dlim_0_;tmp_dlim_1_])];
stairs(D_trnUp99_tstAp99_nix_p01_thr_t_,min(10,-D_trnUp99_tstAp99_nix_p01_logp_auc_t_),'r-','LineWidth',2);
stairs(D_Ap99_from_trnUp99_tstAp99_nix_p01_thr_t_,min(10,-D_Ap99_from_trnUp99_tstAp99_nix_p01_logp_auc_t_),'g-','LineWidth',2);
%stairs(D_Ap99_from_trnUp99_tstAp99_nix_p01_thr_t_,min(10,10*D_Ap99_from_trnUp99_tstAp99_nix_p01_auc_t_),'-','Color',0.85*[1,1,1],'LineWidth',2);
plot(tmp_dlim_,-log(0.05)*ones(2,1),'k:');
hold off;
xlim(tmp_dlim_);
ylim([0,10.1]); grid on;
xlabel('pc1 threshold');
ylabel('-log(p(AUC))','Interpreter','none');
legend({'train','test'},'Location','NorthWest');
title('AZnV_D_trnUp99_tstAp99_nix_p01_','Interpreter','none');
fname_fig_pre = sprintf('%s/test_Up99_vs_Ap99_18b_FIGA',dir_jpg); %<-- note that this matches the previous test_Up99_vs_Ap99_18b.m ;
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
end;%if flag_disp;

dir_mat_replication = sprintf('%s/dir_mat_replication',dir_trunk);
if ~exist(dir_mat_replication,'dir'); disp(sprintf(' %% mkdir %s',dir_mat_replication')); mkdir(dir_mat_replication); end;
fname_mat = sprintf('%s/trnUp99_tst_Ap99_ncontinent%d.mat',dir_mat_replication,ncontinent);
if ~exist(fname_mat,'file');
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_mat ...
     ,'AZnV_D_Ap99_from_trnUp99_tstAp99_nix_p01_' ...
     ,'AZnV_D_trnUp99_tstAp99_nix_p01_' ...
     ,'mr_A_alt_trn_' ...
     ,'mr_A_alt_tst_' ...
     ,'mr_A_rem_trn_' ...
     ,'mr_Z_alt_trn_' ...
     ,'mr_Z_alt_tst_' ...
     ,'mr_dvx_Ap99_from_trnUp99_tstAp99_nix_' ...
     ,'mr_dvx_trnUp99_tstAp99_nix_' ...
     );
end;%if ~exist(fname_mat,'file');
if  exist(fname_mat,'file');
disp(sprintf(' %% %s found, not creating',fname_mat));
end;%if  exist(fname_mat,'file');

