function test_Up05_vs_Ap05_23(ncontinent_use);

%%%%%%%%;
% which continent to use (0based). ;
%%%%%%%%;
if (nargin<1); ncontinent_use=[]; end;
if isempty(ncontinent_use); ncontinent_use = 0; end;

platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); str_home = 'data'; end;
if (strcmp(platform,'OptiPlex')); str_home = 'home'; end;
if (strcmp(platform,'eval1')); str_home = 'home'; end;
if (strcmp(platform,'rusty')); str_home = 'mnt/home'; end;
%%%%%%%%;
run(sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev/dir_m/setup_0',str_home)); %<-- set up the paths. ;
addpath(sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev/dir_m_dependencies/',str_home));
addpath(sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev/dir_m/',str_home));
addpath(sprintf('/%s/rangan/dir_bcc/dir_jelman/dir_m/',str_home));
flag_verbose = 1;
flag_disp = 1+flag_verbose; nf=0;
flag_replot = 0;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Assume path is set using dir_lakcluster_c/dir_m/setup_0.m. ;')); end;

dir_code = sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev',str_home);
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman/dir_AD_Biclustering_20230315/dir_ADNI_vs_UKB',str_home);
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

%%%%%%%%;
% first we load the 2 datasets. ;
%%%%%%%%;
n_dataset = 2;
str_dataset_ = {'Up05','Ap05'};
dataset_ = cell(n_dataset,1);
for ndataset=0:n_dataset-1;
dataset = struct('type','dataset');
dataset.str_dataset = str_dataset_{1+ndataset};
dataset.dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman/dir_AD_Biclustering_20230315/dir_ADNI_vs_UKB/dir_%s',str_home,dataset.str_dataset);
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
% index_Up05_from_cap_ cross-references the allele_cap and the Up05 dataset. ;
% index_Ap05_from_cap_ cross-references the allele_cap and the Ap05 dataset. ;
% index_Up05_setminus_cap_ lists those allele-combinations in Up05 which are not part of the intersection. ;
% index_Ap05_setminus_cap_ lists those allele-combinations in Ap05 which are not part of the intersection. ;
%%%%%%%%;
ndataset_Up05 = 0;
ndataset_Ap05 = 1;
n_patient_Up05 = dataset_{1+ndataset_Up05}.n_patient;
n_patient_Ap05 = dataset_{1+ndataset_Ap05}.n_patient;
n_snp_Up05 = dataset_{1+ndataset_Up05}.n_snp;
n_snp_Ap05 = dataset_{1+ndataset_Ap05}.n_snp;
[allele_cap_,ij_Up05_from_cap_,ij_Ap05_from_cap_] = intersect(dataset_{1+ndataset_Up05}.bim_name_,dataset_{1+ndataset_Ap05}.bim_name_,'stable');
index_Up05_from_cap_ = ij_Up05_from_cap_ - 1;
index_Ap05_from_cap_ = ij_Ap05_from_cap_ - 1;
ij_Up05_from_Ap05__ = sparse(ij_Up05_from_cap_,ij_Ap05_from_cap_,1,n_snp_Up05,n_snp_Ap05);
ij_Ap05_from_Up05__ = sparse(ij_Ap05_from_cap_,ij_Up05_from_cap_,1,n_snp_Ap05,n_snp_Up05);
index_Up05_setminus_cap_ = setdiff(efind(dataset_{1+ndataset_Up05}.mx__.mc_A_),index_Up05_from_cap_);
index_Ap05_setminus_cap_ = setdiff(efind(dataset_{1+ndataset_Ap05}.mx__.mc_A_),index_Ap05_from_cap_);
if (flag_verbose); disp(sprintf(' %% numel(index_Up05_from_cap_): %d',numel(index_Up05_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% numel(index_Ap05_from_cap_): %d',numel(index_Ap05_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% numel(index_Up05_setminus_cap_): %d',numel(index_Up05_setminus_cap_))); end;
if (flag_verbose); disp(sprintf(' %% numel(index_Ap05_setminus_cap_): %d',numel(index_Ap05_setminus_cap_))); end;

%%%%%%%%;
% The number of pcs we will calculate is pca_rank_use. ;
% The minor-allele-frequency we will use is maf_lo_threshold. ;
%%%%%%%%;
pca_rank_use = 2;
maf_lo_threshold = 0.01;

%%%%%%%%;
% Over the course of the next two paragraphs we define Up05_A_p_ and Ap05_A_p_. ;
% These are the 'baseline' values expected for each allele-combination in each dataset. ;
%%%%%%%%;

%%%%%%%%;
% Set Up05_A_p_. ;
%%%%%%%%;
ndataset=ndataset_Up05;
mx_tst__ = dataset_{1+ndataset}.mx__;
tmp_parameter_tst = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.dir_0in: %s',tmp_parameter_tst.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.str_prefix: %s',tmp_parameter_tst.str_prefix)); end;
tmp_parameter_tst.dir_code = sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev',str_home);
tmp_parameter_tst.maf_lo_threshold = maf_lo_threshold;
tmp_parameter_tst.maf_hi_threshold = 0.50;
tmp_parameter_tst.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter_tst.flag_verbose = 0;
tmp_parameter_tst.flag_force_create = 0;
pca_rank = pca_rank_use;
mr_A_ori_tst_ = mx_tst__.mr_A_full_;
mr_Z_ori_tst_ = mx_tst__.mr_Z_full_;
mc_A_ori_tst_ = mx_tst__.mc_A_;
pca_str_infix_tst='D_trnUp05_tstAp05_nix_p01';
%%%%;
parameter_A_p_p01 = tmp_parameter_tst;
pca_A_p_str_infix = 'full';
parameter_A_p_p01.str_name_s0000 = 'A_p_p01';
parameter_A_p_p01.flag_force_create = 0;
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
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot(tmp_Up05_A_p_p01_);
xlabel('pcol');ylabel('A_p','Interpreter','none');
title(str_Up05_A_p_p01,'Interpreter','none');
end;%if (flag_disp);

%%%%%%%%;
% Set Ap05_A_p_. ;
%%%%%%%%;
ndataset=ndataset_Ap05;
mx_tst__ = dataset_{1+ndataset}.mx__;
tmp_parameter_tst = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.dir_0in: %s',tmp_parameter_tst.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.str_prefix: %s',tmp_parameter_tst.str_prefix)); end;
tmp_parameter_tst.dir_code = sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev',str_home);
tmp_parameter_tst.maf_lo_threshold = maf_lo_threshold;
tmp_parameter_tst.maf_hi_threshold = 0.50;
tmp_parameter_tst.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter_tst.flag_verbose = 0;
tmp_parameter_tst.flag_force_create = 0;
pca_rank = pca_rank_use;
mr_A_ori_tst_ = mx_tst__.mr_A_full_;
mr_Z_ori_tst_ = mx_tst__.mr_Z_full_;
mc_A_ori_tst_ = mx_tst__.mc_A_;
pca_str_infix_tst='D_Ap05_from_trnUp05_tstAp05_nix_p01';
%%%%;
parameter_A_p_p01 = tmp_parameter_tst;
pca_A_p_str_infix = 'full';
parameter_A_p_p01.str_name_s0000 = 'A_p_p01';
parameter_A_p_p01.flag_force_create = 0;
[ ...
 parameter_A_p_p01 ...
,str_Ap05_A_p_p01 ...
] = ...
xxxcluster_fromdisk_uADZSZDA_A_p_from_mx_ver16( ...
parameter_A_p_p01 ...
,{ mr_A_ori_tst_ } ...
,{ mr_Z_ori_tst_ } ...
,mc_A_ori_tst_ ...
,pca_A_p_str_infix ...
);
tmp_Ap05_A_p_p01_ = mda_read_r8(str_Ap05_A_p_p01); 
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot(tmp_Ap05_A_p_p01_);
xlabel('pcol');ylabel('A_p','Interpreter','none');
title(str_Ap05_A_p_p01,'Interpreter','none');
end;%if (flag_disp);

%%%%%%%%;
% Here we load the continent-labels for the datasets. ;
% Note that the new labels (from dir_AD_Biclustering_20230315/dir_ADNI_vs_UKB) ;
% differ in length and structure from the old labels ;
% (stored in mr_Up05_p01_continent_.txt and mr_Ap05_from_0UKB_cap_Ap05_p01_continent_.txt). ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% We divide Up05 data into three continents. ;')); end;
dir_mat_replication_old = sprintf('%s/dir_mat_replication',dir_trunk);
if ~exist(dir_mat_replication_old,'dir'); disp(sprintf(' %% mkdir %s',dir_mat_replication_old)); mkdir(dir_mat_replication_old); end;
ncontinent=0; ncontinent_1based = ncontinent + 1;
str_datafile = sprintf('trnUp05_tst_Ap05_ncontinent%d',ncontinent_1based); %<-- note that in this filename ncontinent is 1-based rather than 0-based. ;
mr_ncontinent0_ = load(sprintf('%s/%s.mat',dir_mat_replication_old,str_datafile));
ncontinent=1; ncontinent_1based = ncontinent + 1;
str_datafile = sprintf('trnUp05_tst_Ap05_ncontinent%d',ncontinent_1based); %<-- note that in this filename ncontinent is 1-based rather than 0-based. ;
mr_ncontinent1_ = load(sprintf('%s/%s.mat',dir_mat_replication_old,str_datafile));
clear ncontinent ncontinent_1based str_datafile ;
%%%%%%%%;
k_gamma_use_ncontinent0 = 0.50; k_gamma_use_ncontinent1 = 0.15; %<-- Jeremy: Continent1 converges at .50, continent2 must be lowered to 0.15 ;
%%%%%%%%;
% ensure that the continent-labels do not overlap. ;
%%%%%%%%;
assert(sum(mr_ncontinent0_.mr_dvx_trnUp05_tstAp05_nix_ & mr_ncontinent1_.mr_dvx_trnUp05_tstAp05_nix_)==0);
assert(numel(mr_ncontinent0_.mr_dvx_trnUp05_tstAp05_nix_)==n_patient_Up05);
assert(numel(mr_ncontinent1_.mr_dvx_trnUp05_tstAp05_nix_)==n_patient_Up05);
assert(sum(mr_ncontinent0_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_ & mr_ncontinent1_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_)==0);
assert(numel(mr_ncontinent0_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_)==n_patient_Ap05);
assert(numel(mr_ncontinent1_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_)==n_patient_Ap05);
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Up05 ncontinent 1+0=1: label==1: %0.4d',sum(mr_ncontinent0_.mr_dvx_trnUp05_tstAp05_nix_==1))); end;%if (flag_verbose>0); 
if (flag_verbose>0); disp(sprintf(' %% Up05 ncontinent 1+0=1: label==2: %0.4d',sum(mr_ncontinent0_.mr_dvx_trnUp05_tstAp05_nix_==2))); end;%if (flag_verbose>0); 
if (flag_verbose>0); disp(sprintf(' %% Up05 ncontinent 1+1=1: label==1: %0.4d',sum(mr_ncontinent1_.mr_dvx_trnUp05_tstAp05_nix_==1))); end;%if (flag_verbose>0); 
if (flag_verbose>0); disp(sprintf(' %% Up05 ncontinent 1+1=1: label==2: %0.4d',sum(mr_ncontinent1_.mr_dvx_trnUp05_tstAp05_nix_==2))); end;%if (flag_verbose>0); 
if (flag_verbose>0); disp(sprintf(' %% Ap05 ncontinent 1+0=1: label==1: %0.4d',sum(mr_ncontinent0_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1))); end;%if (flag_verbose>0); 
if (flag_verbose>0); disp(sprintf(' %% Ap05 ncontinent 1+0=1: label==2: %0.4d',sum(mr_ncontinent0_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2))); end;%if (flag_verbose>0); 
if (flag_verbose>0); disp(sprintf(' %% Ap05 ncontinent 1+1=1: label==1: %0.4d',sum(mr_ncontinent1_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1))); end;%if (flag_verbose>0); 
if (flag_verbose>0); disp(sprintf(' %% Ap05 ncontinent 1+1=1: label==2: %0.4d',sum(mr_ncontinent1_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2))); end;%if (flag_verbose>0); 
%%%%%%%%;
mr_ori_Up05_p01_continent_pc__ = [ (mr_ncontinent0_.mr_dvx_trnUp05_tstAp05_nix_>0) , (mr_ncontinent1_.mr_dvx_trnUp05_tstAp05_nix_>0) , zeros(n_patient_Up05,1) ];
mr_ori_Up05_p01_continent_pc__(:,1+2) = ~(mr_ori_Up05_p01_continent_pc__(:,1+0) | mr_ori_Up05_p01_continent_pc__(:,1+1));
mr_ori_Ap05_p01_continent_pc__ = [ (mr_ncontinent0_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_>0) , (mr_ncontinent1_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_>0) , zeros(n_patient_Ap05,1) ];
mr_ori_Ap05_p01_continent_pc__(:,1+2) = ~(mr_ori_Ap05_p01_continent_pc__(:,1+0) | mr_ori_Ap05_p01_continent_pc__(:,1+1));
%%%%%%%%;

%%%%%%%%;
% Now we specify the biclustering run to use. ;
% To begin with, we limit ourselves to the training dataset Up05, and a specific ncontinent. ;
% This continent will be used to define: ;
% mr_A_alt_trn_, the list of case-patients in that continent (within the training dataset Up05). ;
% mr_Z_alt_trn_, the list of ctrl-patients in that continent (within the training dataset Up05). ;
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
dir_trunk_Up05 = sprintf('%s/dir_Up05',dir_trunk);
ncontinent = ncontinent_use; %<-- here we set the continent to use below. ;
ncontinent_1based = ncontinent + 1; %<-- This is used to name files. ;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% maf_lo_threshold %0.2f; ncontinent %d;',maf_lo_threshold,ncontinent)); end;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% setting parameters. ;')); end;
str_lak_vs_dex = 'dex';
str_prefix = 'test2mds_maf01';
gamma = 0.05;
n_mds_0in = 2; n_mds_repl = 1; ij_mds_use_ = [1:2];
parameter_trn.dir_code = dir_code;
parameter_trn.dir_trunk = dir_trunk_Up05;
parameter_trn.str_lak_vs_dex = str_lak_vs_dex;
parameter_trn.str_prefix = str_prefix;
parameter_trn.gamma = gamma;
parameter_trn.n_mds = n_mds_0in;
parameter_trn.n_mds_repl = n_mds_repl;
parameter_trn.ij_mds_use_ = ij_mds_use_;
parameter_trn.flag_force_create = 0; %<-- reload previous run. ;
parameter_trn.flag_verbose = max(0,flag_verbose-1);
parameter_trn.maf_lo_threshold = maf_lo_threshold;
parameter_trn.n_shuffle = 128;
parameter_trn.dir_0in = sprintf('%s/dir_%s',parameter_trn.dir_trunk,parameter_trn.str_prefix);
fname_mr_A_full = sprintf('%s/%s_mr_A_full.b16',parameter_trn.dir_0in,parameter_trn.str_prefix);
fname_mr_Z_full = sprintf('%s/%s_mr_Z_full.b16',parameter_trn.dir_0in,parameter_trn.str_prefix);
fname_mc_A = sprintf('%s/%s_mc_A.b16',parameter_trn.dir_0in,parameter_trn.str_prefix);
mr_A_ori_trn_ = binary_uncompress(fname_mr_A_full)>0; %<-- thresholded to avoid negative values. ;
mr_Z_ori_trn_ = binary_uncompress(fname_mr_Z_full)>0; %<-- thresholded to avoid negative values. ;
mc_A_ori_trn_ = binary_uncompress(fname_mc_A)>0; %<-- thresholded to avoid negative values. ;
mr_A_alt_trn_ = mr_A_ori_trn_.*mr_ori_Up05_p01_continent_pc__(:,1+ncontinent);
mr_Z_alt_trn_ = mr_Z_ori_trn_.*mr_ori_Up05_p01_continent_pc__(:,1+ncontinent);
parameter_trn.str_mr_0in = sprintf('continent%d',ncontinent_1based); %<-- note that this is now 1-based. ;
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
figure(1+nf);nf=nf+1;clf;figmed;
linewidth_sml = 0.5;
linewidth_big = 3;
markersize_big = 16;
subplot(1,2,1);
hold on;
plot(trace_trn____.niter_s0000_,trace_trn____.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(trace_trn____.niter_s0000_,trace_trn____.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(trace_trn____.niter_s0000_(tmp_ij_nlpR),trace_trn____.ZR_s0000_(tmp_ij_nlpR),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
hold off;
xlim([min(trace_trn____.niter_s0000_),max(trace_trn____.niter_s0000_)]); xlabel('iteration');
ylabel('negative-log-p');
subplot(1,2,2);
hold on;
plot(max(trace_trn____.r_rem_s0000_)-trace_trn____.r_rem_s0000_,trace_trn____.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(max(trace_trn____.r_rem_s0000_)-trace_trn____.r_rem_s0000_,trace_trn____.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(max(trace_trn____.r_rem_s0000_)-trace_trn____.r_rem_s0000_(tmp_ij_nlpR),trace_trn____.ZR_s0000_(tmp_ij_nlpR),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
hold off;
xlim([min(trace_trn____.r_rem_s0000_),max(trace_trn____.r_rem_s0000_)]); xlabel('patients eliminated');
ylabel('negative-log-p');
end;%if flag_disp;
%%%%%%%%;

dir_mat_replication_new = sprintf('%s/dir_mat_replication_new',dir_trunk);
if ~exist(dir_mat_replication_new,'dir'); disp(sprintf(' %% mkdir %s',dir_mat_replication_new)); mkdir(dir_mat_replication_new); end;
dir_jpg_replication_new = sprintf('%s/dir_jpg_replication_new',dir_trunk);
if ~exist(dir_jpg_replication_new,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg_replication_new)); mkdir(dir_jpg_replication_new); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now as a test we use the bicluster to generate random-projections. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% First we reload the data. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ncontinent %d',ncontinent)); end;
if (flag_verbose>0); disp(sprintf(' %% ncontinent_1based %d',ncontinent_1based)); end;
if (flag_verbose>0); disp(sprintf(' %% dir_code: %s dir_trunk_Up05 %s;',dir_code,dir_trunk_Up05)); end;
if (flag_verbose>0); disp(sprintf(' %% str_lak_vs_dex: %s str_prefix %s;',str_lak_vs_dex,str_prefix)); end;
if (flag_verbose>0); disp(sprintf(' %% gamma %0.2f n_mds_0in %d n_mds_repl %d ij_mds_use_ [%d,%d];',gamma,n_mds_0in,n_mds_repl,ij_mds_use_)); end;
if (flag_verbose>0); disp(sprintf(' %% maf_lo_threshold %0.2f;',maf_lo_threshold)); end;
if (flag_verbose>0); disp(sprintf(' %% pca_rank_use %d;',pca_rank_use)); end;
if (flag_verbose>0); disp(sprintf(' %% ;')); end;
if (flag_verbose>0); disp(sprintf(' %% setting parameter_randproj_trn. ;')); end;
parameter_randproj_trn = struct('type','parameter');
parameter_randproj_trn.dir_code = dir_code;
parameter_randproj_trn.dir_trunk = dir_trunk_Up05;
parameter_randproj_trn.str_lak_vs_dex = str_lak_vs_dex;
parameter_randproj_trn.str_prefix = str_prefix;
parameter_randproj_trn.gamma = gamma;
parameter_randproj_trn.n_mds = n_mds_0in;
parameter_randproj_trn.n_mds_repl = n_mds_repl;
parameter_randproj_trn.ij_mds_use_ = ij_mds_use_;
parameter_randproj_trn.flag_force_create = 0; %<-- reload previous run. ;
parameter_randproj_trn.flag_verbose = max(0,flag_verbose-1);
parameter_randproj_trn.maf_lo_threshold = maf_lo_threshold;
parameter_randproj_trn.n_shuffle = 0; %<-- number of shuffles used to determine original trace: loaded below. ;
parameter_randproj_trn.dir_0in = sprintf('%s/dir_%s',parameter_randproj_trn.dir_trunk,parameter_randproj_trn.str_prefix);
fname_mr_A_full = sprintf('%s/%s_mr_A_full.b16',parameter_randproj_trn.dir_0in,parameter_randproj_trn.str_prefix);
fname_mr_Z_full = sprintf('%s/%s_mr_Z_full.b16',parameter_randproj_trn.dir_0in,parameter_randproj_trn.str_prefix);
fname_mc_A = sprintf('%s/%s_mc_A.b16',parameter_randproj_trn.dir_0in,parameter_randproj_trn.str_prefix);
mr_A_ori_trn_ = binary_uncompress(fname_mr_A_full)>0; %<-- thresholded to avoid negative values. ;
mr_Z_ori_trn_ = binary_uncompress(fname_mr_Z_full)>0; %<-- thresholded to avoid negative values. ;
mc_A_ori_trn_ = binary_uncompress(fname_mc_A)>0; %<-- thresholded to avoid negative values. ;
mr_A_alt_trn_ = mr_A_ori_trn_.*mr_ori_Up05_p01_continent_pc__(:,1+ncontinent); %<-- 'alternative' row-masks indicating continent-specific cases. ;
mr_Z_alt_trn_ = mr_Z_ori_trn_.*mr_ori_Up05_p01_continent_pc__(:,1+ncontinent); %<-- 'alternative' row-masks indicating continent-specific ctrls. ;
parameter_randproj_trn.str_mr_0in = sprintf('continent%d',ncontinent_1based); %<-- note that this is now 1based. ;
parameter_randproj_trn.nshuffle = 0;
parameter_randproj_trn = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_randproj_trn); %<-- ensure that bicluster has been calculated and fill out remaining parameters. ;
%%%%%%%%;
parameter_randproj_trn.str_out_xdrop_a_s0000 = sprintf('%s/out_xdrop_a.txt',parameter_randproj_trn.dir_out_s0000);
xdrop_trn_ = load_out_xdrop_from_str_ver0(parameter_randproj_trn.str_out_xdrop_a_s0000);
index_rdrop_trn_ = xdrop_trn_.index_rdrop_;
index_cdrop_trn_ = xdrop_trn_.index_cdrop_;
parameter_randproj_trn.dir_out_s0000_jpg = sprintf('%s/dir_out_jpg',parameter_randproj_trn.dir_out_s0000);
if ~exist(parameter_randproj_trn.dir_out_s0000_jpg,'dir'); disp(sprintf(' %% mkdir %s',parameter_randproj_trn.dir_out_s0000_jpg)); mkdir(parameter_randproj_trn.dir_out_s0000_jpg); end;
trace_trn____ = load_trace__from_dir_ver0(parameter_randproj_trn.dir_out_trace); %<-- load all the shuffled-traces. ;
n_shuffle = trace_trn____.n_shuffle;
if (flag_verbose>0); disp(sprintf(' %% found n_shuffle %d',n_shuffle)); end;
parameter_randproj_trn.n_shuffle = n_shuffle ;
%%%%%%%%;

%%%%%%%%;
if ncontinent==0; tmp_ij_nlpR = 1+21; end;
if ncontinent==1; tmp_ij_nlpR = 1+30; end;
tmp_ij_maxR = numel(trace_trn____.niter_s0000_); assert(tmp_ij_maxR==1+trace_trn____.niter_s0000_(end));
if (flag_verbose>0); disp(sprintf(' %% setting tmp_ij_nlpR 1+%d',tmp_ij_nlpR-1)); end;
r_rem_trn = trace_trn____.r_rem_s0000_(tmp_ij_nlpR);
r_max_trn = trace_trn____.r_rem_s0000_(1+0);
c_rem_trn = trace_trn____.c_rem_s0000_(tmp_ij_nlpR);
c_max_trn = trace_trn____.c_rem_s0000_(1+0);
index_r_rem_trn_ = xdrop_trn_.index_rkeep_(1:r_rem_trn);
index_r_max_trn_ = xdrop_trn_.index_rkeep_(1:r_max_trn);
index_c_rem_trn_ = xdrop_trn_.index_ckeep_(1:c_rem_trn);
index_c_max_trn_ = xdrop_trn_.index_ckeep_(1:c_max_trn);
assert(numel(index_r_rem_trn_)==numel(intersect(efind(mr_A_alt_trn_),index_r_rem_trn_)));
assert(numel(index_r_max_trn_)==numel(intersect(efind(mr_A_alt_trn_),index_r_max_trn_)));
assert(numel(index_c_rem_trn_)==numel(intersect(efind(mc_A_ori_trn_),index_c_rem_trn_)));
assert(numel(index_c_max_trn_)==numel(intersect(efind(mc_A_ori_trn_),index_c_max_trn_)));
%%%%%%%%;

flag_force_create = 0;
flag_disp=1;
n_shuffle_randproj = 200;
if ncontinent==0; n_shuffle_randproj = 200; end;
if ncontinent==1; n_shuffle_randproj = 500; end; %<-- just to be sure. ;
n_shuffle_interior = 64;
if (flag_verbose>0); disp(sprintf(' %% n_shuffle_randproj %d;',n_shuffle_randproj)); end;
parameter_randproj_trn.n_shuffle_randproj = n_shuffle_randproj; %<-- number of shuffles used for randproj (i.e., a postprocessing step). ;

str_datafile = sprintf('trnUp05_tst_Ap05_ncontinent%d',ncontinent_1based); %<-- note that in this filename ncontinent is 1-based rather than 0-based. ;
tmp_dir_new = sprintf('%s/dir_%s',dir_mat_replication_new,str_datafile);
if ~exist(tmp_dir_new,'dir'); disp(sprintf(' %% mkdir %s',tmp_dir_new)); mkdir(tmp_dir_new); end;
n_shuffle_randproj = parameter_randproj_trn.n_shuffle_randproj;
for nshuffle_randproj=0:n_shuffle_randproj-1+1;
tmp_fname_pre = sprintf('%s/randproj_trnUp05_tstAp05_ni%d_s%.4d',tmp_dir_new,tmp_ij_nlpR-1,nshuffle_randproj);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre);
if flag_force_create | ~tmp_flag_skip;
tmp_parameter_randproj_trn = parameter_randproj_trn;
pca_rank = pca_rank_use;
mr_A_rem_trn_ = 0.0*mr_A_alt_trn_; mr_A_rem_trn_(1+index_r_rem_trn_) = 1;
mr_Z_rem_trn_ = mr_Z_alt_trn_; %<-- limit the ctrls to those in ncontinent. ;
assert(sum(mr_A_rem_trn_)==sum(mr_A_rem_trn_.*mr_A_alt_trn_)); %<-- ensure that we only turned on rows in ncontinent. ;
assert(sum(mr_Z_rem_trn_)==sum(mr_Z_rem_trn_.*mr_Z_alt_trn_)); %<-- ensure that we only turned on rows in ncontinent. ;
mc_A_rem_trn_ = 0.0*mc_A_ori_trn_; mc_A_rem_trn_(1+index_c_rem_trn_) = 1;
assert(numel(efind(mc_A_rem_trn_))==numel(intersect(efind(mc_A_rem_trn_),index_c_max_trn_))); %<-- ensure that we only turned on cols allowed. ;
%%%%%%%%%%%%%%%%;
% If nshuffle_randproj> 0, define random randproj. ;
%%%%%%%%%%%%%%%%;
mr_A_randproj_trn_ = mr_A_rem_trn_;
mr_Z_randproj_trn_ = mr_Z_rem_trn_;
mc_A_randproj_trn_ = mc_A_rem_trn_;
if nshuffle_randproj>0;
rng(nshuffle_randproj);
tmp_mr_AZ_ = mr_A_alt_trn_ + mr_Z_alt_trn_;
tmp_index_AZ_ = efind(tmp_mr_AZ_);
tmp_p_ = randperm(numel(tmp_index_AZ_));
tmp_index_r_A_ = tmp_index_AZ_(tmp_p_(1:numel(index_r_rem_trn_)));
tmp_index_r_Z_ = tmp_index_AZ_(tmp_p_(numel(index_r_rem_trn_)+[1:numel(efind(mr_Z_rem_trn_))]));
mr_A_randproj_trn_ = 0.0*mr_A_alt_trn_; mr_A_randproj_trn_(1+tmp_index_r_A_) = 1;
mr_Z_randproj_trn_ = 0.0*mr_Z_alt_trn_; mr_Z_randproj_trn_(1+tmp_index_r_Z_) = 1;
assert(sum(mr_A_randproj_trn_)==sum(mr_A_randproj_trn_.*tmp_mr_AZ_)); %<-- ensure that we only turned on rows in ncontinent. ;
assert(sum(mr_Z_randproj_trn_)==sum(mr_Z_randproj_trn_.*tmp_mr_AZ_)); %<-- ensure that we only turned on rows in ncontinent. ;
assert(sum(mr_A_randproj_trn_.*mr_Z_randproj_trn_)==0); %<-- ensure that the cases and ctrls are distinct. ;
clear tmp_mr_AZ_ tmp_index_AZ_ tmp_p_ tmp_index_r_A_ tmp_index_r_Z_ ;
tmp_mc_ = mc_A_rem_trn_;
tmp_index_ = efind(tmp_mc_);
tmp_p_ = randperm(numel(tmp_index_));
tmp_index_c_ = tmp_index_(tmp_p_(1:numel(index_c_rem_trn_)));
mc_A_randproj_trn_ = 0.0*mc_A_ori_trn_; mc_A_randproj_trn_(1+tmp_index_c_) = 1;
assert(numel(efind(mc_A_randproj_trn_))==numel(intersect(efind(mc_A_randproj_trn_),index_c_max_trn_))); %<-- ensure that we only turned on cols allowed. ;
clear tmp_mc_ tmp_index_ tmp_p_ tmp_index_c_ ;
end;%if nshuffle_randproj>0;
if (flag_verbose>0);
disp(sprintf(' %% mr_A_rem_trn_: %0.4d mr_A_randproj_trn_: %0.4d intersection: %0.4d',numel(efind(mr_A_rem_trn_)),numel(efind(mr_A_randproj_trn_)),numel(efind(mr_A_rem_trn_ & mr_A_randproj_trn_))));
disp(sprintf(' %% mr_Z_rem_trn_: %0.4d mr_Z_randproj_trn_: %0.4d intersection: %0.4d',numel(efind(mr_Z_rem_trn_)),numel(efind(mr_Z_randproj_trn_)),numel(efind(mr_Z_rem_trn_ & mr_Z_randproj_trn_))));
tmp_cap__ = [ ...
                               0 , numel(efind(mr_A_rem_trn_))                      ,                      numel(efind(mr_A_randproj_trn_)) ; ...
     numel(efind(mr_Z_rem_trn_)) , numel(efind(mr_Z_rem_trn_ & mr_A_rem_trn_))      ,      numel(efind(mr_Z_rem_trn_ & mr_A_randproj_trn_)) ; ...
numel(efind(mr_Z_randproj_trn_)) , numel(efind(mr_Z_randproj_trn_ & mr_A_rem_trn_)) , numel(efind(mr_Z_randproj_trn_ & mr_A_randproj_trn_)) ; ...
];
disp(sprintf(' %% intersection matrix: (rows = mr_Z_rem_trn_ and mr_Z_randproj_trn_) (cols = mr_A_rem_trn_ and mr_A_randproj_trn_)')); disp(num2str(tmp_cap__,'%0.4d '));
end;%if (flag_verbose>0);
%%%%%%%%%%%%%%%%;
% First find principal-components within Up05. ;
%%%%%%%%%%%%%%%%;
pca_mr_A_randproj_trn_ = { mr_A_randproj_trn_ };
pca_mr_Z_randproj_trn_ = { mr_Z_randproj_trn_ };
pca_mc_A_randproj_trn = mc_A_randproj_trn_; pca_mc_A_randproj_trn(1+index_Up05_setminus_cap_) = 0; %<-- limited to allele-combinations in allele-cap. ;
pca_str_infix_randproj_trn=sprintf('randproj_trnUp05_tstAp05_ni%d_s%.4d_p01',tmp_ij_nlpR-1,nshuffle_randproj);
tmp_parameter_randproj_trn.flag_force_create = 0;
parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01 = tmp_parameter_randproj_trn;
parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01.str_A_p = str_Up05_A_p_p01;
parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01.str_name_s0000 = sprintf('pca_D_randproj_trnUp05_tstAp05_ni%d_s%.4d_p01',tmp_ij_nlpR-1,nshuffle_randproj);
[ ...
 parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01 ...
,AZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_ ...
,AnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_ ...
,ZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_ ...
,V_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01 ...
,pca_rank ...
,pca_mr_A_randproj_trn_ ...
,pca_mr_Z_randproj_trn_ ...
,pca_mc_A_randproj_trn ...
,pca_str_infix_randproj_trn ...
);
%%%%%%%%;
mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_ = 0.0*mr_A_ori_trn_;
mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_(1+efind(mr_Z_randproj_trn_)) = 1; %<-- color ctrl-patients. ;
mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_(1+efind(mr_A_randproj_trn_)) = 2; %<-- color case-patients. ;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_sml = 12;
hold on;
tmp_index_ = efind( (mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_==1) | (mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_==2) );
scatter(AZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_(1+tmp_index_,1+0),AZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_(1+tmp_index_,1+1),markersize_use,mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_(1+tmp_index_),'filled','MarkerEdgeColor','k');
axis equal; grid on;
title('AZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
end;%if (flag_disp>1);
%%%%%%%%;
% Extract principal-components from above calculation, then prepare for Ap05 projection . ;
%%%%%%%%;
pca_str_infix_randproj_trn=sprintf('randproj_trnUp05_tstAp05_ni%d_s%.4d_p01',tmp_ij_nlpR-1,nshuffle_randproj); %<-- this defines the pca file to use. ;
parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01.str_V_ = sprintf('%s/dir_pca/dir_pca_mda/pca_D_%s_k%d_B44_V_.mda',parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01.dir_out_s0000,pca_str_infix_randproj_trn,pca_rank);
if ~exist(parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01.str_V_,'file');
disp(sprintf(' %% Warning, %s not found',parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01.str_V_));
end;%if ~exist(parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01.str_V_,'file');
D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_ = mda_read_r8(parameter_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01.str_V_);
if (flag_verbose); disp(sprintf(' %% size of intersection: %d',numel(index_Up05_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% Support of V_: %d',sum(abs(D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_(:,1))>0))); end;
if (flag_verbose); disp(sprintf(' %% intersection of supp(V_) with cap: %d',numel(intersect(index_Up05_from_cap_,efind(abs(D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_(:,1))>0))))); end;
D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_ = ij_Ap05_from_Up05__*D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_;
if (flag_verbose); disp(sprintf(' %% fnorm(D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_(ij_Ap05_from_cap_,:) - D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_(ij_Up05_from_cap_,:)): %0.16f',fnorm(D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_(ij_Ap05_from_cap_,:) - D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_(ij_Up05_from_cap_,:)))); end;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1); hist(D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_(1+efind(D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_~=0)),128);
title('D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_','Interpreter','none');
subplot(1,2,2); hist(D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_(1+efind(D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_~=0)),128);
title('D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_','Interpreter','none');
end;%if (flag_disp>1);
%%%%%%%%;
% Now project Ap05 data onto D_Ap05_from_randproj_trnUp05_tstAp05_nix_p01_V_. ;
%%%%%%%%;
ndataset=ndataset_Ap05;
mx_tst__ = dataset_{1+ndataset}.mx__;
tmp_parameter_tst = dataset_{1+ndataset}.parameter;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.dir_0in: %s',tmp_parameter_tst.dir_0in)); end;
if (flag_verbose); disp(sprintf(' %% tmp_parameter_tst.str_prefix: %s',tmp_parameter_tst.str_prefix)); end;
tmp_parameter_tst.dir_code = sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev',str_home);
tmp_parameter_tst.maf_lo_threshold = maf_lo_threshold;
tmp_parameter_tst.maf_hi_threshold = 0.50;
tmp_parameter_tst.gamma = 0.05; %<-- not used for pca calculation. ;
tmp_parameter_tst.flag_verbose = 0;
tmp_parameter_tst.flag_force_create = 0;
pca_rank = pca_rank_use;
mr_A_ori_tst_ = mx_tst__.mr_A_full_;
mr_Z_ori_tst_ = mx_tst__.mr_Z_full_;
mc_A_ori_tst_ = mx_tst__.mc_A_;
mr_A_alt_tst_ = mr_A_ori_tst_.*mr_ori_Ap05_p01_continent_pc__(:,1+ncontinent);
mr_Z_alt_tst_ = mr_Z_ori_tst_.*mr_ori_Ap05_p01_continent_pc__(:,1+ncontinent);
mc_A_alt_tst_ = mc_A_ori_tst_;
pca_mr_A_tst_ = { mr_A_alt_tst_ };
pca_mr_Z_tst_ = { mr_Z_alt_tst_ };
pca_mc_A_tst = mc_A_alt_tst_; pca_mc_A_tst(1+index_Ap05_setminus_cap_) = 0;
pca_str_infix_tst=sprintf('D_Ap05_from_randproj_trnUp05_tstAp05_ni%d_s%.4d_p01',tmp_ij_nlpR-1,nshuffle_randproj);
parameter_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01 = tmp_parameter_tst;
parameter_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01.str_A_p = str_Ap05_A_p_p01;
parameter_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01.str_name_s0000 = ...
  sprintf('pca_D_Ap05c%d_from_trnUp05c%d_tstAp05c%d_nix_p01',ncontinent,ncontinent,ncontinent);
[ ...
 parameter_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01 ...
,AZnV_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_ ...
,AnV_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_ ...
,ZnV_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_ ...
,tmp_V_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_proj_from_mx_ver16( ...
 parameter_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01 ...
,pca_rank ...
,D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_V_ ...
,pca_mr_A_tst_ ...
,pca_mr_Z_tst_ ...
,pca_mc_A_tst ...
,pca_str_infix_tst ...
);
%%%%;
mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_ = 0.0*mr_A_ori_tst_;
mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_(1+efind(mr_Z_alt_tst_)) = 1; %<-- color ctrl-patients. ;
mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_(1+efind(mr_A_alt_tst_)) = 2; %<-- color case-patients. ;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_use = 12;
c_ctrl_ = [0,1,1];
c_case_ = [1,0,1];
hold on;
tmp_index_ = efind( (mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_==1) );
plot(AZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_(1+tmp_index_,1+0),AZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_(1+tmp_index_,1+1),'o','MarkerSize',markersize_use,'MarkerFaceColor',c_ctrl_,'MarkerEdgeColor','k');
tmp_index_ = efind( (mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_==2) );
plot(AZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_(1+tmp_index_,1+0),AZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_(1+tmp_index_,1+1),'o','MarkerSize',markersize_use,'MarkerFaceColor',c_case_,'MarkerEdgeColor','k');
tmp_index_ = efind(mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_==1);
plot(AZnV_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_(1+tmp_index_,1+0),AZnV_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_(1+tmp_index_,1+1),'^','MarkerSize',markersize_use,'MarkerFaceColor',c_ctrl_,'MarkerEdgeColor','k');
tmp_index_ = efind(mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_==2);
plot(AZnV_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_(1+tmp_index_,1+0),AZnV_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_(1+tmp_index_,1+1),'^','MarkerSize',markersize_use,'MarkerFaceColor',c_case_,'MarkerEdgeColor','k');
axis equal; grid on;
title('AZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
end;%if (flag_disp>1);
%%%%%%%%;
% Now align. ;
%%%%%%%%;
tmp_index_Up05_ = efind( (mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_==1) | (mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_==2) | (mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_==3) );
tmp_prm_Up05_ = transpose(1:numel(tmp_index_Up05_));
tmp_index_Ap05_ = efind( (mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_==1) | (mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_==2) );
tmp_prm_Ap05_ = transpose(1:numel(tmp_index_Ap05_));
if (nshuffle_randproj> 0);
if (flag_verbose>0); disp(sprintf(' %% shuffling with rng(%.4d)',nshuffle_randproj)); end;
rng(nshuffle_randproj); tmp_prm_Up05_ = randperm(numel(tmp_prm_Up05_)); tmp_prm_Ap05_ = randperm(numel(tmp_prm_Ap05_));
end;%if (nshuffle_randproj> 0);
tmp_mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_ = mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_ ;
tmp_mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_(1+tmp_index_Up05_) = tmp_mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_(1+tmp_index_Up05_(tmp_prm_Up05_)) ;
tmp_mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_ = mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_;
tmp_mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_(1+tmp_index_Ap05_) = tmp_mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_(1+tmp_index_Ap05_(tmp_prm_Ap05_)) ;
parameter_apm = struct('type','parameter_apm');
parameter_apm.k_use = 32; %<-- This is the initial number of nearest neighbors to use. ;
if ncontinent==0; k_gamma_use = k_gamma_use_ncontinent0; end;
if ncontinent==1; k_gamma_use = k_gamma_use_ncontinent1; end;
parameter_apm.k_gamma = k_gamma_use ;
parameter_apm.flag_disp = 0;
parameter_apm.n_shuffle = n_shuffle_interior;
[ ...
 parameter_apm ...
,tmp_z_z_ ...
,tmp_z_zp__ ...
,tmp_z_x_ ...
,tmp_z_xp__ ...
,tmp_z_y_ ...
,tmp_z_yp__ ...
,tmp_f_z_ ...
,tmp_f_zp__ ...
,tmp_f_x_ ...
,tmp_f_xp__ ...
,tmp_f_y_ ...
,tmp_f_yp__ ...
,tmp_a_est_ ...
,tmp_A_est__ ...
] = ...
test_Up05_vs_Ap05_18c_helper_0_helper_0( ...
 parameter_apm ...
,tmp_mr_dvx_randproj_trnUp05_tstAp05_nix_sxxxx_ ...
,AZnV_D_randproj_trnUp05_tstAp05_nix_sxxxx_p01_ ...
,tmp_mr_dvx_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_ ...
,AZnV_D_Ap05_from_randproj_trnUp05_tstAp05_nix_sxxxx_p01_ ...
);
save(tmp_fname_mat ...
     ,'nshuffle_randproj' ...
     ,'parameter_apm' ...
     ,'tmp_z_z_' ...
     ,'tmp_z_zp__' ...
     ,'tmp_z_x_' ...
     ,'tmp_z_xp__' ...
     ,'tmp_z_y_' ...
     ,'tmp_z_yp__' ...
     ,'tmp_f_z_' ...
     ,'tmp_f_zp__' ...
     ,'tmp_f_x_' ...
     ,'tmp_f_xp__' ...
     ,'tmp_f_y_' ...
     ,'tmp_f_yp__' ...
     ,'tmp_a_est_' ...
     ,'tmp_A_est__' ...
     );
close_fname_tmp(tmp_fname_pre);
end;%if flag_force_create | ~tmp_flag_skip;
end;%for nshuffle_randproj=0:n_shuffle_randproj-1+1;

%%%%%%%%;
% Now collate randproj data and determine typical variance as function of fraction nearest. ;
%%%%%%%%;
mr_A_alt_tst_ = mr_A_ori_tst_.*mr_ori_Ap05_p01_continent_pc__(:,1+ncontinent);
mr_Z_alt_tst_ = mr_Z_ori_tst_.*mr_ori_Ap05_p01_continent_pc__(:,1+ncontinent);
n_x = sum(mr_A_alt_tst_ + mr_Z_alt_tst_);
n_y = numel(index_r_rem_trn_) + sum(mr_Z_alt_trn_>0);
n_z = max(n_x,n_y);
%z_zr__ = zeros(n_z,1+n_shuffle_randproj);
%z_zpr___ = zeros(n_z,n_shuffle_interior,1+n_shuffle_randproj);
%z_xr__ = zeros(n_x,1+n_shuffle_randproj);
%z_xpr___ = zeros(n_x,n_shuffle_interior,1+n_shuffle_randproj);
%z_yr__ = zeros(n_y,1+n_shuffle_randproj);
%z_ypr___ = zeros(n_y,n_shuffle_interior,1+n_shuffle_randproj);
f_zr__ = zeros(n_z,1+n_shuffle_randproj);
f_zpr___ = zeros(n_z,n_shuffle_interior,1+n_shuffle_randproj);
f_xr__ = zeros(n_x,1+n_shuffle_randproj);
f_xpr___ = zeros(n_x,n_shuffle_interior,1+n_shuffle_randproj);
f_yr__ = zeros(n_y,1+n_shuffle_randproj);
f_ypr___ = zeros(n_y,n_shuffle_interior,1+n_shuffle_randproj);
found_flag_r_ = zeros(1+n_shuffle_randproj);
%%%%;
nl=0;
for nshuffle_randproj=0:n_shuffle_randproj;
tmp_fname_pre = sprintf('%s/randproj_trnUp05_tstAp05_ni%d_s%.4d',tmp_dir_new,tmp_ij_nlpR-1,nshuffle_randproj);
tmp_fname_mat = sprintf('%s.mat',tmp_fname_pre);
if ~exist(tmp_fname_mat,'file'); if (flag_verbose>0); disp(sprintf(' %% %s not found, skipping',tmp_fname_mat)); end; end;
if  exist(tmp_fname_mat,'file');
if (flag_verbose>1); disp(sprintf(' %% %s found, not skipping',tmp_fname_mat)); end;
found_flag_r_(1+nshuffle_randproj) = 1;
tmp_randproj_ = load(tmp_fname_mat);
%z_zr__(:,1+nl) = tmp_randproj_.tmp_z_z_;
%z_zpr___(:,:,1+nl) = tmp_randproj_.tmp_z_zp__;
%z_xr__(:,1+nl) = tmp_randproj_.tmp_z_x_;
%z_xpr___(:,:,1+nl) = tmp_randproj_.tmp_z_xp__;
%z_yr__(:,1+nl) = tmp_randproj_.tmp_z_y_;
%z_ypr___(:,:,1+nl) = tmp_randproj_.tmp_z_yp__;
f_zr__(:,1+nl) = tmp_randproj_.tmp_f_z_;
f_zpr___(:,:,1+nl) = tmp_randproj_.tmp_f_zp__;
f_xr__(:,1+nl) = tmp_randproj_.tmp_f_x_;
f_xpr___(:,:,1+nl) = tmp_randproj_.tmp_f_xp__;
f_yr__(:,1+nl) = tmp_randproj_.tmp_f_y_;
f_ypr___(:,:,1+nl) = tmp_randproj_.tmp_f_yp__;
nl=nl+1;
clear tmp_randproj_ ;
end;%if  exist(tmp_fname_mat,'file');
end;%for nshuffle_randproj=0:n_shuffle_randproj;
n_l = nl;
f_zr__ = f_zr__(:,1:n_l);;
f_zpr___ = f_zpr___(:,:,1:n_l);;
f_xr__ = f_xr__(:,1:n_l);;
f_xpr___ = f_xpr___(:,:,1:n_l);;
f_yr__ = f_yr__(:,1:n_l);;
f_ypr___ = f_ypr___(:,:,1:n_l);;
%%%%%%%%;
f_avg_z_ = mean(f_zr__(:,2:end),2);
f_std_z_ = std(f_zr__(:,2:end),1,2);
z_zr__ = bsxfun(@rdivide,bsxfun(@minus,f_zr__,mean(f_zr__(:,2:end),2)),max(1e-12,std(f_zr__(:,2:end),1,2)));
[~,tmp_ij_zr__] = sort(f_zr__,2,'ascend');
[~,prctile_zr__] = sort(tmp_ij_zr__,2,'ascend');
prctile_zr__ = prctile_zr__/size(prctile_zr__,2);
%%%%;
f_avg_x_ = mean(f_xr__(:,2:end),2);
f_std_x_ = std(f_xr__(:,2:end),1,2);
z_xr__ = bsxfun(@rdivide,bsxfun(@minus,f_xr__,mean(f_xr__(:,2:end),2)),max(1e-12,std(f_xr__(:,2:end),1,2)));
[~,tmp_ij_xr__] = sort(f_xr__,2,'ascend');
[~,prctile_xr__] = sort(tmp_ij_xr__,2,'ascend');
prctile_xr__ = prctile_xr__/size(prctile_xr__,2);
%%%%;
f_avg_y_ = mean(f_yr__(:,2:end),2);
f_std_y_ = std(f_yr__(:,2:end),1,2);
z_yr__ = bsxfun(@rdivide,bsxfun(@minus,f_yr__,mean(f_yr__(:,2:end),2)),max(1e-12,std(f_yr__(:,2:end),1,2)));
[~,tmp_ij_yr__] = sort(f_yr__,2,'ascend');
[~,prctile_yr__] = sort(tmp_ij_yr__,2,'ascend');
prctile_yr__ = prctile_yr__/size(prctile_yr__,2);
%%%%%%%%;
f_avg_z_ = mean(f_zpr___(:,2:end,2:end),[2,3]);
f_std_z_ = std(f_zpr___(:,2:end,2:end),1,[2,3]);
z_zr__ = bsxfun(@rdivide,bsxfun(@minus,f_zr__,f_avg_z_),max(1e-12,f_std_z_));
z_zpr___ = bsxfun(@rdivide,bsxfun(@minus,f_zpr___,f_avg_z_),max(1e-12,f_std_z_));
%%%%;
f_avg_x_ = mean(f_xpr___(:,2:end,2:end),[2,3]);
f_std_x_ = std(f_xpr___(:,2:end,2:end),1,[2,3]);
z_xr__ = bsxfun(@rdivide,bsxfun(@minus,f_xr__,f_avg_x_),max(1e-12,f_std_x_));
z_xpr___ = bsxfun(@rdivide,bsxfun(@minus,f_xpr___,f_avg_x_),max(1e-12,f_std_x_));
%%%%;
f_avg_y_ = mean(f_ypr___(:,2:end,2:end),[2,3]);
f_std_y_ = std(f_ypr___(:,2:end,2:end),1,[2,3]);
z_yr__ = bsxfun(@rdivide,bsxfun(@minus,f_yr__,f_avg_y_),max(1e-12,f_std_y_));
z_ypr___ = bsxfun(@rdivide,bsxfun(@minus,f_ypr___,f_avg_y_),max(1e-12,f_std_y_));
%%%%;
tmp_envelope_randproj_z_ = squeeze(mean(std(f_zpr___,1,2),3));
%%%%%%%%;

%%%%%%%%;
% Now collate the results. ;
%%%%%%%%;
mr_A_alt_tst_ = mr_A_ori_tst_.*mr_ori_Ap05_p01_continent_pc__(:,1+ncontinent);
mr_Z_alt_tst_ = mr_Z_ori_tst_.*mr_ori_Ap05_p01_continent_pc__(:,1+ncontinent);
n_x_old = sum(mr_A_alt_tst_ + mr_Z_alt_tst_);
n_y_old = sum(mr_A_alt_trn_>0) + sum(mr_Z_alt_trn_>0);
n_z_old = max(n_x_old,n_y_old);
n_shuffle_exterior = 500; %<-- used by jeremy. ;
f_ze__ = zeros(n_z_old,1+n_shuffle_exterior);
f_xe__ = zeros(n_x_old,1+n_shuffle_exterior);
f_ye__ = zeros(n_y_old,1+n_shuffle_exterior);
tmp_dir_old = sprintf('%s/dir_%s',dir_mat_replication_old,str_datafile);
if ~exist(tmp_dir_old,'dir'); disp(sprintf(' %% mkdir %s',tmp_dir_old)); mkdir(tmp_dir_old); end;
for nshuffle_exterior=0:n_shuffle_exterior;
tmp_fname_pre = sprintf('%s/%s_s%.4d',tmp_dir_old,str_datafile,nshuffle_exterior);
tmp_fname_mat = sprintf('%s.mat',tmp_fname_pre);
if ~exist(tmp_fname_mat,'file');
if (flag_verbose); disp(sprintf(' %% Warning, %s not found',tmp_fname_mat)); end;
end;%if ~exist(tmp_fname_mat,'file');
if  exist(tmp_fname_mat,'file');
tmp_tmp_ = load(tmp_fname_mat);
f_ze__(:,1+nshuffle_exterior) = tmp_tmp_.tmp_f_z_;
f_xe__(:,1+nshuffle_exterior) = tmp_tmp_.tmp_f_x_;
f_ye__(:,1+nshuffle_exterior) = tmp_tmp_.tmp_f_y_;
clear tmp_tmp_;
end;%if  exist(tmp_fname_mat,'file');
end;%for nshuffle_exterior=0:n_shuffle_exterior;
%%%%%%%%;
f_avg_z_ = mean(f_ze__(:,2:end),2);
f_std_z_ = std(f_ze__(:,2:end),1,2);
z_ze__ = bsxfun(@rdivide,bsxfun(@minus,f_ze__,mean(f_ze__(:,2:end),2)),max(1e-12,std(f_ze__(:,2:end),1,2)));
[~,tmp_ij_ze__] = sort(f_ze__,2,'ascend');
[~,prctile_ze__] = sort(tmp_ij_ze__,2,'ascend');
prctile_ze__ = prctile_ze__/size(prctile_ze__,2);
%%%%;
f_avg_x_ = mean(f_xe__(:,2:end),2);
f_std_x_ = std(f_xe__(:,2:end),1,2);
z_xe__ = bsxfun(@rdivide,bsxfun(@minus,f_xe__,mean(f_xe__(:,2:end),2)),max(1e-12,std(f_xe__(:,2:end),1,2)));
[~,tmp_ij_xe__] = sort(f_xe__,2,'ascend');
[~,prctile_xe__] = sort(tmp_ij_xe__,2,'ascend');
prctile_xe__ = prctile_xe__/size(prctile_xe__,2);
%%%%;
f_avg_y_ = mean(f_ye__(:,2:end),2);
f_std_y_ = std(f_ye__(:,2:end),1,2);
z_ye__ = bsxfun(@rdivide,bsxfun(@minus,f_ye__,mean(f_ye__(:,2:end),2)),max(1e-12,std(f_ye__(:,2:end),1,2)));
[~,tmp_ij_ye__] = sort(f_ye__,2,'ascend');
[~,prctile_ye__] = sort(tmp_ij_ye__,2,'ascend');
prctile_ye__ = prctile_ye__/size(prctile_ye__,2);
%%%%%%%%;
fname_fig_pre = sprintf('%s/%s_FIGD',dir_jpg_replication_new,str_datafile);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,768]);
fontsize_use = 14;
p_row = 1; p_col = 3; ntab = 0;
linewidth_sml = 0.5; linewidth_big = 2.0; ylim_ = 10*[-1,+1];
tmp_x_old_ = linspace(0,1,n_x_old); tmp_y_old_ = linspace(0,1,n_y_old); tmp_z_old_ = linspace(0,1,n_z_old);
tmp_w_ = tmp_z_old_; f_we__ = f_ze__; z_we__ = z_ze__; prctile_we__ = prctile_ze__; str_f_w = 'f_z_'; str_z_w = 'z_z_'; str_p_w = 'p_z_'; 
%%%%;
subplot(p_row,p_col,1+ntab+0);
hold on;
plot(tmp_w_(tmp_w_<=0.5),f_we__(tmp_w_<=0.5,2:end),'k','LineWidth',linewidth_sml,'Color',[0 0 0 0.1]);
plot(tmp_w_(tmp_w_<=0.5),f_we__(tmp_w_<=0.5,1+0)  ,'r','LineWidth',linewidth_big);
hold off;
% title(str_f_w,'Interpreter','none');
xlim([0,.50]); 
set(gca,'XTick',0:0.05:.50,'FontSize',fontsize_use); grid on; xtickangle(90);ytickformat('%.2f'),xtickformat('%.2f')
xlabel('fraction nearest','FontSize',18);
ylabel('f fraction','FontSize',18);
%%%%;
subplot(p_row,p_col,1+ntab+1);
hold on;
plot(tmp_w_(tmp_w_<=0.5),max(min(ylim_),min(max(ylim_),z_we__(tmp_w_<=0.5,2:end))),'k','LineWidth',linewidth_sml,'Color',[0 0 0 0.1]);
plot(tmp_w_(tmp_w_<=0.5),max(min(ylim_),min(max(ylim_),z_we__(tmp_w_<=0.5,1+0)  )),'r','LineWidth',linewidth_big);
hold off;
% title(str_z_w,'Interpreter','none');
xlim([0,.50]); 
set(gca,'XTick',0:0.05:.50,'FontSize',fontsize_use); grid on; xtickangle(90);ytickformat('%.2f'),xtickformat('%.2f')
xlabel('fraction nearest','FontSize',18);
ylim(ylim_); ylabel('z score','FontSize',18);
%%%%;
subplot(p_row,p_col,1+ntab+2);
hold on;
yline(0.95, '--','LineWidth',2,'Color',[0.3 0.3 0.3 .3]) % add a horizontal line .95
plot(tmp_w_(tmp_w_<=0.5),prctile_we__(tmp_w_<=0.5,1+0)  ,'r','LineWidth',linewidth_big);
hold off;
% title(str_p_w,'Interpreter','none');
xlim([0,.50]); 
set(gca,'YTick',0:0.05:1.0,'FontSize',fontsize_use); set(gca,'XTick',0:0.05:.50,'FontSize',fontsize_use); grid on; xtickangle(90); ytickformat('%.2f'),xtickformat('%.2f')
set(gca,'TickLength',[0,0]);
ylim([0,1]); 
xlabel('fraction nearest','FontSize',18);
ylabel('1-p empirical','FontSize',18);
%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
if flag_replot | ~exist(fname_fig_jpg) ; print('-djpeg',fname_fig_jpg); end;
if flag_replot | ~exist(fname_fig_eps) ; print('-depsc',fname_fig_eps); end;
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

tmp_z_ = linspace(0,1,n_z);
tmp_z_old_ = linspace(0,1,n_z_old);
tmp_envelope_randproj_z_ = squeeze(mean(std(f_zpr___(:,2:end,2:end),1,2),3));
tmp_envelope_randproj_old_z_ = reshape(interp1(tmp_z_,tmp_envelope_randproj_z_,tmp_z_old_),[n_z_old,1]);
z_corrected_ze__ = bsxfun(@rdivide,bsxfun(@minus,z_ze__,mean(z_ze__(:,2:end),2)),tmp_envelope_randproj_old_z_);
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
hold on;
index_old_cut_ = efind(tmp_z_old_>=0.05 & tmp_z_old_<=0.50);
tmp_z_old_cut_ = tmp_z_old_(1+index_old_cut_);
for ns=flip(0:500);
c_use_ = [0,0,0]; if (ns==0); c_use_ = [1,0,0]; end;
tmp_z_z_ = f_ze__(1+index_old_cut_,1+ns);
[tmp_z,tmp_ij] = max(tmp_z_z_);
plot(tmp_z_old_cut_,tmp_z_z_,'-','LineWidth',0.25,'Color',c_use_);
plot(tmp_z_old_cut_(tmp_ij),tmp_z,'o','Color',c_use_);
end;%for ns=0:500;
hold off;
xlabel('fraction nearest');
ylabel('z-score');
title('not corrected');
subplot(1,2,2);
hold on;
index_old_cut_ = efind(tmp_z_old_>=0.05 & tmp_z_old_<=0.50);
tmp_z_old_cut_ = tmp_z_old_(1+index_old_cut_);
for ns=flip(0:500);
c_use_ = [0,0,0]; if (ns==0); c_use_ = [1,0,0]; end;
tmp_z_z_ = z_corrected_ze__(1+index_old_cut_,1+ns);
[tmp_z,tmp_ij] = max(tmp_z_z_);
plot(tmp_z_old_cut_,tmp_z_z_,'-','LineWidth',0.25,'Color',c_use_);
plot(tmp_z_old_cut_(tmp_ij),tmp_z,'o','Color',c_use_);
end;%for ns=0:500;
hold off;
xlabel('fraction nearest');
ylabel('z-score');
title('yes corrected');
end;%if flag_disp;
%%%%%%%%;

index_old_cut_ = efind(tmp_z_old_>=0.05 & tmp_z_old_<=0.50);
tmp_z_max_e_ = max(z_ze__(1+index_old_cut_,:),[],1);
p_max_not_corrected = numel(efind(tmp_z_max_e_(2:end)>tmp_z_max_e_(1+0)))/numel(tmp_z_max_e_);
tmp_z_avg_e_ = mean(z_ze__(1+index_old_cut_,:),1);
p_avg_not_corrected = numel(efind(tmp_z_avg_e_(2:end)>tmp_z_avg_e_(1+0)))/numel(tmp_z_avg_e_);
tmp_z_corrected_max_e_ = max(z_corrected_ze__(1+index_old_cut_,:),[],1);
p_max_yes_corrected = numel(efind(tmp_z_corrected_max_e_(2:end)>tmp_z_corrected_max_e_(1+0)))/numel(tmp_z_corrected_max_e_);
tmp_z_corrected_avg_e_ = mean(z_corrected_ze__(1+index_old_cut_,:),1);
p_avg_yes_corrected = numel(efind(tmp_z_corrected_avg_e_(2:end)>tmp_z_corrected_avg_e_(1+0)))/numel(tmp_z_corrected_avg_e_);

if (flag_verbose>0);
disp(sprintf(' %% ncontinent %d ncontinent_1based %d',ncontinent,ncontinent_1based));
disp(sprintf(' %% p_max_not_corrected: %0.6f',p_max_not_corrected));
disp(sprintf(' %% p_avg_not_corrected: %0.6f',p_avg_not_corrected));
disp(sprintf(' %% p_max_yes_corrected: %0.6f',p_max_yes_corrected));
disp(sprintf(' %% p_avg_yes_corrected: %0.6f',p_avg_yes_corrected));
end;%if (flag_verbose>0);


