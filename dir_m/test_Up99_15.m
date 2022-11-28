%%%%%%%%;
% testing Up99. ;
% setting ent_cutoff and mss_cutoff to 0.15. ;
% This should be sufficiently large that all snps and patients are retained. ;
%%%%%%%%;

clear;
% platform = 'rusty';
% if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
% if (strcmp(platform,'access1')); str_home = 'data'; end;
% if (strcmp(platform,'OptiPlex')); str_home = 'home'; end;
% if (strcmp(platform,'eval1')); str_home = 'home'; end;
% if (strcmp(platform,'rusty')); str_home = 'mnt/home'; end;
%%%%%%%%;
% run(sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev/dir_m/setup_0',str_home));
verbose = 1;
flag_disp = 1+verbose; nf=0;
flag_replot = 0;
if (verbose); disp(sprintf(' %% ;')); end;
if (verbose); disp(sprintf(' %% Assume path is set using dir_lakcluster_c/dir_m/setup_0.m. ;')); end;

dir_trunk = '/home/jelman/Projects/AD_Biclustering/data/UKB/ukb_pca_p05-p1';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

parameter = struct('type','parameter');
parameter.str_output_prefix = 'Up99';
parameter.dir_trunk = dir_trunk;
parameter.flag_verbose = 1;
parameter.flag_crop = 0;
parameter.maf_cutoff = 0.01; %<-- standard maf cutoff. ;
parameter.ent_cutoff = 0.03; %<-- standard entropy cutoff. ;
parameter.mss_cutoff = 0.10; %<-- standard missingness cutoff. ;
parameter.flag_force_create = 1;
parameter.flag_stop = 0;
parameter.slurm_memdecl = 128; 
n_study = 1;
study_name_of_branch_s_ = {'dir_Up99'};
study_name_without_extension_s_ = {'ukb_imp_allchr_ADCaseControls_p1'};
n_mds_p = 8217;
n_mds_v = 2;
fname_mds_without_extension = sprintf('ukb_imp_allchr_ADCaseControls_MDS2.tsv');
tmp_dir = sprintf('%s',dir_trunk);
mds_pv__ = textread(sprintf('%s/%s',tmp_dir,fname_mds_without_extension));
mds_fidandiid_p_ = [];
n_famex = [];
famex_fidandiid_p_ = [];

[ ...
 parameter ...
] = ...
bed_to_b16_flip_ver7( ...
 parameter ...
,n_study ...
,study_name_of_branch_s_ ...
,study_name_without_extension_s_ ...
,n_mds_p ...
,n_mds_v ...
,mds_pv__ ...
,mds_fidandiid_p_ ...
,n_famex ...
,famex_fidandiid_p_ ...
);
