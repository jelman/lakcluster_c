clear;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); str_home = 'data'; end;
if (strcmp(platform,'OptiPlex')); str_home = 'home'; end;
if (strcmp(platform,'eval1')); str_home = 'home'; end;
if (strcmp(platform,'rusty')); str_home = 'mnt/home'; end;
%%%%%%%%;
run(sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev/dir_m/setup_0',str_home));
flag_verbose = 1;
flag_disp = 1+flag_verbose; nf=0;
flag_replot = 0;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Assume path is set using dir_lakcluster_c/dir_m/setup_0.m. ;')); end;

if (flag_verbose); disp(sprintf(' %% Comparing Up05 with Ap05 data. ;')); end;
if (flag_verbose); disp(sprintf(' %% ncontinent = 0. ;')); end;

dir_code = sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev',str_home);
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman',str_home);
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
dataset.dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman/dir_%s',str_home,dataset.str_dataset);
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
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% We divide Up05 data into three continents. ;')); end;
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman',str_home);
mr_Up05_p01_continent_ = textread(sprintf('%s/dir_Up05/mr_Up05_p01_continent_.txt',dir_trunk));
n_continent = max(mr_Up05_p01_continent_)+1;
mr_Up05_p01_continent_pc__ = zeros(n_patient_Up05,n_continent);
for ncontinent=0:n_continent-1;
mr_Up05_p01_continent_pc__(:,1+ncontinent) = (mr_Up05_p01_continent_==ncontinent);
end;%for ncontinent=0:n_continent-1;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% We divide Ap05 data into three continents. ;')); end;
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman',str_home);
mr_Ap05_p01_continent_ = textread(sprintf('%s/dir_Ap05/mr_Ap05_from_0UKB_cap_Ap05_p01_continent_.txt',dir_trunk));
n_continent = max(mr_Ap05_p01_continent_)+1;
mr_Ap05_p01_continent_pc__ = zeros(n_patient_Ap05,n_continent);
for ncontinent=0:n_continent-1;
mr_Ap05_p01_continent_pc__(:,1+ncontinent) = (mr_Ap05_p01_continent_==ncontinent);
end;%for ncontinent=0:n_continent-1;
%%%%%%%%;

%%%%%%%%;
% Now we specify the biclustering run to use. ;
% To begin with, we limit ourselves to the training dataset Up05, and ncontinent=0. ;
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
assert(numel(mr_A_ori_trn_)==numel(mr_Up05_p01_continent_));
mr_A_alt_trn_ = mr_A_ori_trn_.*mr_Up05_p01_continent_pc__(:,1+ncontinent);
mr_Z_alt_trn_ = mr_Z_ori_trn_.*mr_Up05_p01_continent_pc__(:,1+ncontinent);
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
% to exclude those elements of Up05 which are not in the allele_cap ;
% defined above (i.e., the allele-combination-intersection of Up05 and Ap05). ;
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
pca_mc_A_trn = mc_A_rem_trn_; pca_mc_A_trn(1+index_Up05_setminus_cap_) = 0;
%%%%%%%%;
% Note that pca_mc_A_trn is 1 only for the allele-combinations in the bicluster ;
% that are also in the allele_cap. ;
%%%%%%%%;
pca_str_infix_trn=sprintf('trnUp05_tstAp05_ni%d_p01',tmp_ij_nlpR-1);
tmp_parameter_trn.flag_force_create = 0;
parameter_D_trnUp05_tstAp05_nix_p01 = tmp_parameter_trn;
parameter_D_trnUp05_tstAp05_nix_p01.str_A_p = str_Up05_A_p_p01;
parameter_D_trnUp05_tstAp05_nix_p01.str_name_s0000 = 'pca_D_trnUp05_tstAp05_nix_p01';
[ ...
 parameter_D_trnUp05_tstAp05_nix_p01 ...
,AZnV_D_trnUp05_tstAp05_nix_p01_ ...
,AnV_D_trnUp05_tstAp05_nix_p01_ ...
,ZnV_D_trnUp05_tstAp05_nix_p01_ ...
,V_D_trnUp05_tstAp05_nix_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_D_trnUp05_tstAp05_nix_p01 ...
,pca_rank ...
,pca_mr_A_trn_ ...
,pca_mr_Z_trn_ ...
,pca_mc_A_trn ...
,pca_str_infix_trn ...
);
%%%%;
mr_dvx_trnUp05_tstAp05_nix_ = 0.0*mr_A_ori_trn_;
mr_dvx_trnUp05_tstAp05_nix_(1+efind(mr_Z_alt_trn_)) = 1; %<-- color case-patients. ;
mr_dvx_trnUp05_tstAp05_nix_(1+efind(mr_A_alt_trn_)) = 2; %<-- color ctrl-patients. ;
mr_dvx_trnUp05_tstAp05_nix_(1+efind(mr_A_rem_trn_)) = 3; %<-- color bicluster-patients (a subset of cases). ;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_sml = 12;
markersize_big = 16;
hold on;
tmp_index_ = efind( (mr_dvx_trnUp05_tstAp05_nix_==1) | (mr_dvx_trnUp05_tstAp05_nix_==2) );
scatter(AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,1+0),AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,1+1),markersize_sml,mr_dvx_trnUp05_tstAp05_nix_(1+tmp_index_),'filled','MarkerEdgeColor','k');
tmp_index_ = efind( (mr_dvx_trnUp05_tstAp05_nix_==3) );
scatter(AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,1+0),AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,1+1),markersize_big,mr_dvx_trnUp05_tstAp05_nix_(1+tmp_index_),'filled','MarkerEdgeColor',0.85*[1,1,1]);
axis equal; grid on;
title('AZnV_D_trnUp05_tstAp05_nix_p01_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;
%%%%;

%%%%%%%%;
% save down first bicluster. ;
%%%%%%%%;
tmp_fam_iid_ = dataset_{1+ndataset_Up05}.fam_iid_(1+index_r_rem_trn_);
tmp_fname = sprintf('%s/%s_ni%d_fam_iid_.txt',parameter_trn.dir_out_s0000,parameter_trn.str_prefix,tmp_ij_nlpR-1);
if ~exist(tmp_fname,'file'); fp = fopen(tmp_fname,'w'); for nl=0:numel(tmp_fam_iid_)-1; fprintf(fp,'%s\n',tmp_fam_iid_{1+nl}); end; fclose(fp); end;
tmp_bim_vid_ = dataset_{1+ndataset_Up05}.bim_vid_(1+index_c_rem_trn_);
tmp_fname = sprintf('%s/%s_ni%d_bim_vid_.txt',parameter_trn.dir_out_s0000,parameter_trn.str_prefix,tmp_ij_nlpR-1);
if ~exist(tmp_fname,'file'); fp = fopen(tmp_fname,'w'); for nl=0:numel(tmp_bim_vid_)-1; fprintf(fp,'%s\n',tmp_bim_vid_{1+nl}); end; fclose(fp); end;
tmp_bim_name_ = dataset_{1+ndataset_Up05}.bim_name_(1+index_c_rem_trn_);
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
dir_trunk_Up05 = sprintf('%s/dir_Up05',dir_trunk);
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
parameter_0b.dir_trunk = dir_trunk_Up05;
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
assert(numel(mr_A_ori_0b_)==numel(mr_Up05_p01_continent_));
mr_A_alt_0b_ = mr_A_ori_0b_.*mr_Up05_p01_continent_pc__(:,1+ncontinent);
mr_A_alt_0b_(1+index_r_rem_trn_) = 0; %<-- Here we explicitly exclude the case-rows in the first bicluster. ;
mr_Z_alt_0b_ = mr_Z_ori_0b_.*mr_Up05_p01_continent_pc__(:,1+ncontinent);
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
% Now we calculate many _D_ principal-components across the entire Up05 data-set. ;
% The '_D_' refers to the situation where we calculate standard pcs using the cases ;
% (i.e., maximizing the variance across all the case-patients). ;
% First we set up the basic parameters. ;
%%%%%%%%;
str_lak_vs_dex = 'dex';
str_prefix = 'test2mds_maf01';
gamma = 0.05;
n_mds_0in = 2; n_mds_repl = 1; ij_mds_use_ = [1:2];
parameter_Up05.dir_code = dir_code;
parameter_Up05.dir_trunk = dir_trunk_Up05;
parameter_Up05.str_lak_vs_dex = str_lak_vs_dex;
parameter_Up05.str_prefix = str_prefix;
parameter_Up05.gamma = gamma;
parameter_Up05.n_mds = n_mds_0in;
parameter_Up05.n_mds_repl = n_mds_repl;
parameter_Up05.ij_mds_use_ = ij_mds_use_;
parameter_Up05.flag_force_create = 0; %<-- reload previous run. ;
parameter_Up05.flag_verbose = max(0,flag_verbose-1);
parameter_Up05.maf_lo_threshold = maf_lo_threshold;
parameter_Up05.n_shuffle = 128;
parameter_Up05.dir_0in = sprintf('%s/dir_%s',parameter_Up05.dir_trunk,parameter_Up05.str_prefix);
parameter_Up05.nshuffle = 0;
parameter_Up05 = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_Up05);
%%%%%%%%;
fname_mr_A_full = sprintf('%s/%s_mr_A_full.b16',parameter_Up05.dir_0in,parameter_Up05.str_prefix);
fname_mr_Z_full = sprintf('%s/%s_mr_Z_full.b16',parameter_Up05.dir_0in,parameter_Up05.str_prefix);
fname_mc_A = sprintf('%s/%s_mc_A.b16',parameter_Up05.dir_0in,parameter_Up05.str_prefix);
mr_A_ori_Up05_ = binary_uncompress(fname_mr_A_full)>0; %<-- thresholded to avoid negative values. ;
mr_Z_ori_Up05_ = binary_uncompress(fname_mr_Z_full)>0; %<-- thresholded to avoid negative values. ;
mc_A_ori_Up05_ = binary_uncompress(fname_mc_A)>0; %<-- thresholded to avoid negative values. ;
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
pca_str_infix_trn=sprintf('Up05_p01');
parameter_D_Up05_p01 = parameter_Up05;
parameter_D_Up05_p01.flag_force_create = 0;
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
,pca_mr_A_trn_ ...
,pca_mr_Z_trn_ ...
,pca_mc_A_trn ...
,pca_str_infix_trn ...
);
%%%%%%%%;
mx__ = load_mx__from_parameter_ver0(parameter_D_Up05_p01);
mr_dvx_ = 0.0*mx__.mr_A_full_;
mr_dvx_(1+efind(mx__.mr_A_full_))=2;
mr_dvx_(1+efind(mx__.mr_Z_full_))=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;fig80s;
%%%%;
subplot(1,2,1);
scatter(AZnV_D_Up05_p01_(:,1),AZnV_D_Up05_p01_(:,2),8,mr_dvx_,'filled');
xlabel('pc0'); ylabel('pc1'); title('case magenta, ctrl cyan');
axis equal; axis vis3d;
%%%%;
subplot(1,2,2);
scatter(AZnV_D_Up05_p01_(:,1),AZnV_D_Up05_p01_(:,2),8,mr_Up05_p01_continent_,'filled');
xlabel('pc0'); ylabel('pc1'); title('continent 0,1,2');
axis equal; axis vis3d;
%%%%;
fname_fig = sprintf('/%s/rangan/dir_bcc/dir_jelman/dir_Up05/Up05_D_p01_FIGA',str_home);
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if flag_disp;
%%%%%%%%;
end;%if flag_calc;

flag_calc = 1;
%%%%%%%%;
% Now calculate dominant _DandX_ principal-components across the entire Up05 data-set. ;
% The '_DandX_' refers to the situation where we calculate standard pcs using all the patients ;
% (i.e., maximizing the variance across all the case- and ctrl-patients provided). ;
% Note that we modify the parameters used above. ;
%%%%%%%%;
if flag_calc;
pca_rank = 2;
pca_mr_A_Up05_ = { 1*mr_A_ori_Up05_ + 1*mr_Z_ori_Up05_ };
pca_mr_Z_Up05_ = { 0*mr_A_ori_Up05_ + 0*mr_Z_ori_Up05_ };
pca_mc_A_Up05 = mc_A_ori_Up05_;
pca_str_infix_Up05=sprintf('Up05_DandX_p01');
parameter_DandX_Up05_p01 = parameter_Up05;
parameter_DandX_Up05_p01.flag_force_create = 0;
parameter_DandX_Up05_p01.str_A_p = str_Up05_A_p_p01;
parameter_DandX_Up05_p01.str_name_s0000 = 'pca_DandX_Up05_p01';
parameter_DandX_Up05_p01.slurm_memdecl = 48;
[ ...
 parameter_DandX_Up05_p01 ...
,AZnV_DandX_Up05_p01_ ...
,AnV_DandX_Up05_p01_ ...
,ZnV_DandX_Up05_p01_ ...
,V_DandX_Up05_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_DandX_Up05_p01 ...
,pca_rank ...
,pca_mr_A_Up05_ ...
,pca_mr_Z_Up05_ ...
,pca_mc_A_Up05 ...
,pca_str_infix_Up05 ...
);
%%%%%%%%;
mx__ = load_mx__from_parameter_ver0(parameter_DandX_Up05_p01);
mr_dvx_ = 0.0*mx__.mr_A_full_;
mr_dvx_(1+efind(mx__.mr_A_full_))=2;
mr_dvx_(1+efind(mx__.mr_Z_full_))=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;fig80s;
%%%%;
subplot(1,2,1);
%scatter3(AZnV_DandX_Up05_p01_(:,1),AZnV_DandX_Up05_p01_(:,2),AZnV_DandX_Up05_p01_(:,3),8,mr_dvx_,'filled');
scatter(AZnV_DandX_Up05_p01_(:,1),AZnV_DandX_Up05_p01_(:,2),8,mr_dvx_,'filled');
xlabel('pc0'); ylabel('pc1'); title('case magenta, ctrl cyan');
axis equal; axis vis3d;
%%%%;
subplot(1,2,2);
scatter(AZnV_DandX_Up05_p01_(:,1),AZnV_DandX_Up05_p01_(:,2),8,mr_Up05_p01_continent_,'filled');
xlabel('pc0'); ylabel('pc1'); title('continent 0,1,2');
axis equal; axis vis3d;
%%%%;
fname_fig = sprintf('/%s/rangan/dir_bcc/dir_jelman/dir_Up05/Up05_DandX_p01_FIGA',str_home);
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if flag_disp;
%%%%%%%%;
end;%if flag_calc;

flag_calc = 0;
if flag_calc;
%%%%%%%%;
% Now calculate many _DvX_ principal-components across the entire Up05 data-set. ;
% The '_DvX_' refers to the situation where we calculate 'category-corrected' pcs. ;
% These are designed to maximize the variance between cases and ctrls. ;
% Note that we modify the parameters used above. ;
%%%%%%%%;
pca_rank = 2;
pca_mr_A_trn_ = { 1*mr_A_ori_trn_ };
pca_mr_Z_trn_ = { 1*mr_Z_ori_trn_ };
pca_mc_A_trn = mc_A_ori_trn_;
pca_str_infix_trn=sprintf('Up05_p01');
parameter_DvX_Up05_p01 = parameter_Up05;
parameter_DvX_Up05_p01.flag_force_create = 0;
parameter_DvX_Up05_p01.str_A_p = str_Up05_A_p_p01;
parameter_DvX_Up05_p01.str_name_s0000 = 'pca_DvX_Up05_p01';
[ ...
 parameter_DvX_Up05_p01 ...
,AZnV_DvX_Up05_p01_ ...
,AnV_DvX_Up05_p01_ ...
,ZnV_DvX_Up05_p01_ ...
,V_DvX_Up05_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_DvX_from_mx_ver16( ...
 parameter_DvX_Up05_p01 ...
,pca_rank ...
,pca_mr_A_trn_ ...
,pca_mr_Z_trn_ ...
,pca_mc_A_trn ...
,pca_str_infix_trn ...
);
%%%%%%%%;
mx__ = load_mx__from_parameter_ver0(parameter_DvX_Up05_p01);
mr_dvx_ = 0.0*mx__.mr_A_full_;
mr_dvx_(1+efind(mx__.mr_A_full_))=2;
mr_dvx_(1+efind(mx__.mr_Z_full_))=1;
scatter(AZnV_DvX_Up05_p01_(:,1),AZnV_DvX_Up05_p01_(:,2),8,mr_dvx_);
%%%%%%%%;
end;%if flag_calc;

%%%%%%%%;
% Here we can take a quick peek at the histogram of cases (magenta) as ;
% projected onto the dominant principal-component of the bicluster. ;
% The ctrls are shown in blue. ;
% Note that this projection is limited to the training-set (Up05). ;
%%%%%%%%;
x_ = transpose(linspace(-60,+80,128)); 
xlim_ = [min(x_),max(x_)];
tmp_h_Z_ = hist(AZnV_D_trnUp05_tstAp05_nix_p01_(1+efind(mr_Z_alt_trn_),1+0),x_);
tmp_h_Z_ = tmp_h_Z_/sum(tmp_h_Z_);
tmp_h_A_ = hist(AZnV_D_trnUp05_tstAp05_nix_p01_(1+efind(mr_A_alt_trn_),1+0),x_);
tmp_h_A_ = tmp_h_A_/sum(tmp_h_A_);
figure(1+nf);nf=nf+1;figsml;
hold on;
stairs(x_,tmp_h_Z_,'b-');
stairs(x_,tmp_h_A_,'m-');
hold off;
xlim(xlim_);
%%%%%%%%;

%%%%%%%%;
% Now, for the training dataset Up05, ;
% We define the 'projected value' of each patient as the ;
% position of that patient after projection onto the ;
% dominant pc of the first bicluster ;
% (defined using only the allele-combinations in the bicluster which lie within allele_cap). ;
% Then we calculate the auc between cases and ctrls ;
% for those patients above a certain threshold of this projected value. ;
% We store this in: ;
% D_trnUp05_tstAp05_nix_p01_auc_t_, auc values. ;
% D_trnUp05_tstAp05_nix_p01_logp_auc_t_, log p-values for the auc. ;
% D_trnUp05_tstAp05_nix_p01_thr_t_, threshold values. ;
%%%%%%%%
tmp_Z_ = AZnV_D_trnUp05_tstAp05_nix_p01_(1+efind(mr_Z_alt_trn_),1+0);
tmp_A_ = AZnV_D_trnUp05_tstAp05_nix_p01_(1+efind(mr_A_alt_trn_),1+0);
tmp_threshold_ = sort([tmp_Z_;tmp_A_],'ascend'); tmp_n_threshold = numel(tmp_threshold_);
D_trnUp05_tstAp05_nix_p01_auc_t_ = zeros(tmp_n_threshold,1);
D_trnUp05_tstAp05_nix_p01_logp_auc_t_ = zeros(tmp_n_threshold,1);
for tmp_nthreshold=0:tmp_n_threshold-1;
tmp_threshold = tmp_threshold_(1+tmp_nthreshold);
tmp_Z_thr_ = tmp_Z_(1+efind(tmp_Z_>=tmp_threshold));
tmp_A_thr_ = tmp_A_(1+efind(tmp_A_>=tmp_threshold));
tmp_auc = auc_0(tmp_Z_thr_,tmp_A_thr_);
tmp_n = min(numel(tmp_Z_thr_),numel(tmp_A_thr_));
tmp_logp_auc = logp_auc_0(tmp_auc,tmp_n);
if (flag_verbose>1); disp(sprintf(' %% thr %0.2f: tmp_auc: %0.3f <-- nlp %0.2f <-- p %0.6f',tmp_threshold,tmp_auc,-tmp_logp_auc,exp(tmp_logp_auc))); end;
D_trnUp05_tstAp05_nix_p01_auc_t_(1+tmp_nthreshold) = tmp_auc;
D_trnUp05_tstAp05_nix_p01_logp_auc_t_(1+tmp_nthreshold) = tmp_logp_auc;
end;%for tmp_nthreshold=0:tmp_n_threshold-1;
D_trnUp05_tstAp05_nix_p01_thr_t_ = tmp_threshold_;

%%%%%%%%;
% Now project Ap05 data onto the principal-components defined using the bicluster. ;
% Note that we have to use those principal-components defined using the intersection ;
% of the bicluster with the allele_cap (see pca_str_infix_trn on the next line). ;
%%%%%%%%;
pca_str_infix_trn=sprintf('trnUp05_tstAp05_ni%d_p01',tmp_ij_nlpR-1); %<-- this defines the pca file to use. ;
parameter_D_trnUp05_tstAp05_nix_p01.str_V_ = sprintf('%s/dir_pca/dir_pca_mda/pca_D_%s_k%d_B44_V_.mda',parameter_D_trnUp05_tstAp05_nix_p01.dir_out_s0000,pca_str_infix_trn,pca_rank);
if ~exist(parameter_D_trnUp05_tstAp05_nix_p01.str_V_,'file');
disp(sprintf(' %% Warning, %s not found',parameter_D_trnUp05_tstAp05_nix_p01.str_V_));
end;%if ~exist(parameter_D_trnUp05_tstAp05_nix_p01.str_V_,'file');
D_trnUp05_tstAp05_nix_p01_V_ = mda_read_r8(parameter_D_trnUp05_tstAp05_nix_p01.str_V_);
if (flag_verbose); disp(sprintf(' %% size of intersection: %d',numel(index_Up05_from_cap_))); end;
if (flag_verbose); disp(sprintf(' %% Support of V_: %d',sum(abs(D_trnUp05_tstAp05_nix_p01_V_(:,1))>0))); end;
if (flag_verbose); disp(sprintf(' %% intersection of supp(V_) with cap: %d',numel(intersect(index_Up05_from_cap_,efind(abs(D_trnUp05_tstAp05_nix_p01_V_(:,1))>0))))); end;
D_Ap05_from_trnUp05_tstAp05_nix_p01_V_ = ij_Ap05_from_Up05__*D_trnUp05_tstAp05_nix_p01_V_;
if (flag_verbose); disp(sprintf(' %% fnorm(D_Ap05_from_trnUp05_tstAp05_nix_p01_V_(ij_Ap05_from_cap_,:) - D_trnUp05_tstAp05_nix_p01_V_(ij_Up05_from_cap_,:)): %0.16f',fnorm(D_Ap05_from_trnUp05_tstAp05_nix_p01_V_(ij_Ap05_from_cap_,:) - D_trnUp05_tstAp05_nix_p01_V_(ij_Up05_from_cap_,:)))); end;

%%%%%%%%;
% Now project Ap05 data onto D_Ap05_from_trnUp05_tstAp05_nix_p01_V_. ;
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
assert(numel(mr_A_ori_tst_)==numel(mr_Ap05_p01_continent_));
mr_A_alt_tst_ = mr_A_ori_tst_.*mr_Ap05_p01_continent_pc__(:,1+ncontinent);
mr_Z_alt_tst_ = mr_Z_ori_tst_.*mr_Ap05_p01_continent_pc__(:,1+ncontinent);
mc_A_alt_tst_ = mc_A_ori_tst_;
pca_mr_A_tst_ = { mr_A_alt_tst_ };
pca_mr_Z_tst_ = { mr_Z_alt_tst_ };
pca_mc_A_tst = mc_A_alt_tst_; pca_mc_A_tst(1+index_Ap05_setminus_cap_) = 0;
pca_str_infix_tst='D_Ap05_from_trnUp05_tstAp05_nix_p01';
parameter_D_Ap05_from_trnUp05_tstAp05_nix_p01 = tmp_parameter_tst;
parameter_D_Ap05_from_trnUp05_tstAp05_nix_p01.str_A_p = str_Ap05_A_p_p01;
parameter_D_Ap05_from_trnUp05_tstAp05_nix_p01.str_name_s0000 = ...
  sprintf('pca_D_Ap05c%d_from_trnUp05c%d_tstAp05c%d_nix_p01',ncontinent,ncontinent,ncontinent);
[ ...
 parameter_D_Ap05_from_trnUp05_tstAp05_nix_p01 ...
,AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_ ...
,AnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_ ...
,ZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_ ...
,tmp_V_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_proj_from_mx_ver16( ...
 parameter_D_Ap05_from_trnUp05_tstAp05_nix_p01 ...
,pca_rank ...
,D_Ap05_from_trnUp05_tstAp05_nix_p01_V_ ...
,pca_mr_A_tst_ ...
,pca_mr_Z_tst_ ...
,pca_mc_A_tst ...
,pca_str_infix_tst ...
);
%%%%;
mr_dvx_Ap05_from_trnUp05_tstAp05_nix_ = 0.0*mr_A_ori_tst_;
mr_dvx_Ap05_from_trnUp05_tstAp05_nix_(1+efind(mr_Z_alt_tst_)) = 1;
mr_dvx_Ap05_from_trnUp05_tstAp05_nix_(1+efind(mr_A_alt_tst_)) = 2;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
markersize_use = 12;
tmp_index_ = efind(mr_dvx_Ap05_from_trnUp05_tstAp05_nix_);
scatter(AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,1+0),AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,1+1),markersize_use,mr_dvx_Ap05_from_trnUp05_tstAp05_nix_(1+tmp_index_),'filled','MarkerEdgeColor','k');
axis equal; grid on;
title('AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
%%%%;

%%%%%%%%;
% Here we can take a quick peek at the histogram of cases (magenta) as ;
% projected onto the dominant principal-component of the bicluster. ;
% The ctrls are shown in blue. ;
% Note that this projection is limited to the testing-set (Ap05). ;
%%%%%%%%;
x_ = transpose(linspace(-60,+80,128)); 
xlim_ = [min(x_),max(x_)];
tmp_h_Z_ = hist(AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+efind(mr_Z_alt_tst_),1+0),x_);
tmp_h_Z_ = tmp_h_Z_/sum(tmp_h_Z_);
tmp_h_A_ = hist(AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+efind(mr_A_alt_tst_),1+0),x_);
tmp_h_A_ = tmp_h_A_/sum(tmp_h_A_);
figure(1+nf);nf=nf+1;figsml;
hold on;
stairs(x_,tmp_h_Z_,'b-');
stairs(x_,tmp_h_A_,'m-');
hold off;
xlim(xlim_);

%%%%%%%%;
% Now, for the training dataset Ap05, ;
% We define the 'projected value' of each patient as the ;
% position of that patient after projection onto the ;
% dominant pc of the first bicluster ;
% (defined using only the allele-combinations in the bicluster which lie within allele_cap). ;
% Then we calculate the auc between cases and ctrls ;
% for those patients above a certain threshold of this projected value. ;
% We store this in: ;
% D_Ap05_from_trnUp05_tstAp05_nix_p01_auc_t_, auc values. ;
% D_Ap05_from_trnUp05_tstAp05_nix_p01_logp_auc_t_, log p-values for the auc. ;
% D_Ap05_from_trnUp05_tstAp05_nix_p01_thr_t_, threshold values. ;
%%%%%%%%
tmp_Z_ = AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+efind(mr_Z_alt_tst_),1+0);
tmp_A_ = AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+efind(mr_A_alt_tst_),1+0);
tmp_threshold_ = sort([tmp_Z_;tmp_A_],'ascend'); tmp_n_threshold = numel(tmp_threshold_);
D_Ap05_from_trnUp05_tstAp05_nix_p01_auc_t_ = zeros(tmp_n_threshold,1);
D_Ap05_from_trnUp05_tstAp05_nix_p01_logp_auc_t_ = zeros(tmp_n_threshold,1);
for tmp_nthreshold=0:tmp_n_threshold-1;
tmp_threshold = tmp_threshold_(1+tmp_nthreshold);
tmp_Z_thr_ = tmp_Z_(1+efind(tmp_Z_>=tmp_threshold));
tmp_A_thr_ = tmp_A_(1+efind(tmp_A_>=tmp_threshold));
tmp_auc = auc_0(tmp_Z_thr_,tmp_A_thr_);
tmp_n = min(numel(tmp_Z_thr_),numel(tmp_A_thr_));
tmp_logp_auc = logp_auc_0(tmp_auc,tmp_n);
if (flag_verbose>1); disp(sprintf(' %% thr %0.2f: tmp_auc: %0.3f <-- nlp %0.2f <-- p %0.6f',tmp_threshold,tmp_auc,-tmp_logp_auc,exp(tmp_logp_auc))); end;
D_Ap05_from_trnUp05_tstAp05_nix_p01_auc_t_(1+tmp_nthreshold) = tmp_auc;
D_Ap05_from_trnUp05_tstAp05_nix_p01_logp_auc_t_(1+tmp_nthreshold) = tmp_logp_auc;
end;%for tmp_nthreshold=0:tmp_n_threshold-1;
D_Ap05_from_trnUp05_tstAp05_nix_p01_thr_t_ = tmp_threshold_;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Now we plot the log-p-values for the auc of cases-vs-ctrls ;
% for each of the thresholds used above.
% This is after we project onto the dominant pc of the bicluster, ;
% defined using only the shared allele-combinations. ;
% For this plot the training-data is shown in red, ;
% while the testing-data is shown in green. ;
% Note that there is a large range of thresholds for which ;
% the testing-data is very significant. ;
% (significance value of 0.05 is marked with dashed line). ;
% (everything above dashed line is significant). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
ymax = 10;
tmp_dlim_0_ = [min(D_trnUp05_tstAp05_nix_p01_thr_t_);max(D_trnUp05_tstAp05_nix_p01_thr_t_)];
tmp_dlim_1_ = [min(D_Ap05_from_trnUp05_tstAp05_nix_p01_thr_t_);max(D_Ap05_from_trnUp05_tstAp05_nix_p01_thr_t_)];
tmp_dlim_ = [min([tmp_dlim_0_;tmp_dlim_1_]),max([tmp_dlim_0_;tmp_dlim_1_])];
stairs(D_trnUp05_tstAp05_nix_p01_thr_t_,min(10,-D_trnUp05_tstAp05_nix_p01_logp_auc_t_),'r-','LineWidth',2);
stairs(D_Ap05_from_trnUp05_tstAp05_nix_p01_thr_t_,min(10,-D_Ap05_from_trnUp05_tstAp05_nix_p01_logp_auc_t_),'g-','LineWidth',2);
%stairs(D_Ap05_from_trnUp05_tstAp05_nix_p01_thr_t_,min(10,10*D_Ap05_from_trnUp05_tstAp05_nix_p01_auc_t_),'-','Color',0.85*[1,1,1],'LineWidth',2);
plot(tmp_dlim_,-log(0.05)*ones(2,1),'k:');
hold off;
xlim(tmp_dlim_);
ylim([0,10.1]); grid on;
xlabel('pc1 threshold');
ylabel('-log(p(AUC))','Interpreter','none');
legend({'train','test'},'Location','NorthWest');
title('AZnV_D_trnUp05_tstAp05_nix_p01_','Interpreter','none');
end;%if flag_disp;

