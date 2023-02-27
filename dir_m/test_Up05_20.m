clear;
%%%%%%%%;
% Here we set the matlab path and define str_home for naming. ;
%%%%%%%%;
% platform = 'rusty';
% if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
% if (strcmp(platform,'access1')); str_home = 'data'; end;
% if (strcmp(platform,'OptiPlex')); str_home = 'home'; end;
% if (strcmp(platform,'eval1')); str_home = 'home'; end;
% if (strcmp(platform,'rusty')); str_home = 'mnt/home'; end;
% %%%%%%%%;
run('/home/jelman/Github/lakcluster_c/dir_m/setup_0'); %<-- set up the paths. ;
flag_verbose = 1;
flag_disp = 1+flag_verbose; nf=0;
flag_replot = 0;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Assume path is set using dir_lakcluster_c/dir_m/setup_0.m. ;')); end;

if (flag_verbose); disp(sprintf(' %% Comparing Up99 with Up05 data. ;')); end;
memory_GB = 128; %<-- maybe this should be increased? ;
if (flag_verbose); disp(sprintf(' %% trying with memory_GB %d ;',memory_GB)); end;

dir_code = '/home/jelman/Github/lakcluster_c';
dir_trunk = '/home/jelman/Projects/AD_Biclustering/data/UKB/ukb_pca_p05-p1';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;


dataset_Up05 = struct('type','dataset');
dataset_Up05.str_dataset = 'Up05';
dataset_Up05.dir_trunk = sprintf('%s/dir_%s',dir_trunk,dataset_Up05.str_dataset);
dataset_Up05.str_prefix = 'test2mds_maf01';
dataset_Up05.dir_0in = sprintf('%s/dir_%s',dataset_Up05.dir_trunk,dataset_Up05.str_prefix);
%%%%%%%%;
dataset_Up05.fname_bimext = sprintf('%s/%s_bim.ext',dataset_Up05.dir_0in,dataset_Up05.str_prefix);
if ~exist(dataset_Up05.fname_bimext,'file');
dataset_Up05.fname_bimext = sprintf('%s/%s_bim.ext',dataset_Up05.dir_0in,dataset_Up05.str_prefix);
end;%if ~exist(dataset_Up05.fname_bimext,'file');
if ~exist(dataset_Up05.fname_bimext,'file');
dataset_Up05.fname_bimext = sprintf('%s/test2mds_maf01_bim.ext',dataset_Up05.dir_0in);
end;%if ~exist(dataset_Up05.fname_bimext,'file');
if ~exist(dataset_Up05.fname_bimext,'file');
dataset_Up05.fname_bimext = sprintf('%s/%s_bim.ext',dataset_Up05.dir_0in,dataset_Up05.str_dataset);
end;%if ~exist(dataset_Up05.fname_bimext,'file');
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% load_bimext_ver1 ...')); end;
[ ...
,dataset_Up05.n_snp ...
,dataset_Up05.bim_khr_ ...
,dataset_Up05.bim_vid_ ...
,dataset_Up05.bim_gdi_ ...
,dataset_Up05.bim_pdi_ ...
,dataset_Up05.bim_al1_ ...
,dataset_Up05.bim_al2_ ...
,dataset_Up05.bim_alt_ ...
,dataset_Up05.bim_ent_ ...
,dataset_Up05.bim_frq_ ...
,dataset_Up05.bim_mss_ ...
,dataset_Up05.bim_maf_ ...
,dataset_Up05.bim_name_ ...
,~ ...
] = ...
load_bimext_ver1( ...
dataset_Up05.fname_bimext ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% load_bimext_ver1: %0.6fs',tmp_t)); end;
%%%%%%%%;
dataset_Up05.fname_famext = sprintf('%s/%s_fam.ext',dataset_Up05.dir_0in,dataset_Up05.str_prefix);
if ~exist(dataset_Up05.fname_famext,'file');
dataset_Up05.fname_famext = sprintf('%s/%s_fam.ext',dataset_Up05.dir_0in,dataset_Up05.str_prefix);
end;%if ~exist(dataset_Up05.fname_famext,'file');
if ~exist(dataset_Up05.fname_famext,'file');
dataset_Up05.fname_famext = sprintf('%s/test2mds_maf01_fam.ext',dataset_Up05.dir_0in);
end;%if ~exist(dataset_Up05.fname_famext,'file');
if ~exist(dataset_Up05.fname_famext,'file');
dataset_Up05.fname_famext = sprintf('%s/%s_fam.ext',dataset_Up05.dir_0in,dataset_Up05.str_dataset_Up05);
end;%if ~exist(dataset_Up05.fname_famext,'file');
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% load_famext_ver1 ...')); end;
[ ...
,dataset_Up05.n_patient ...
,dataset_Up05.fam_fid_ ...
,dataset_Up05.fam_iid_ ...
,dataset_Up05.fam_yid_ ...
,dataset_Up05.fam_xid_ ...
,dataset_Up05.fam_sex_ ...
,dataset_Up05.fam_dvx_ ...
,dataset_Up05.fam_dir_ ...
,dataset_Up05.fam_fidandiid_ ...
,~ ...
  ] = ...
load_famext_ver1( ...
 dataset_Up05.fname_famext ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% load_famext_ver1: %0.6fs',tmp_t)); end;
%%%%%%%%;

[~,n_row,n_col] = binary_getsize(sprintf('%s/%s_A_full_n.b16',dataset_Up05.dir_0in,dataset_Up05.str_prefix));
if (flag_verbose); disp(sprintf(' %% %% A_full_n_.b16 [%d,%d]',n_row,n_col)); end;
if (flag_verbose); disp(sprintf(' %% %% n_patient [%d]',dataset_Up05.n_patient)); end;
if (flag_verbose); disp(sprintf(' %% %% n_snp [%d]',dataset_Up05.n_snp)); end;

parameter_Up05 = struct('type','parameter');
parameter_Up05.slurm_memdecl = memory_GB; %<-- maybe this should be increased? ;
parameter_Up05.dir_0in = dataset_Up05.dir_0in;
parameter_Up05.str_prefix = dataset_Up05.str_prefix;
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% load_mx__from_parameter_ver0 ...')); end;
mx__ = load_mx__from_parameter_ver0(parameter_Up05);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% load_mx__from_parameter_ver0: %0.6fs',tmp_t)); end;
dataset_Up05.parameter = parameter_Up05;
dataset_Up05.mx__ = mx__;


dir_trunk_Up05 = dataset_Up05.dir_trunk;
n_patient_Up05 = dataset_Up05.n_patient;
n_snp_Up05 = dataset_Up05.n_snp;

pca_rank_use = 2; %<-- number of principal-components. ;
maf_lo_threshold = 0.01; %<-- minor-allele-frequency lower bound. ;

%%%%%%%%;
% Set Up05_A_p_. ;
%%%%%%%%;
mx_tst__ = dataset_Up05.mx__;
tmp_parameter_tst = dataset_Up05.parameter;
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


if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% We divide Up05 data into three continents. ;')); end;
mr_Up05_p01_continent_ = textread(sprintf('%s/dir_Up05/mr_Up05_p01_continent_.txt',dir_trunk));
n_continent = max(mr_Up05_p01_continent_)+1; %<-- continent index is zero-based. ;
mr_Up05_p01_continent_pc__ = zeros(n_patient_Up05,n_continent);
for ncontinent=0:n_continent-1;
mr_Up05_p01_continent_pc__(:,1+ncontinent) = (mr_Up05_p01_continent_==ncontinent);
end;%for ncontinent=0:n_continent-1;



flag_calc=1;
if flag_calc;
%%%%%%%%;
% Additionally, we can load another run (with more traces) ;
% and estimate the p-value associated with the original data. ;
% Note that this estimation (p_all_out) uses both the average and the maximum. ;
% This might not be what you want. ;
% For example, you might only want to consider the p-value as estimated from the 'maximum' ;
% of each trace. In this case you could use the p_top_out instead of p_all_out. ;
% Note also that the default setting is to ignore the bottom 5% and the top 5% of iterations. ;
% when calculating the maximum.
% The default is set this way because (for sc-RNA-seq data) these iterations are usually not informative. ;
%%%%%%%%;
str_lak_vs_dex = 'dex';
str_prefix = 'test2mds_maf01';
str_mr_0in = 'continent3';
str_mc_0in = '';
gamma = 0.05;
n_mds_0in = 2; n_mds_repl = 1; ij_mds_use_ = [1:2];
parameter_Up05 = struct('type','parameter');
parameter_Up05.slurm_memdecl = memory_GB; %<-- maybe this should be increased?
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
parameter_Up05.str_mr_0in = str_mr_0in;
parameter_Up05.n_shuffle = 500;
parameter_Up05.dir_0in = sprintf('%s/dir_%s',parameter_Up05.dir_trunk,parameter_Up05.str_prefix);
for nshuffle=1:parameter_Up05.n_shuffle;
parameter_Up05.nshuffle = nshuffle;
parameter_Up05 = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_Up05);
end;%for nshuffle=0:parameter_Up05.n_shuffle;
%%%%%%%%;
end;%if flag_calc;

% Load trace, get exclude mask for p-value calculation and Z-max
% Pval: exclude bicluster with <5% of case population
% Zmax: Exclude biclusters with <5% or >95% of case population
trace_Up05__ = load_trace__from_dir_ver0(parameter_Up05.dir_out_trace);
tmp_r_eli_ = trace_Up05__.r_rem_s0000_(1+0) - trace_Up05__.r_rem_s0000_ ;
ncases = trace_Up05__.r_rem_s0000_(1+0);
pval_mask = trace_Up05__.r_rem_s0000_/ncases > .05;
Z_max_mask = (trace_Up05__.r_rem_s0000_/ncases > .05) & (trace_Up05__.r_rem_s0000_/ncases < .95);

% Calculate p-value for trace
parameter_Up05.flag_verbose = 1; 
parameter_Up05.niteration_alo = 1;
parameter_Up05.niteration_ahi = max(trace_Up05__.niter_s0000_(pval_mask))+1;
[~,p_all_out,p_avg_out,p_top_out] = p_from_trace__ver0(parameter_Up05,trace_Up05__);


%%%%%%%%;
% Now run Z_imax_zerobased to pick out an internal iteration of interest. ;
% This particular function should be within /dir_lakcluster_c/dir_m_dependencies. ;
%%%%%%%%;
tmp_Z_ = trace_Up05__.ZR_s0000_(Z_max_mask); %<-- this is the z-score for the traces of the original data (i.e., nshuffle==0). ;
tmp_Z_ = trace_Up05__.ZR_is__(Z_max_mask,1+0); %<-- this is the same as the previous definition. ;
tmp_Z_min = -Inf; %<-- The lowest z-score to consider, putting in -Inf will consider all z-scores. ;
[tmp_Z_max,tmp_Z_max_index] = Z_imax_zerobased(flag_verbose,tmp_Z_,tmp_Z_min);
%%%%;

%%%%;
% Here we visualize the traces. ;
% We put a little red circle at the maximum value (picked out by Z_imax_zerobased). ;
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
linewidth_big = 3;
linewidth_sml = 1;
subplot(1,2,1);
hold on;
plot(trace_Up05__.niter_s0000_(pval_mask),trace_Up05__.ZR_is__(pval_mask,2:end),'b-','LineWidth',linewidth_sml,'Color', [0 0 0 0.1]);
plot(trace_Up05__.niter_s0000_(pval_mask),trace_Up05__.ZR_is__(pval_mask,1 + 0),'r-','LineWidth',linewidth_big);
plot(trace_Up05__.niter_s0000_(1+tmp_Z_max_index),trace_Up05__.ZR_is__(1+tmp_Z_max_index,1+0),'ro','MarkerFaceColor','r','MarkerSize',10);
hold off;
xlim([0,max(trace_Up05__.niter_s0000_(pval_mask))+1]); xlabel('iteration'); ylim([-10,+10]); ylabel('Z');
subplot(1,2,2);
hold on;
plot(tmp_r_eli_(pval_mask)/ncases,trace_Up05__.ZR_is__(pval_mask,2:end),'b-','LineWidth',linewidth_sml,'Color', [0 0 0 0.1]);
plot(tmp_r_eli_(pval_mask)/ncases,trace_Up05__.ZR_is__(pval_mask,1 + 0),'r-','LineWidth',linewidth_big);
plot(tmp_r_eli_(1+tmp_Z_max_index)/ncases,trace_Up05__.ZR_is__(1+tmp_Z_max_index,1+0),'ro','MarkerFaceColor','r','MarkerSize',10);
hold off;
xlim([min(tmp_r_eli_(pval_mask))/ncases,max(tmp_r_eli_(pval_mask))/ncases]); xlabel('proportion of cases eliminated'); ylim([-10,+10]); ylabel('Z');
%%%%;

%%%%;
% Here we visualize the scatter plots of cases and controls. ;
% First plot shows iteration 0, second shows bicluster. ;
%%%%;
fname_mr_A_full = sprintf('%s/%s_mr_A_full.b16',parameter_Up05.dir_0in,parameter_Up05.str_prefix);
fname_mr_Z_full = sprintf('%s/%s_mr_Z_full.b16',parameter_Up05.dir_0in,parameter_Up05.str_prefix);
fname_mc_A = sprintf('%s/%s_mc_A.b16',parameter_Up05.dir_0in,parameter_Up05.str_prefix);
mr_A_ori_Up05_ = binary_uncompress(fname_mr_A_full)>0; %<-- thresholded to avoid negative values. ;
mr_Z_ori_Up05_ = binary_uncompress(fname_mr_Z_full)>0; %<-- thresholded to avoid negative values. ;
mc_A_ori_Up05_ = binary_uncompress(fname_mc_A)>0; %<-- thresholded to avoid negative values. ;

mr_dvx_ = zeros(n_patient_Up05,1);
mr_dvx_(find(mr_A_ori_Up05_))=2;
mr_dvx_(find(mr_Z_ori_Up05_))=1;
index_rsub_ = efind(mr_dvx_);
parameter_Up05.str_out_xdrop_a_s0000 = sprintf('%s/out_xdrop_a.txt',parameter_Up05.dir_out_s0000);
[ ...
 xdrop_ ...
] = ...
load_out_xdrop_from_str_ver0( ...
parameter_Up05.str_out_xdrop_a_s0000 ...
);
index_rdrop_ = xdrop_.index_rdrop_;
index_cdrop_ = xdrop_.index_cdrop_;
if flag_disp;
figure(1+nf);nf=nf+1;clf;
figmed;fig80s;
end;%if flag_disp;
p_row = 1; p_col = 2; np=0;
markersize_sml = 4;
markersize_med = 8;
markersize_big = 16;
xdrop_ij_rkeep_ = xdrop_.ij_rkeep_;
xdrop_ij_ckeep_ = xdrop_.ij_ckeep_;
%%%%;
[A_p_c_,A_p_0_,AZ_rsum_] = load_A_p_c_from_dir_0(parameter_Up05.dir_out_s0000);
alpha_c_ = A_p_c_ - (1-A_p_c_);
D_c_ = sqrt(1./max(0.01,4.0*A_p_c_.*(1-A_p_c_)));
mx__ = load_mx__from_parameter_ver0(parameter_Up05);
%%%%%%%%;
ni=0;
[ ...
 AZnV_ni0_driver_pd__ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_ni_ver16( ...
 parameter_Up05 ...
,ni ...
,xdrop_ ...
,trace_Up05__ ...
,mx__ ...
);
if flag_disp;
subplot(p_row,p_col,1+np);np=np+1;
hold on
scatter(AZnV_ni0_driver_pd__(mr_dvx_(1+index_rsub_)==1,1),AZnV_ni0_driver_pd__(mr_dvx_(1+index_rsub_)==1,2),16,'m','filled','MarkerEdgeColor','k','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);
scatter(AZnV_ni0_driver_pd__(mr_dvx_(1+index_rsub_)==2,1),AZnV_ni0_driver_pd__(mr_dvx_(1+index_rsub_)==2,2),16,'c','filled','MarkerEdgeColor','k','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);
legend('Controls','Cases','Location','Best');
axisnotick; title('Full dataset'); xlabel('PC1'); ylabel('PC2');
hold off
end;%if flag_disp;
%%%%%%%%;
ni=tmp_Z_max_index;
[ ...
 AZnV_nix_driver_pd__ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_ni_ver16( ...
 parameter_Up05 ...
,ni ...
,xdrop_ ...
,trace_Up05__ ...
,mx__ ...
);
if flag_disp;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
scatter(AZnV_nix_driver_pd__(mr_dvx_(1+index_rsub_)==1,1),AZnV_nix_driver_pd__(mr_dvx_(1+index_rsub_)==1,2),16,'m','filled','MarkerEdgeColor','k','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);
scatter(AZnV_nix_driver_pd__(mr_dvx_(1+index_rsub_)==2,1),AZnV_nix_driver_pd__(mr_dvx_(1+index_rsub_)==2,2),16,'c','filled','MarkerEdgeColor','k','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);
legend('Controls','Cases','Location','Best');
axisnotick; title('Bicluster'); xlabel('PC1'); ylabel('PC2');
hold off;
end;%if flag_disp;
%%%%%%%%%;



flag_calc = 1;
if flag_calc;    
%%%%%%%%;
% Now load the out_xdrop_a.txt file. ;
%%%%%%%%;
xdrop_Up05_ = load_out_xdrop_from_str_ver0(parameter_Up05.str_out_xdrop_a_s0000);
%%%%%%%%;
% Note that 'index_rdrop_' contains the row indices (zero-based) in the order they were eliminated. ;
% The 'ij_rdrop_' contains the row indices (one-based) in the order they were eliminated. ;
% The 'index_rkeep_' and 'ij_rkeep_' are the reversed lists (respectively). ;
%%%%%%%%;
% Now extract a bicluster from iteration 16 (you can change this if you like). ;
%%%%%%%%;
tmp_niteration = tmp_Z_max_index; %<-- this is the iteration index (zero-based) at which you want to define the bicluster. ;
tmp_rkeep = trace_Up05__.r_rem_s0000_(1+tmp_niteration); %<-- number of rows to keep. ;
tmp_rdrop = trace_Up05__.r_rem_s0000_(1+0) - tmp_rkeep; %<-- number of rows to drop. ;
tmp_ckeep = trace_Up05__.c_rem_s0000_(1+tmp_niteration); %<-- number of cols to keep. ;
tmp_cdrop = trace_Up05__.c_rem_s0000_(1+0) - tmp_ckeep; %<-- number of cols to drop. ;
%%%%%%%%;
% We store the indices for this bicluster as an array with the following format: ;
% tmp_xdrop__ = [ row_index ,    -1     ] ;
%               [ row_index ,    -1     ] ;
%               [    ...    ,    -1     ] ;
%               [ row_index ,    -1     ] ;
%               [    -1     , col_index ] ;
%               [    -1     , col_index ] ;
%               [    -1     ,    ...    ] ;
%               [    -1     , col_index ] ;
% Where 'row_index' lists the various different (zero-based) row-indices of the bicluster, ;
% and 'col_index' lists the various different (zero-based) col-indices of the bicluster. ;
%%%%%%%%;
tmp_xdrop__ = [ xdrop_Up05_.index_rkeep_(1:tmp_rkeep) , -ones(tmp_rkeep,1) ; -ones(tmp_ckeep , 1) , xdrop_Up05_.index_ckeep_(1:tmp_ckeep) ];
tmp_fname_xdrop = sprintf('%s/%s_BC0_xdrop_ni%d.txt',parameter_Up05.dir_out_s0000,parameter_Up05.str_prefix,tmp_niteration);
if ~exist(tmp_fname_xdrop,'file');
disp(sprintf(' %% %s not found, creating',tmp_fname_xdrop));
tmp_fp = fopen(tmp_fname_xdrop,'w');
fprintf(tmp_fp,'%d %d\n',transpose(tmp_xdrop__));
fclose(tmp_fp);
end;%if ~exist(tmp_fname_xdrop,'file');
%%%%%%%%;
% Save out list of snp IDs corresponding to first bicluster
%%%%%%%%;
tmp_xdrop_vid_ = dataset_Up05.bim_vid_(xdrop_Up05_.index_ckeep_(1:tmp_ckeep));
tmp_xdrop_al1_ = dataset_Up05.bim_al1_(xdrop_Up05_.index_ckeep_(1:tmp_ckeep));
tmp_xdrop_al2_ = dataset_Up05.bim_al2_(xdrop_Up05_.index_ckeep_(1:tmp_ckeep));
tmp_xdrop_alt_ = dataset_Up05.bim_alt_(xdrop_Up05_.index_ckeep_(1:tmp_ckeep));
xdrop_vid_fname = sprintf('%s/%s_BC0_xdrop_ni%d_vid.txt',parameter_Up05.dir_out_s0000,parameter_Up05.str_prefix,tmp_niteration);
writecell([tmp_xdrop_vid_ cellstr(tmp_xdrop_al1_)  cellstr(tmp_xdrop_al2_) tmp_xdrop_alt_], xdrop_vid_fname,'Delimiter','\t');

%%%%%%%%;
% Now we rerun the biclustering once more, instructing the driver to 'scramble' this bicluster. ;
% Note that you can actually define the bicluster to scramble any sub-array you like ;
% by setting the row_index and col_index entries in the tmp_xdrop__ array to be whatever you want. ;
% Just remember that these are all zero-based. ;
%%%%%%%%;
parameter_Up05_scramble_BC0 = struct('type','parameter');
parameter_Up05_scramble_BC0.slurm_memdecl = memory_GB;
parameter_Up05_scramble_BC0.dir_code = dir_code;
parameter_Up05_scramble_BC0.dir_trunk = dir_trunk_Up05;
parameter_Up05_scramble_BC0.str_lak_vs_dex = str_lak_vs_dex;
parameter_Up05_scramble_BC0.str_prefix = str_prefix;
parameter_Up05_scramble_BC0.str_mr_0in = str_mr_0in;
parameter_Up05_scramble_BC0.str_mc_0in = str_mc_0in;
parameter_Up05_scramble_BC0.gamma = gamma;
parameter_Up05_scramble_BC0.n_mds = n_mds_0in;
parameter_Up05_scramble_BC0.n_mds_repl = n_mds_repl;
parameter_Up05_scramble_BC0.ij_mds_use_ = ij_mds_use_;
parameter_Up05_scramble_BC0.flag_force_create = 0; %<-- reload previous run. ;
parameter_Up05_scramble_BC0.flag_verbose = max(0,flag_verbose-1);
parameter_Up05_scramble_BC0.maf_lo_threshold = maf_lo_threshold;
parameter_Up05_scramble_BC0.n_shuffle = 64;
parameter_Up05_scramble_BC0.dir_0in = sprintf('%s/dir_%s',parameter_Up05_scramble_BC0.dir_trunk,parameter_Up05_scramble_BC0.str_prefix);
%%%%;
parameter_Up05_scramble_BC0.nshuffle = 0; %<-- this should be set to 0 when scrambling a previous bicluster. ;
parameter_Up05_scramble_BC0.n_scramble = 1; %<-- this is 0 by default, but is now 1 since we are scrambling 1 bicluster. ;
parameter_Up05_scramble_BC0.scramble_out_xdrop_ = {tmp_fname_xdrop}; %<-- this is empty by default, but is now a cell-array of length n_scramble. ;
parameter_Up05_scramble_BC0.scramble_rseed_ = [1]; %<-- this is empty by default, but is now an array of length n_scramble. ;
%%%%;
parameter_Up05_scramble_BC0 = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_Up05_scramble_BC0);
%%%%%%%%;
% Now we can visualize the results once more. ;
% This time the trace from the scrambled run (after scrambling BC0) is shown in green. ;
% Note that we compare the scrambled run to the shuffled samples from the original run. ;
% (see the definition of ZR_scramble_BC0_ below). ;
%%%%;
trace_Up05__ = load_trace__from_dir_ver0(parameter_Up05.dir_out_trace);
tmp_r_eli_ = trace_Up05__.r_rem_s0000_(1+0) - trace_Up05__.r_rem_s0000_ ;
trace_Up05_scramble_BC0__ = load_trace__from_dir_ver0(parameter_Up05_scramble_BC0.dir_out_trace);
ZR_scramble_BC0_ = (trace_Up05_scramble_BC0__.QR_s0000_ - trace_Up05__.QR_avg_i_)./max(1e-12,trace_Up05__.QR_std_i_);
figure(1+nf);nf=nf+1;clf;figmed;
linewidth_big = 3;
linewidth_sml = 1;
subplot(1,2,1);
hold on;
plot(trace_Up05__.niter_s0000_,trace_Up05__.ZR_is__(:,2:end),'b-','LineWidth',linewidth_sml,'Color', [0 0 0 0.3]);
plot(trace_Up05__.niter_s0000_,trace_Up05__.ZR_is__(:,1 + 0),'r-','LineWidth',linewidth_big);
plot(trace_Up05__.niter_s0000_,ZR_scramble_BC0_,'g-','LineWidth',linewidth_big);
hold off;
xlim([0,trace_Up05__.niter_is__(end,1+0)+1]); xlabel('niteration'); ylim([-10,+10]); ylabel('ZR');
subplot(1,2,2);
hold on;
plot(tmp_r_eli_,trace_Up05__.ZR_is__(:,2:end),'b-','LineWidth',linewidth_sml,'Color', [0 0 0 0.3]);
plot(tmp_r_eli_,trace_Up05__.ZR_is__(:,1 + 0),'r-','LineWidth',linewidth_big);
plot(tmp_r_eli_,ZR_scramble_BC0_,'g-','LineWidth',linewidth_big);
hold off;
xlim([min(tmp_r_eli_),max(tmp_r_eli_)]); xlabel('rows eliminated'); ylim([-10,+10]); ylabel('ZR');
%%%%;
%%%%%%%%;
end;%if flag_calc;



% Pval: exclude bicluster with <5% of case population
% Zmax: Exclude biclusters with <5% or >95% of case population
ncases_scramble_BC0 = trace_Up05_scramble_BC0__.r_rem_s0000_(1+0);
pval_mask_scramble_B0 = trace_Up05_scramble_BC0__.r_rem_s0000_/ncases_scramble_BC0 > .05;
Z_max_mask_scramble_B0 = (trace_Up05_scramble_BC0__.r_rem_s0000_/ncases_scramble_BC0 > .05) & (trace_Up05_scramble_BC0__.r_rem_s0000_/ncases_scramble_BC0 < .95);

% Calculate p-value for trace
parameter_Up05_scramble_BC0.flag_verbose = 1; 
parameter_Up05_scramble_BC0.niteration_alo = 1;
parameter_Up05_scramble_BC0.niteration_ahi = max(trace_Up05_scramble_BC0__.niter_s0000_(pval_mask_scramble_B0))+1;
[~,p_all_out_scramble_B0,p_avg_out_scramble_B0,p_top_out_scramble_B0] = p_from_trace__ver0(parameter_Up05_scramble_BC0,trace_Up05_scramble_BC0__);


%%%%%%%%;
% Now run Z_imax_zerobased to pick out an internal iteration of interest. ;
% This particular function should be within /dir_lakcluster_c/dir_m_dependencies. ;
%%%%%%%%;
tmp_Z_scramble_B0_ = trace_Up05_scramble_BC0__.ZR_s0000_(Z_max_mask_scramble_B0); %<-- this is the z-score for the traces of the original data (i.e., nshuffle==0). ;
tmp_Z_scramble_B0_ = trace_Up05_scramble_BC0__.ZR_is__(Z_max_mask_scramble_B0,1+0); %<-- this is the same as the previous definition. ;
tmp_Z_min_scramble_B0 = -Inf; %<-- The lowest z-score to consider, putting in -Inf will consider all z-scores. ;
[tmp_Z_max,tmp_Z_max_index] = Z_imax_zerobased(flag_verbose,tmp_Z_scramble_B0_,tmp_Z_min_scramble_B0);
%%%%;




flag_calc = 1;
%%%%%%%%;
% Now calculate dominant _DandX_ principal-components across the entire Up05 data-set. ;
%%%%%%%%;
if flag_calc;

pca_rank = 2;
p_threshold_ = 0.05:0.05:1.00; n_p_threshold = numel(p_threshold_);
%p_threshold_ = 0.60:0.005:0.65; n_p_threshold = numel(p_threshold_);
AZnV_DandX_Up05_p01_pnt___ = zeros(dataset_Up05.n_patient,pca_rank,n_p_threshold);
%%%%%%%%;
for np_threshold=0:n_p_threshold-1;
%%%%%%%%;
p_threshold = p_threshold_(1+np_threshold);
str_p_threshold = sprintf('%.2d',min(99,floor(100*p_threshold)));
pca_mr_A_Up05_ = { 1*mr_A_ori_Up05_ + 1*mr_Z_ori_Up05_ };
pca_mr_Z_Up05_ = { 0*mr_A_ori_Up05_ + 0*mr_Z_ori_Up05_ };
pca_mc_A_Up05 = mc_A_ori_Up05_;
tmp_mc_ = zeros(dataset_Up05.n_snp,1);
tmp_mc_(1+efind(Up05_bim_ADp_<=p_threshold))=1;
pca_mc_A_Up05 = tmp_mc_;
pca_str_infix_Up05=sprintf('Up05t%s_DandX_p01',str_p_threshold);
parameter_DandX_Up05_p01 = parameter_Up05;
parameter_DandX_Up05_p01.flag_force_create = 0;
parameter_DandX_Up05_p01.str_A_p = str_Up05_A_p_p01;
parameter_DandX_Up05_p01.str_name_s0000 = sprintf('pca_Up05t%s_DandX_p01',str_p_threshold);
tmp_t = tic(); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16 ...')); end;
[ ...
 parameter_DandX_Up05_p01 ...
,tmp_AZnV_DandX_Up05_p01_ ...
,tmp_AnV_DandX_Up05_p01_ ...
,tmp_ZnV_DandX_Up05_p01_ ...
,tmp_V_DandX_Up05_p01_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16( ...
 parameter_DandX_Up05_p01 ...
,pca_rank ...
,pca_mr_A_Up05_ ...
,pca_mr_Z_Up05_ ...
,pca_mc_A_Up05 ...
,pca_str_infix_Up05 ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_pca_D_from_mx_ver16: %0.6fs',tmp_t)); end;
AZnV_DandX_Up05_p01_pnt___(:,:,1+np_threshold) = tmp_AZnV_DandX_Up05_p01_;
%%%%%%%%;
end;%for np_threshold=0:n_p_threshold-1;
%%%%%%%%;

nf=0;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); str_home = 'data'; end;
if (strcmp(platform,'OptiPlex')); str_home = 'home'; end;
if (strcmp(platform,'eval1')); str_home = 'home'; end;
if (strcmp(platform,'rusty')); str_home = 'mnt/home'; end;
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman',str_home);
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
mr_Up05_p01_continent_ = textread(sprintf('%s/dir_Up05/mr_Up05_p01_continent_.txt',dir_trunk));

fname_fig_pre = sprintf('%s/pca_proj_D_Up05txx_DandX_p01_k22_B44_AZnV_',dir_jpg); %<-- fname_fig_pre for coarse resolution of p-value-threshold. ;
%fname_fig_pre = sprintf('%s/pca_proj_D_Up05t6x_DandX_p01_k22_B44_AZnV_',dir_jpg); %<-- fname_fig_pre for finer resolution. ;
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
%%%%%%%%;
dir_pca_tmp = sprintf('%s/dir_Up05/dir_test2mds_maf01_analyze/dir_test2mds_maf01_dex_p01_D_m2r1_g050/dir_pca/dir_pca_mda',dir_trunk);
p_threshold_ = 0.05:0.05:1.00; n_p_threshold = numel(p_threshold_); %<-- coarse resolution of p-value-threshold. ;
%p_threshold_ = 0.60:0.005:0.65; n_p_threshold = numel(p_threshold_); %<-- finer resolution. ;
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
markersize_use = 8;
p_row = 4; p_col = ceil(n_p_threshold/p_row); np=0;
for np_threshold=0:n_p_threshold-1;
p_threshold = p_threshold_(1+np_threshold);
str_p_threshold = sprintf('%.2d',min(99,floor(100*p_threshold)));
subplot(p_row,p_col,1+np);np=np+1;
tmp_fname_AnV_txx = sprintf('%s/pca_proj_D_Up05t%s_DandX_p01_k2_B44_AnV_.mda',dir_pca_tmp,str_p_threshold);
tmp_fname_ZnV_txx = sprintf('%s/pca_proj_D_Up05t%s_DandX_p01_k2_B44_ZnV_.mda',dir_pca_tmp,str_p_threshold);
if exist(tmp_fname_AnV_txx,'file') & exist(tmp_fname_AnV_txx,'file');
tmp_AnV_txx__ = mda_read_r8(tmp_fname_AnV_txx);
tmp_ZnV_txx__ = mda_read_r8(tmp_fname_ZnV_txx);
tmp_AZnV_txx__ = tmp_AnV_txx__ + tmp_ZnV_txx__;
scatter(tmp_AZnV_txx__(:,1+0),tmp_AZnV_txx__(:,1+1),markersize_use,mr_Up05_p01_continent_,'filled');
end;%if exist(tmp_fname_AnV_txx,'file') & exist(tmp_fname_AnV_txx,'file');
xlabel('pc0');ylabel('pc1');title(sprintf('p<=%+0.4f',p_threshold));
end;%for np_threshold=0:n_p_threshold-1;
%%%%%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close gcf;
%%%%%%%%;

%{
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
 %}


end;%if flag_calc;



