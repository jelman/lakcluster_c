function ...
[ ...
 parameter ...
] = ...
xxxcluster_fromdisk_uADZSZDA_ver16( ...
 parameter ...
)
% Assumes that we never shuffle a scrambled data-set. ;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if isempty(parameter); parameter = struct('type','parameter'); end;
%%%%%%%%;
if ~isfield(parameter,'dir_trunk'); parameter.dir_trunk = pwd; end;
if ~isfield(parameter,'dir_code'); parameter.dir_code = pwd; end;
if ~isfield(parameter,'str_prefix'); parameter.str_prefix = 'test'; end;
if ~isfield(parameter,'maf_lo_threshold'); parameter.maf_lo_threshold = 0.1; end;
if ~isfield(parameter,'maf_hi_threshold'); parameter.maf_hi_threshold = 0.5; end;
if ~isfield(parameter,'str_mr_0in'); parameter.str_mr_0in = ''; end;
if ~isfield(parameter,'str_mc_0in'); parameter.str_mc_0in = ''; end;
if ~isfield(parameter,'str_lak_vs_dex'); parameter.str_lak_vs_dex = 'dex'; end;
if ~isfield(parameter,'flag_reverse'); parameter.flag_reverse = 0; end;
if ~isfield(parameter,'flag_sparse_0in'); parameter.flag_sparse_0in = []; end;
if ~isfield(parameter,'kappa_squared_0in'); parameter.kappa_squared_0in = []; end;
if ~isfield(parameter,'QR_strategy'); parameter.QR_strategy = 'YnWt condense'; end;
if ~isfield(parameter,'QC_strategy'); parameter.QC_strategy = 'YnWt store one'; end;
if ~isfield(parameter,'n_bin'); parameter.n_bin = 1; end;
if ~isfield(parameter,'n_mds'); parameter.n_mds = 2; end;
if ~isfield(parameter,'ij_mds_use_'); parameter.ij_mds_use_ = [1:2]; end;
if ~isfield(parameter,'n_mds_repl'); parameter.n_mds_repl = 1; end;
if ~isfield(parameter,'gamma'); parameter.gamma = 0.002; end;
if ~isfield(parameter,'B_MLT'); parameter.B_MLT = 32; end;
if ~isfield(parameter,'Ireq'); parameter.Ireq = 0; end;
if ~isfield(parameter,'n_scramble'); parameter.n_scramble = 0; end;
if ~isfield(parameter,'scramble_out_xdrop_'); parameter.scramble_out_xdrop_ = {}; end;
if ~isfield(parameter,'scramble_rseed_'); parameter.scramble_rseed_ = []; end;
if ~isfield(parameter,'nshuffle'); parameter.nshuffle = 0; end;
if ~isfield(parameter,'flag_rerun'); parameter.flag_rerun = 0; end;
if ~isfield(parameter,'slurm_walltime'); parameter.slurm_walltime = 0; end;
if ~isfield(parameter,'slurm_nnodes'); parameter.slurm_nnodes = 1; end;
if ~isfield(parameter,'slurm_tpn'); parameter.slurm_tpn = 15; end;
if ~isfield(parameter,'slurm_memdecl'); parameter.slurm_memdecl = 32; end;
if ~isfield(parameter,'row_factor'); parameter.row_factor = 1.0; end;
if ~isfield(parameter,'col_factor'); parameter.col_factor = 1.0; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
dir_trunk = parameter.dir_trunk;
dir_code = parameter.dir_code;
str_prefix = parameter.str_prefix;
maf_lo_threshold = parameter.maf_lo_threshold;
maf_hi_threshold = parameter.maf_hi_threshold;
str_mr_0in = parameter.str_mr_0in;
str_mc_0in = parameter.str_mc_0in;
str_lak_vs_dex = parameter.str_lak_vs_dex;
flag_reverse = parameter.flag_reverse;
flag_sparse_0in = parameter.flag_sparse_0in;
kappa_squared_0in = parameter.kappa_squared_0in;
QR_strategy = parameter.QR_strategy;
QC_strategy = parameter.QC_strategy;
n_bin = parameter.n_bin;
n_mds = parameter.n_mds;
ij_mds_use_ = parameter.ij_mds_use_;
n_mds_repl = parameter.n_mds_repl;
gamma = parameter.gamma;
B_MLT = parameter.B_MLT;
Ireq = parameter.Ireq;
n_scramble = parameter.n_scramble;
scramble_out_xdrop_ = parameter.scramble_out_xdrop_;
scramble_rseed_ = parameter.scramble_rseed_;
nshuffle = parameter.nshuffle;
flag_rerun = parameter.flag_rerun;
slurm_walltime = parameter.slurm_walltime;
slurm_nnodes = parameter.slurm_nnodes;
slurm_tpn = parameter.slurm_tpn;
slurm_memdecl = parameter.slurm_memdecl;
row_factor = parameter.row_factor;
col_factor = parameter.col_factor;
flag_verbose = parameter.flag_verbose;
%%%%%%%%;
%{
parameter.dir_trunk = dir_trunk;
parameter.dir_code = dir_code;
parameter.str_prefix = str_prefix;
parameter.maf_lo_threshold = maf_lo_threshold;
parameter.maf_hi_threshold = maf_hi_threshold;
parameter.str_mr_0in = str_mr_0in;
parameter.str_mc_0in = str_mc_0in;
parameter.str_lak_vs_dex = str_lak_vs_dex;
parameter.flag_reverse = flag_reverse;
parameter.flag_sparse_0in = flag_sparse_0in;
parameter.kappa_squared_0in = kappa_squared_0in;
parameter.QR_strategy = QR_strategy;
parameter.QC_strategy = QC_strategy;
parameter.n_bin = n_bin;
parameter.n_mds = n_mds;
parameter.ij_mds_use_ = ij_mds_use_;
parameter.n_mds_repl = n_mds_repl;
parameter.gamma = gamma;
parameter.B_MLT = B_MLT;
parameter.Ireq = Ireq;
parameter.n_scramble = n_scramble;
parameter.scramble_out_xdrop_ = scramble_out_xdrop_;
parameter.scramble_rseed_ = scramble_rseed_;
parameter.nshuffle = nshuffle;
parameter.flag_rerun = flag_rerun;
parameter.slurm_walltime = slurm_walltime;
parameter.slurm_nnodes = slurm_nnodes;
parameter.slurm_tpn = slurm_tpn;
parameter.slurm_memdecl = slurm_memdecl;
parameter.row_factor = row_factor;
parameter.col_factor = col_factor;
parameter.flag_verbose = flag_verbose;
 %}

GLOBAL_memory_gb = slurm_memdecl;

if (nshuffle>0 & n_scramble>0);
disp(sprintf(' %% Warning, n_scramble %d nshuffle %d in xxxcluster_fromdisk_uADZSZDA_ver16. Not yet implemented appropriately.',n_scramble,nshuffle));
end;%if (nshuffle>0 & n_scramble>0);

bitj = 16;
if (n_bin==1); Ireq = 0; end;

parameter_s0000 = parameter; parameter_s0000.nshuffle = 0;
str_xfix_s0000 = xxxcluster_fromdisk_uADZSZDA_xfix_gen_ver16(parameter_s0000);
str_name_s0000 = sprintf('%s_%s',str_prefix,str_xfix_s0000);
str_suffix = sprintf('%s','analyze');
disp(sprintf(' %% str_name_s0000: %s',str_name_s0000));
dir_0in = sprintf('%s/dir_%s',dir_trunk,str_prefix);
dir_tmp = sprintf('%s_%s',dir_0in,str_suffix);
if ~exist(dir_tmp,'dir'); disp(sprintf(' %% mkdir %s',dir_tmp)); mkdir(dir_tmp); end;
dir_out_s0000 = sprintf('%s_%s/dir_%s',dir_0in,str_suffix,str_name_s0000);
if ~exist(dir_out_s0000,'dir'); disp(sprintf(' %% mkdir %s',dir_out_s0000)); mkdir(dir_out_s0000); end;
dir_out_trace = sprintf('%s/dir_trace',dir_out_s0000);
if ~exist(dir_out_trace,'dir'); disp(sprintf(' %% mkdir %s',dir_out_trace)); mkdir(dir_out_trace); end;
str_timing_s0000 = sprintf('%s/timing.m',dir_out_s0000);
flag_exist_timing = 0;
if (exist(str_timing_s0000,'file')); run(str_timing_s0000); elrt_s0000 = elrt; flag_exist_timing = 1; end;

str_xfix_sxxxx = xxxcluster_fromdisk_uADZSZDA_xfix_gen_ver16(parameter);
str_name_sxxxx = sprintf('%s_%s',str_prefix,str_xfix_sxxxx);
str_suffix = sprintf('%s','analyze');
disp(sprintf(' %% str_name_sxxxx: %s',str_name_sxxxx));
dir_0in = sprintf('%s/dir_%s',dir_trunk,str_prefix);
dir_tmp = sprintf('%s_%s',dir_0in,str_suffix);
if ~exist(dir_tmp,'dir'); disp(sprintf(' %% mkdir %s',dir_tmp)); mkdir(dir_tmp); end;
dir_out = sprintf('%s_%s/dir_%s',dir_0in,str_suffix,str_name_sxxxx);
if ~exist(dir_out,'dir'); disp(sprintf(' %% mkdir %s',dir_out)); mkdir(dir_out); end;
flag_found = 0; 
str_timing = sprintf('%s/timing.m',dir_out);
str_trace = sprintf('%s/out_trace_s%0.4d.txt',dir_out_trace,nshuffle);
str_xdrop = sprintf('%s/out_xdrop_a_s%0.4d.txt',dir_out_trace,nshuffle);
flag_exist_all = exist(str_timing,'file') | (exist(str_trace,'file') & exist(str_xdrop,'file'));
if (~flag_rerun & flag_exist_all); 
if (exist(str_timing,'file')); disp(sprintf(' %% found %s, not rerunning.',str_timing)); end;
if (exist(str_trace,'file') & exist(str_xdrop,'file')) disp(sprintf(' %% found %s and\n %% found %s, not rerunning.',str_trace,str_xdrop)); end;
flag_found = 1; 
end; %if (~flag_rerun & flag_exist_all);
if (flag_rerun | ~flag_exist_all);
disp(sprintf(' %% could not find one of: \n %% %s\n %% %s\n %% %s, rerunning.',str_timing,str_trace,str_xdrop)); 
flag_found = 0; 
end;%if (flag_rerun | ~flag_exist_all);

if ~flag_found

dir_0in_plus_prefix = sprintf('%s/%s',dir_0in,str_prefix); 
dir_out_plus_prefix = sprintf('%s/%s',dir_out,str_name_sxxxx);

Y_n_cols=0; %<-- not using genetic controls. ;

% reading original row-masks for A and Z ;
mr_A_ori_ = cell(n_bin,1);
mr_Z_ori_ = cell(n_bin,1);
A_n_rij_ = cell(n_bin,1);
Z_n_rij_ = cell(n_bin,1);
for nb=0:n_bin-1;
if (numel(str_mr_0in)>0);
if (n_bin==1); str_tmp = sprintf('%s_mr_A_%s_full.b16',dir_0in_plus_prefix,str_mr_0in); else; str_tmp = sprintf('%s_mr_A_%s_%0.2d.b16',dir_0in_plus_prefix,str_mr_0in,1+nb); end;
 else;
if (n_bin==1); str_tmp = sprintf('%s_mr_A_full.b16',dir_0in_plus_prefix); else; str_tmp = sprintf('%s_mr_A_%0.2d.b16',dir_0in_plus_prefix,1+nb); end;
end;%if (numel(str_mr_0in)>0);
[bitj,nrows,ncols] = binary_getsize(str_tmp); 
mr_A_ori_{1+nb} = (binary_uncompress(str_tmp,1:nrows,1:ncols)>0);
A_n_rij_{1+nb} = find(mr_A_ori_{1+nb});
if (numel(str_mr_0in)>0);
if (n_bin==1); str_tmp = sprintf('%s_mr_Z_%s_full.b16',dir_0in_plus_prefix,str_mr_0in); else; str_tmp = sprintf('%s_mr_Z_%s_%0.2d.b16',dir_0in_plus_prefix,str_mr_0in,1+nb); end;
 else;
if (n_bin==1); str_tmp = sprintf('%s_mr_Z_full.b16',dir_0in_plus_prefix); else; str_tmp = sprintf('%s_mr_Z_%0.2d.b16',dir_0in_plus_prefix,1+nb); end;
end;%if (numel(str_mr_0in)>0);
[bitj,nrows,ncols] = binary_getsize(str_tmp); 
mr_Z_ori_{1+nb} = (binary_uncompress(str_tmp,1:nrows,1:ncols)>0);
Z_n_rij_{1+nb} = find(mr_Z_ori_{1+nb});
end;%for nb=0:n_bin-1;

% reading original col-masks for A ;
if (numel(str_mc_0in)>0);
str_tmp = sprintf('%s_mc_A_%s.b16',dir_0in_plus_prefix,str_mc_0in);[bitj,nrows,ncols] = binary_getsize(str_tmp); 
mc_A_pre = (binary_uncompress(str_tmp,1:nrows,1:ncols)>0);
else;
str_tmp = sprintf('%s_mc_A.b16',dir_0in_plus_prefix);[bitj,nrows,ncols] = binary_getsize(str_tmp); 
mc_A_pre = (binary_uncompress(str_tmp,1:nrows,1:ncols)>0);
end;%if (numel(str_mc_0in)>0);
A_n_cij = find(mc_A_pre);

M_n_rows_ = zeros(n_bin,1);
A_n_rows_ = zeros(n_bin,1);
Z_n_rows_ = zeros(n_bin,1);
for nb=0:n_bin-1;
M_n_rows_(1+nb) = numel(mr_A_ori_{1+nb});
A_n_rows_(1+nb) = sum(mr_A_ori_{1+nb});
Z_n_rows_(1+nb) = sum(mr_Z_ori_{1+nb});
end;%for nb=0:n_bin-1;
M_n_cols = numel(mc_A_pre);
disp(sprintf(' %% n_bin %d; total vs cases vs controls',n_bin));
disp(num2str([M_n_rows_ , A_n_rows_ , Z_n_rows_]));

flag_T = 0;

% compressing T_n_ in case where numel(ij_mds_use_)~=2 or n_mds_repl<1;
if (n_mds_repl<1 | numel(ij_mds_use_)~=2);
flag_T = 0;
T_n_ = cell(n_bin,1);
T_n_crop_ = cell(n_bin,1);
T_t_crop_ = cell(n_bin,1);
T_n_crop_cols = 1+numel(ij_mds_use_);
for nb=0:n_bin-1;
if (n_bin==1); str_tmp = sprintf('%s_T_full_n.b16',dir_0in_plus_prefix); else; str_tmp = sprintf('%s_T_%0.2d_n.b16',dir_0in_plus_prefix,1+nb); end;
[bitj,nrows,ncols] = binary_getsize(str_tmp); disp(sprintf(' %% reading %s = (%d,%d)',str_tmp,nrows,ncols));
T_n_{1+nb} = binary_uncompress(str_tmp,1:nrows,1:ncols); 
T_n_crop_{1+nb} = T_n_{1+nb}(:,[1,(1+(ij_mds_use_))]); 
T_t_crop_{1+nb} = transpose(T_n_crop_{1+nb});
if (n_bin==1); str_tmp = sprintf('%s_T_crop_full_n.b16',dir_out_plus_prefix); else; str_tmp = sprintf('%s_T_crop_%0.2d_n.b16',dir_out_plus_prefix,1+nb); end;
disp(sprintf(' %% writing %s = (%d,%d)',str_tmp,size(T_n_crop_{1+nb})));
binary_compress(bitj,T_n_crop_{1+nb}>0,str_tmp);
if (n_bin==1); str_tmp = sprintf('%s_T_crop_full_t.b16',dir_out_plus_prefix); else; str_tmp = sprintf('%s_T_crop_%0.2d_t.b16',dir_out_plus_prefix,1+nb); end;
disp(sprintf(' %% writing %s = (%d,%d)',str_tmp,size(T_t_crop_{1+nb})));
binary_compress(bitj,T_t_crop_{1+nb}>0,str_tmp);
end;%for nb=0:n_bin-1;
mc_T_crop = ones(T_n_crop_cols,1);
str_tmp = sprintf('%s_mc_T_crop.b16',dir_out_plus_prefix); 
disp(sprintf(' %% writing %s = (%d,%d)',str_tmp,size(mc_T_crop))); binary_compress(bitj,mc_T_crop(:)>0,str_tmp);
disp(sprintf('mc_T_crop:'));disp(num2str(transpose(mc_T_crop)));
T_n_crop_cij = 1:T_n_crop_cols;
end;%if (n_mds_repl<1 | numel(ij_mds_use_)~=2);

% copying T_n_ from T_m2rx_n in case where numel(ij_mds_use_)==2 and n_mds_repl>0;
if (n_mds_repl>=1 & numel(ij_mds_use_)==2);
flag_T = 1;
mds_str = sprintf('m%dr%d',numel(ij_mds_use_),n_mds_repl);
if (isempty(kappa_squared_0in) | kappa_squared_0in<=0);
mds_kappa_squared = textread(sprintf('%s_T_%s_kappa_squared.txt',dir_0in_plus_prefix,mds_str));
else;
mds_kappa_squared = kappa_squared_0in;
end;%if (isempty(kappa_squared_0in) | kappa_squared_0in<=0);
T_n_ = cell(n_bin,1);
T_t_ = cell(n_bin,1);
T_n_cols = 1+numel(ij_mds_use_)*n_mds_repl;
for nb=0:n_bin-1;
if (n_bin==1); str_tmp = sprintf('%s_T_%s_full_n.b16',dir_0in_plus_prefix,mds_str); else; str_tmp = sprintf('%s_T_%s_%0.2d_n.b16',dir_0in_plus_prefix,mds_str,1+nb); end;
[bitj,nrows,ncols] = binary_getsize(str_tmp); disp(sprintf(' %% reading %s = (%d,%d)',str_tmp,nrows,ncols));
T_n_{1+nb} = binary_uncompress(str_tmp,1:nrows,1:ncols); 
T_t_{1+nb} = transpose(T_n_{1+nb});
end;%for nb=0:n_bin-1;
mc_T = ones(T_n_cols,1);
str_tmp = sprintf('%s_mc_T_%s.b16',dir_out_plus_prefix,mds_str); 
disp(sprintf(' %% writing %s = (%d,%d)',str_tmp,size(mc_T))); binary_compress(bitj,mc_T(:)>0,str_tmp);
disp(sprintf('mc_T:'));disp(num2str(transpose(mc_T)));
T_n_cij = 1:T_n_cols;
end;%if (n_mds_repl>= | numel(ij_mds_use_)==2);

% checking consistency ;
if (flag_T==0); flag_error = xxxcluster_fromdisk_uADZSZDA_check_ver2(M_n_rows_,M_n_cols,A_n_rij_,A_n_cij,Z_n_rij_,T_n_crop_cols,T_n_crop_,T_n_crop_cij); end;
if (flag_T==1); flag_error = xxxcluster_fromdisk_uADZSZDA_check_ver2(M_n_rows_,M_n_cols,A_n_rij_,A_n_cij,Z_n_rij_,T_n_cols,T_n_,T_n_cij); end;
if (flag_error); disp(sprintf(' %% Warning, incorrect dimensions in xxxcluster_fromdisk_uADZSZDA_ver16')); return; end;

% performing covariate-respecting shuffle ;
if (flag_T==0); [mr_A_prm_,mr_Z_prm_] = xxxcluster_uADZSZDA_shuffle_ver2(nshuffle,M_n_rows_,M_n_cols,A_n_rij_,A_n_cij,Z_n_rij_,T_n_crop_cols,T_n_crop_,T_n_crop_cij); end;
if (flag_T==1); [mr_A_prm_,mr_Z_prm_] = xxxcluster_uADZSZDA_shuffle_ver2(nshuffle,M_n_rows_,M_n_cols,A_n_rij_,A_n_cij,Z_n_rij_,T_n_cols,T_n_,T_n_cij); end;

% writing row-masks ;
A_n_rows_used=0;
Z_n_rows_used=0;
mr_A_use_ = mr_A_ori_; mr_Z_use_ = mr_Z_ori_;
if (nshuffle>0); mr_A_use_ = mr_A_prm_; mr_Z_use_ = mr_Z_prm_; end;%if (nshuffle>0); 
for nb=0:n_bin-1;
if (row_factor<1); mr_A_use_{1+nb} = mr_A_use_{1+nb}.*(rand(size(mr_A_use_{1+nb}))<row_factor); mr_Z_use_{1+nb} = mr_Z_use_{1+nb}.*(rand(size(mr_Z_use_{1+nb}))<row_factor); end;%if (row_factor<1);
disp(sprintf('nb %.2d : mr_A_ori_ npats %.5d ncase %.4d, mr_A_use_ npats %.5d ncase %.4d, overlap %.4d',nb,numel(mr_A_ori_{1+nb}),sum(mr_A_ori_{1+nb}),numel(mr_A_use_{1+nb}),sum(mr_A_use_{1+nb}),sum(mr_A_ori_{1+nb}.*mr_A_use_{1+nb})));
disp(sprintf('nb %.2d : mr_Z_ori_ npats %.5d nctrl %.4d, mr_Z_use_ npats %.5d nctrl %.4d, overlap %.4d',nb,numel(mr_Z_ori_{1+nb}),sum(mr_Z_ori_{1+nb}),numel(mr_Z_use_{1+nb}),sum(mr_Z_use_{1+nb}),sum(mr_Z_ori_{1+nb}.*mr_Z_use_{1+nb})));
if (n_bin==1); str_tmp = sprintf('%s_mr_A_full.b16',dir_out_plus_prefix); else; str_tmp = sprintf('%s_mr_A_%0.2d.b16',dir_out_plus_prefix,1+nb); end;
binary_compress(bitj,mr_A_use_{1+nb}(:)>0,str_tmp);
if (n_bin==1); str_tmp = sprintf('%s_mr_Z_full.b16',dir_out_plus_prefix); else; str_tmp = sprintf('%s_mr_Z_%0.2d.b16',dir_out_plus_prefix,1+nb); end;
binary_compress(bitj,mr_Z_use_{1+nb}(:)>0,str_tmp);
A_n_rows_used = A_n_rows_used + sum(mr_A_use_{1+nb}(:)>0);
Z_n_rows_used = Z_n_rows_used + sum(mr_Z_use_{1+nb}(:)>0);
end;%for nb=0:n_bin-1;

% writing col-mask ;
mc_A_use = mc_A_pre;
fname_bim = sprintf('%s/%s_bim.ex2',dir_0in,str_prefix);
mc_A_bim = mc_from_bim_ext_ver5(fname_bim,maf_lo_threshold,maf_hi_threshold);
mc_A_use = mc_A_use .* mc_A_bim;
disp(sprintf(' %% maf_lo_threshold %0.2f maf_hi_threshold %0.2f, retaining %d mc values, but setting %d mc values to 0',maf_lo_threshold,maf_hi_threshold,sum(mc_A_use),sum(~mc_A_use)));
if (col_factor<1);
mc_A_use = mc_A_use.*(rand(size(mc_A_use))<col_factor); 
disp(sprintf(' %% col_factor %0.2f, retaining %d mc values, but setting %d mc values to 0',col_factor,sum(mc_A_use),sum(~mc_A_use)));
end;%if (col_factor<1);
fname_mc_A = sprintf('%s_mc_A.b16',dir_out_plus_prefix);
disp(sprintf(' %% creating %s',fname_mc_A));
str_tmp_out = fname_mc_A; binary_compress(bitj,mc_A_use(:)>0,str_tmp_out); 
disp(sprintf('mc_A nsnps %.9d/%.9d',sum(mc_A_use(:)>0),numel(mc_A_use)));
A_n_cols_used = sum(mc_A_use(:)>0);

disp(sprintf(' %% A_n_rows_used %d Z_n_rows_used %d A_n_cols_used %d',A_n_rows_used,Z_n_rows_used,A_n_cols_used));

fname_0in = sprintf('%s.in',dir_out_plus_prefix);
fp = fopen(fname_0in,'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fprintf(fp,'GLOBAL_verbose= %d;\n',flag_verbose);
fprintf(fp,'GLOBAL_thread_count= 15;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %scluster_driver;\n',str_lak_vs_dex);
if (strcmp(str_lak_vs_dex,'lak'));
if (isempty(flag_sparse_0in) | flag_sparse_0in~=1);
fprintf(fp,'GLOBAL_TEST_sparse= %d;\n',flag_sparse_0in);
end;%if (isempty(flag_sparse_0in) | flag_sparse_0in~=1);
fprintf(fp,'GLOBAL_QR_strategy= %s;\n',QR_strategy);
fprintf(fp,'GLOBAL_QC_strategy= %s;\n',QC_strategy);
end;%if (strcmp(str_lak_vs_dex,'lak'));
fprintf(fp,'GLOBAL_NBINS= %d;\n',n_bin);
if ( (B_MLT>0) & (B_MLT~=32) );
fprintf(fp,'GLOBAL_B_MLT= %d;\n',B_MLT);
end;%if ( (B_MLT>0) & (B_MLT~=32) );
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
if (flag_T==1); fprintf(fp,'GLOBAL_kappa_squared= %0.16f;\n',mds_kappa_squared); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if (nshuffle==0); A_name_string = '_A_'; end;%if (nshuffle==0);
if (nshuffle> 0); A_name_string = '_A_'; end;%if (nshuffle> 0);
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_name_= ',n_bin,dir_0in_plus_prefix,A_name_string,'_n.b16');
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_A_t_name_= ',n_bin,dir_0in_plus_prefix,A_name_string,'_t.b16');
xxxcluster_fromdisk_uADZSZDA_excerpt_1(fp,'GLOBAL_A_n_rows_= ',n_bin,M_n_rows_);
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',M_n_cols);
if (flag_reverse==1); xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',n_bin,dir_out_plus_prefix,'_mr_Z_','.b16'); end;
if (flag_reverse==0); xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',n_bin,dir_out_plus_prefix,'_mr_A_','.b16'); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A.b16;\n',dir_out_plus_prefix);
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_name_= ',n_bin,dir_0in_plus_prefix,A_name_string,'_n.b16');
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_t_name_= ',n_bin,dir_0in_plus_prefix,A_name_string,'_t.b16');
xxxcluster_fromdisk_uADZSZDA_excerpt_1(fp,'GLOBAL_Z_n_rows_= ',n_bin,M_n_rows_);
if (flag_reverse==1); xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',n_bin,dir_out_plus_prefix,'_mr_A_','.b16'); end;
if (flag_reverse==0); xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',n_bin,dir_out_plus_prefix,'_mr_Z_','.b16'); end;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',Y_n_cols);
if (flag_T==0);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',T_n_crop_cols); 
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_T_n_name_= ',n_bin,dir_out_plus_prefix,'_T_crop_','_n.b16');
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_T_t_name_= ',n_bin,dir_out_plus_prefix,'_T_crop_','_t.b16');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T_crop.b16;\n',dir_out_plus_prefix);
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_S_n_name_= ',n_bin,dir_out_plus_prefix,'_T_crop_','_n.b16');
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_S_t_name_= ',n_bin,dir_out_plus_prefix,'_T_crop_','_t.b16');
end;%if (flag_T==0);
if (flag_T==1);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',T_n_cols); 
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_T_n_name_= ',n_bin,dir_0in_plus_prefix,sprintf('_T_%s_',mds_str),'_n.b16');
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_T_t_name_= ',n_bin,dir_0in_plus_prefix,sprintf('_T_%s_',mds_str),'_t.b16');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T_%s.b16;\n',dir_out_plus_prefix,mds_str);
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_S_n_name_= ',n_bin,dir_0in_plus_prefix,sprintf('_T_%s_',mds_str),'_n.b16');
xxxcluster_fromdisk_uADZSZDA_excerpt_0(fp,'GLOBAL_S_t_name_= ',n_bin,dir_0in_plus_prefix,sprintf('_T_%s_',mds_str),'_t.b16');
end;%if (flag_T==1);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
if (n_scramble>0);
fprintf(fp,'GLOBAL_scramble_num= %d;\n',n_scramble);
fprintf(fp,'GLOBAL_scramble_out_xdrop_= ');
for nscramble=0:n_scramble-1;
fprintf(fp,'%s',scramble_out_xdrop_{1+nscramble});
if (nscramble<n_scramble-1);fprintf(fp,'%s',','); end;
if (nscramble==n_scramble-1);fprintf(fp,'%s',';'); end;
end;%for nscramble=0:n_scramble-1;
fprintf(fp,'\n');
fprintf(fp,'GLOBAL_scramble_rseed_= ');
for nscramble=0:n_scramble-1;
fprintf(fp,'%d',scramble_rseed_(1+nscramble));
if (nscramble<n_scramble-1);fprintf(fp,'%s',','); end;
if (nscramble==n_scramble-1);fprintf(fp,'%s',';'); end;
end;%for nscramble=0:n_scramble-1;
fprintf(fp,'\n');
end;%if (n_scramble>0);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by xxxcluster_fromdisk_uADZSZDA_ver16.m on %s;\n',date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fclose(fp);
type(fname_0in);

flag_call=1*(slurm_walltime<=0);%flag_call = input(' call? 1=yes (default), 0=no:'); if isempty(flag_call); flag_call=1; end;
if flag_call;
disp(sprintf('%s/lakcluster_ver18 --GLOBAL_memory_gb %d < %s',dir_code,GLOBAL_memory_gb,fname_0in));
system(sprintf('%s/lakcluster_ver18 --GLOBAL_memory_gb %d < %s',dir_code,GLOBAL_memory_gb,fname_0in));
disp(sprintf(' %% cleaning up: '));
str_command = sprintf('scp -p %s/out_trace.txt %s/out_trace_s%0.4d.txt; scp -p %s/out_xdrop_a.txt %s/out_xdrop_a_s%0.4d.txt;',dir_out,dir_out_trace,nshuffle,dir_out,dir_out_trace,nshuffle);
disp(sprintf('%s',str_command));
system(str_command);
if (nshuffle>0);
str_command = sprintf('rm -rf %s;',dir_out);
disp(sprintf('%s',str_command));
system(str_command);
end;%if (nshuffle>0);
end;%if flag_call;

flag_slurm = 1*(slurm_walltime>0);
if flag_slurm;
slurm_fname = sprintf('%s.slurm',dir_out_plus_prefix);
slurm_fp = fopen(slurm_fname,'w');
fprintf(slurm_fp,'#!/bin/sh \n');
fprintf(slurm_fp,'#\n');
fprintf(slurm_fp,'#SBATCH --verbose\n');
fprintf(slurm_fp,'#SBATCH --job-name=%s\n',fname_0in);
fprintf(slurm_fp,'#SBATCH --output=%s_output.log\n',dir_out_plus_prefix);
fprintf(slurm_fp,'#SBATCH --error=%s_error.log\n',dir_out_plus_prefix);
slurm_walltime_use = slurm_walltime;
if (flag_exist_timing); slurm_walltime_use = 1.5*elrt_s0000/3600; disp(sprintf(' %% slurm_walltime changed from %0.2f to %0.2f',slurm_walltime,slurm_walltime_use)); end;
slurm_walltime_h = floor(slurm_walltime_use); slurm_walltime_m = min(59,ceil(60*(slurm_walltime_use - slurm_walltime_h)));
sprintf(' %% slurm_walltime=%d:%.2d:59',slurm_walltime_h,slurm_walltime_m);
fprintf(slurm_fp,'#SBATCH --time=%d:%.2d:59\n',slurm_walltime_h,slurm_walltime_m);
fprintf(slurm_fp,'#SBATCH --nodes=%d --ntasks-per-node=%d\n',slurm_nnodes,slurm_tpn);
fprintf(slurm_fp,'#SBATCH --mem=%dGB\n',slurm_memdecl);
fprintf(slurm_fp,'\n');
fprintf(slurm_fp,'/bin/hostname\n');
fprintf(slurm_fp,'/bin/pwd\n');
fprintf(slurm_fp,'module load matlab/2017b\n');
fprintf(slurm_fp,'%s/lakcluster_ver18 --GLOBAL_memory_gb %d < %s\n',dir_code,GLOBAL_memory_gb,fname_0in);
%%%%%%%%%%%%%%%%;
str_command = sprintf('scp -p %s/out_trace.txt %s/out_trace_s%0.4d.txt; scp -p %s/out_xdrop_a.txt %s/out_xdrop_a_s%0.4d.txt;',dir_out,dir_out_trace,nshuffle,dir_out,dir_out_trace,nshuffle);
fprintf(slurm_fp,'%s\n',str_command);
if (nshuffle>0);
str_command = sprintf('rm -rf %s;',dir_out);
fprintf(slurm_fp,'%s\n',str_command);
end;%if (nshuffle>0);
%%%%%%%%%%%%%%%%;
fclose(slurm_fp);
type(slurm_fname);
str_command = sprintf('sbatch %s;',slurm_fname);
fp = fopen(sprintf('%s/log.txt',dir_trunk),'a'); fprintf(fp,'%s # %s\n',str_command,datestr(now)); fclose(fp);
disp(sprintf('%s',str_command));
%system(sprintf('%s\n',str_command));
end;%if flag_slurm;

end;%if ~flag_found;

