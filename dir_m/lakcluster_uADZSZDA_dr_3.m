function lakcluster_uADZSZDA_dr_3(dir_code,dir_trunk,prefix,M_n_,rev_flag,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind,GLOBAL_TEST_sparse,gamma,B_MLT,Ireq,shuffle_num,verbose_flag,force_create_flag,slurm_walltime,slurm_nnodes,slurm_tpn,slurm_memdecl);
% using lakcluster_ver18 and slurm ;
% removing references to Z and S if empty. ;
% test with: ;
%{
  lakcluster_uADZSZDA_dr_3();
  %}

if (nargin<1);
disp('returning'); return;
end;%if (nargin<1);

na=1;
if (nargin<na); dir_code = pwd; end; na=na+1;
if (nargin<na); dir_trunk = pwd; end; na=na+1;
if (nargin<na); prefix = 'test'; end; na=na+1;
if (nargin<na); M_n_{1} = randn(1024); end; na=na+1;
if (nargin<na); rev_flag = 0; end; na=na+1;
if (nargin<na); A_n_rind_ = {1:512}; end; na=na+1;
if (nargin<na); A_n_cind = 1:512; end; na=na+1;
if (nargin<na); Z_n_rind_ = {1:512}; end; na=na+1;
if (nargin<na); T_n_{1} = ones(1024,1); end; na=na+1;
if (nargin<na); T_n_cind = 1; end; na=na+1;
if (nargin<na); GLOBAL_TEST_sparse = 1; end; na=na+1;
if (nargin<na); gamma = 0; end; na=na+1;
if (nargin<na); B_MLT = 32; end; na=na+1;
if (nargin<na); Ireq = 1; end; na=na+1;
if (nargin<na); shuffle_num = 0; end; na=na+1;
if (nargin<na); verbose_flag = 0; end; na=na+1;
if (nargin<na); force_create_flag = 1; end; na=na+1;
if (nargin<na); slurm_walltime = 0; end; na=na+1;
if (nargin<na); slurm_nnodes = 1; end; na=na+1;
if (nargin<na); slurm_tpn = 1; end; na=na+1;
if (nargin<na); slurm_memdecl = 1; end; na=na+1;

nbins=length(M_n_);
[A_n_cols,Y_n_cols,T_n_cols,~,~] = lakcluster_uADZSZDA_check_0(shuffle_num,M_n_,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind);

QR_strategy = 'YnWt';
QC_strategy = 'YnWt store one';
%QC_strategy = 'ZtSWn';
bitj = 16;

Z_bother = 0;
for nb1=0:nbins-1;
if (length(Z_n_rind_{1+nb1})>0); Z_bother = 1; end;
end;%for nb1=0:nbins-1;
if ~Z_bother; rev_flag = 0; end;

test_string = sprintf('%s_%s',prefix,lakcluster_uADZSZDA_xfix_gen_ver1(rev_flag,A_n_rind_,Z_n_rind_,T_n_cind,GLOBAL_TEST_sparse,gamma,B_MLT,Ireq,shuffle_num));
dir__in = sprintf('%s/dir_%s',dir_trunk,prefix);
dir_out = sprintf('%s/dir_%s',dir__in,test_string); 
disp(sprintf(' test_string: %s',test_string));
disp(sprintf(' dir__in: %s',dir__in));
disp(sprintf(' dir_out: %s',dir_out));
if ~exist(dir_out,'dir'); disp(sprintf(' %% creating %s',dir_out)); mkdir(dir_out); 
 else disp(sprintf(' %% directory %s already exists, not creating.',dir_out)); end;
found_trace_flag = 0; 
tmpchar_trace = sprintf('%s/out_trace.txt',dir_out);
if exist(tmpchar_trace,'file');
tmp_trace = textread(tmpchar_trace);
if length(tmp_trace)> 6; disp(sprintf(' %% found %s of length %d, not rerunning.',tmpchar_trace,length(tmp_trace))); found_trace_flag = 1; end;
if (force_create_flag & length(tmp_trace)> 6); disp(sprintf(' %% found %s of length %d, actually, rerunning anyway.',tmpchar_trace,length(tmp_trace))); found_trace_flag = 0; end;
if length(tmp_trace)<=6; disp(sprintf(' %% found %s of length %d,     rerunning.',tmpchar_trace,length(tmp_trace))); found_trace_flag = 0; end;
end;%if exist(tmpchar_trace,'file');

if ~found_trace_flag;

d_inpre = sprintf('%s/%s',dir__in,prefix); 
d_oupre = sprintf('%s/%s',dir_out,test_string);

for nb1=0:nbins-1;
mr_M = zeros(size(M_n_{1+nb1},1),1);
mr_A_ori_{1+nb1} = mr_M; mr_A_ori_{1+nb1}(A_n_rind_{1+nb1})=1;
tmpchar = sprintf('%s_mr_A_%d.b16',d_oupre,0+nb1); 
if (force_create_flag | ~exist(tmpchar,'file'));
disp(sprintf(' %% creating %s of size %d-x-%d',tmpchar,length(mr_A_ori_{1+nb1}),1)); 
tutorial_binary_compress(bitj,mr_A_ori_{1+nb1}(:)>0,tmpchar);
end;% if exist;
if Z_bother;
mr_Z_ori_{1+nb1} = mr_M; mr_Z_ori_{1+nb1}(Z_n_rind_{1+nb1})=1;
tmpchar = sprintf('%s_mr_Z_%d.b16',d_oupre,0+nb1); 
if (force_create_flag | ~exist(tmpchar,'file'));
disp(sprintf(' %% creating %s of size %d-x-%d',tmpchar,length(mr_Z_ori_{1+nb1}),1)); 
tutorial_binary_compress(bitj,mr_Z_ori_{1+nb1}(:)>0,tmpchar);
end;% if exist;
end;%if Z_bother;
end;%for nb1=0:nbins-1;

if (shuffle_num>0); % performing covariate-respecting shuffle ;
if Z_bother;
[mr_A_prm_,mr_Z_prm_] = lakcluster_uADZSZDA_shuffle_1(shuffle_num,M_n_,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind);
for nb1=0:nbins-1;
tmpchar = sprintf('%s_mr_A_%d.b16',d_oupre,0+nb1); tutorial_binary_compress(bitj,mr_A_prm_{1+nb1}(:)>0,tmpchar);
if Z_bother;
tmpchar = sprintf('%s_mr_Z_%d.b16',d_oupre,0+nb1); tutorial_binary_compress(bitj,mr_Z_prm_{1+nb1}(:)>0,tmpchar);
end;%if Z_bother;
end;%for nb1=0:nbins-1;
end;%if Z_bother;
if ~Z_bother;
rng(shuffle_num);
for nb1=0:nbins-1;
[tmp_Q_,~] = qr(randn(size(M_n_{1+nb1},1)));
L_n_{1+nb1} = tmp_Q_*M_n_{1+nb1};
end;%for nb1=0:nbins-1;
end;%if ~Z_bother;
end;% if (shuffle_num>0); % performing covariate-respecting shuffle ;

% write M_n_ ;
for nb1=0:nbins-1;
tmpchar = sprintf('%s_M_%d_n.b16',d_inpre,0+nb1); 
if (force_create_flag | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of length %d-x-%d',tmpchar,size(M_n_{1+nb1}))); tutorial_binary_compress(bitj,M_n_{1+nb1}>0,tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
if ~Z_bother & shuffle_num>0;
tmpchar = sprintf('%s_L_%d_n.b16',d_oupre,0+nb1); 
if (force_create_flag | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of length %d-x-%d',tmpchar,size(L_n_{1+nb1}))); tutorial_binary_compress(bitj,L_n_{1+nb1}>0,tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
end;%if ~Z_bother & shuffle_num>0;
tmpchar = sprintf('%s_M_%d_t.b16',d_inpre,0+nb1); 
if (force_create_flag | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of size %d-x-%d',tmpchar,size(transpose(M_n_{1+nb1})))); tutorial_binary_compress(bitj,transpose(M_n_{1+nb1}>0),tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
if ~Z_bother & shuffle_num>0;
tmpchar = sprintf('%s_L_%d_t.b16',d_oupre,0+nb1); 
if (force_create_flag | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of length %d-x-%d',tmpchar,size(transpose(L_n_{1+nb1})))); tutorial_binary_compress(bitj,transpose(L_n_{1+nb1}>0),tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
end;%if ~Z_bother & shuffle_num>0;
% write T_n_ ;
tmpchar = sprintf('%s_T_%d_n.b16',d_inpre,0+nb1); 
if (force_create_flag | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of size %d-x-%d',tmpchar,size(T_n_{1+nb1}))); tutorial_binary_compress(bitj,T_n_{1+nb1}>0,tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
tmpchar = sprintf('%s_T_%d_t.b16',d_inpre,0+nb1); 
if (force_create_flag | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of size %d-x-%d',tmpchar,size(transpose(T_n_{1+nb1})))); tutorial_binary_compress(bitj,transpose(T_n_{1+nb1}>0),tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
end;%for nb1=0:nbins-1;

% write A_n_cind and T_n_cind ;
mc_A = zeros(A_n_cols,1); mc_A(A_n_cind)=1;
tmpchar = sprintf('%s_mc_A.b16',d_oupre); tutorial_binary_compress(bitj,mc_A(:)>0,tmpchar);
mc_T = zeros(T_n_cols,1); mc_T(T_n_cind)=1;
tmpchar = sprintf('%s_mc_T.b16',d_oupre); tutorial_binary_compress(bitj,mc_T(:)>0,tmpchar);

fname__in = sprintf('%s.in',d_oupre);
fp = fopen(fname__in,'w');
fprintf(fp,'GLOBAL_verbose= %d;\n',verbose_flag);
fprintf(fp,'GLOBAL_thread_count= 8;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','lakcluster_driver');
fprintf(fp,'GLOBAL_QR_strategy= %s;\n',QR_strategy);
fprintf(fp,'GLOBAL_QC_strategy= %s;\n',QC_strategy);
fprintf(fp,'GLOBAL_NBINS= %d;\n',nbins);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',B_MLT);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_TEST_sparse= %d;\n',GLOBAL_TEST_sparse);
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
if Z_bother | shuffle_num==0;
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
end;%if Z_bother | shuffle_num==0;
if ~Z_bother & shuffle_num>0;
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_L_%d_n.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_L_%d_t.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
end;%if ~Z_bother & shuffle_num>0;
fprintf(fp,'GLOBAL_A_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',size(M_n_{1+nb1},1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',A_n_cols);
if (rev_flag==1); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_Z_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_A_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A.b16;\n',d_oupre);
if Z_bother;
fprintf(fp,'GLOBAL_Z_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',size(M_n_{1+nb1},1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
if (rev_flag==1); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_A_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_Z_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
end;%if Z_bother;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',Y_n_cols);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',T_n_cols); 
fprintf(fp,'GLOBAL_T_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T.b16;\n',d_oupre);
if Z_bother;
fprintf(fp,'GLOBAL_S_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_S_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
end;%if Z_bother;
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lakcluster_uADZSZDA_dr_3.m on %s;\n',date);
fclose(fp);

call_flag=1*(slurm_walltime<=0);
if call_flag;
disp(sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in));
system(sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in));
end;%if call_flag;

slurm_flag = 1*(slurm_walltime>0);
if slurm_flag;
slurm_fname = sprintf('%s.slurm',d_oupre);
slurm_fp = fopen(slurm_fname,'w');
fprintf(slurm_fp,'#!/bin/sh \n');
fprintf(slurm_fp,'#\n');
fprintf(slurm_fp,'#SBATCH --verbose\n');
fprintf(slurm_fp,'#SBATCH --job-name=%s\n',fname__in);
fprintf(slurm_fp,'#SBATCH --output=%s_output.log\n',d_oupre);
fprintf(slurm_fp,'#SBATCH --error=%s_error.log\n',d_oupre);
slurm_walltime_h = floor(slurm_walltime);
slurm_walltime_m = min(59,ceil(60*(slurm_walltime - slurm_walltime_h)));
sprintf(' %% slurm_walltime=%d:%.2d:59',slurm_walltime_h,slurm_walltime_m);
fprintf(slurm_fp,'#SBATCH --time=%d:%.2d:59\n',slurm_walltime_h,slurm_walltime_m);
fprintf(slurm_fp,'#SBATCH --nodes=%d --ntasks-per-node=%d\n',slurm_nnodes,slurm_tpn);
fprintf(slurm_fp,'#SBATCH --mem=%dGB\n',slurm_memdecl);
fprintf(slurm_fp,'\n');
fprintf(slurm_fp,'/bin/hostname\n');
fprintf(slurm_fp,'/bin/pwd\n');
fprintf(slurm_fp,'module load matlab/2016b\n');
fprintf(slurm_fp,'%s/lakcluster_ver18 < %s\n',dir_code,fname__in);
fclose(slurm_fp);
type(slurm_fname);
system(sprintf('sbatch %s;\n',slurm_fname));
end;%if slurm_flag;

end;%if ~found_trace_flag;
