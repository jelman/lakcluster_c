function ...
[ ...
 parameter ...
,AZnV_ ...
,AnV_ ...
,ZnV_ ...
,V_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_proj_from_mx_ver16( ...
 parameter ...
,pca_rank ...
,pca_str_V ...
,pca_mr_A_ ...
,pca_mr_Z_ ...
,pca_mc_A ...
,pca_str_infix ...
,mx__ ...
);

str_thisfunction = 'xxxcluster_fromdisk_uADZSZDA_pca_proj_from_mx_ver16';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); pca_rank=[]; end; na=na+1;
if (nargin<1+na); pca_str_V=[]; end; na=na+1;
if (nargin<1+na); pca_mr_A_=[]; end; na=na+1;
if (nargin<1+na); pca_mr_Z_=[]; end; na=na+1;
if (nargin<1+na); pca_mc_A=[]; end; na=na+1;
if (nargin<1+na); pca_str_infix=[]; end; na=na+1;
if (nargin<1+na); mx__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'flag_reverse'); parameter.flag_reverse = 0; end;
if ~isfield(parameter,'str_A_p'); parameter.str_A_p = ''; end;
flag_verbose = parameter.flag_verbose;
flag_reverse = parameter.flag_reverse;
str_A_p = parameter.str_A_p;

if isempty(mx__); mx__ = load_mx__from_parameter_ver0(parameter); end;
if isempty(pca_str_V); pca_str_V = zeros(size(mx__.mc_A_,1),pca_rank); end;
%%%%;
if ~isfield(parameter,'dir_0in'); disp(sprintf(' %% Warning, parameter.dir_0in undefined in %s',str_thisfunction)); end;
dir_0in = parameter.dir_0in;
str_suffix = sprintf('%s','analyze');
dir_tmp = sprintf('%s_%s',dir_0in,str_suffix);
if ~exist(dir_tmp,'dir'); disp(sprintf(' %% mkdir %s',dir_tmp)); mkdir(dir_tmp); end;
%%%%;
if ~isfield(parameter,'str_name_s0000');
parameter.str_name_s0000 = 'default';
end;%if ~isfield(parameter,'str_name_s0000');
str_name_s0000 = parameter.str_name_s0000;
%%%%;
if ~isfield(parameter,'dir_out_s0000');
str_suffix = sprintf('%s','analyze');
dir_out_s0000 = sprintf('%s_%s/dir_%s',dir_0in,str_suffix,str_name_s0000);
if ~exist(dir_out_s0000,'dir'); disp(sprintf(' %% mkdir %s',dir_out_s0000)); mkdir(dir_out_s0000); end;
parameter.dir_out_s0000 = dir_out_s0000;
end;%if ~isfield(parameter,'dir_out_s0000');
%%%%;
str_name_s0000 = parameter.str_name_s0000;
dir_out_s0000 = parameter.dir_out_s0000;
if (flag_verbose); disp(sprintf(' %% str_name_s0000: %s',str_name_s0000)); end;
if (flag_verbose); disp(sprintf(' %% dir_out_s0000: %s',dir_out_s0000)); end;

%%%%;
if  isfield(parameter,'n_study'); n_study = parameter.n_study; end;
if ~isfield(parameter,'n_study'); n_study = mx__.n_study; end;

flag_force_create=0;
if isempty(pca_str_infix);
flag_force_create=1;
str_tmp = sprintf('%.16d',floor(1e16*rem(now,1)));
pca_str_infix = sprintf('deleteme_%s',str_tmp);
end;%if isempty(pca_str_infix);
if (flag_verbose); disp(sprintf(' %% pca_str_infix: %s',pca_str_infix)); end;

flag_reverse = parameter.flag_reverse;

%%%%%%%%;
% write and/or read V_. ;
%%%%%%%%;
if ~isnumeric(pca_str_V) & ~isstring(pca_str_V); 
disp(sprintf(' %% Warning, pca_str_V neither numeric nor string in %s',str_thisfunction));
disp('returning'); return;
end;%if ~isnumeric(pca_str_V) & ~isstring(pca_str_V); 
%%%%;
if ~isnumeric(pca_str_V) &  isstring(pca_str_V); 
if (flag_verbose); disp(sprintf(' %% pca_str_V a string in %s',str_thisfunction)); end;
end;%if ~isnumeric(pca_str_V) &  isstring(pca_str_V); 
%%%%;
if  isnumeric(pca_str_V) & ~isstring(pca_str_V); 
if (flag_verbose); disp(sprintf(' %% pca_str_V numeric in %s',str_thisfunction)); end;
V_ = pca_str_V;
if size(V_,1)~=size(mx__.mc_A_,1);
disp(sprintf(' %% Warning, V_ incompatible with mc_A_ in %s',str_thisfunction));
disp('returning'); return;
end;%if size(V_,1)~=size(mx__.mc_A_,1);
if  isempty(pca_rank); pca_rank = size(V_,2); end;
pca_rank = min(pca_rank,size(V_,2));
pca_str_V = sprintf('%s/%s_k%d_B44_V.mda',dir_out_s0000,str_name_s0000,pca_rank);
mda_write_d3_r8(V_(:,1:pca_rank),pca_str_V);
end;%if  isnumeric(pca_str_V) & ~isstring(pca_str_V); 
%%%%;
parameter.pca_str_V = pca_str_V;
V_ = mda_read_r8(pca_str_V);
%%%%%%%%;

%%%%%%%%;
% projection. ;
%%%%%%%%;
parameter_pca_proj = struct('type','parameter_pca_proj');
parameter_pca_proj.str_driver = 'pca_proj_driver';
parameter_pca_proj.str_V = pca_str_V;
if flag_force_create; parameter_pca_proj.flag_force_create = 1; end;
parameter_pca_proj.pca_mc_T = zeros(size(mx__.mc_T_)); parameter_pca_proj.pca_mc_T(1)=1;
parameter_pca_proj.pca_mc_A = mx__.mc_A_;
if ~isempty(pca_mc_A); parameter_pca_proj.pca_mc_A = pca_mc_A; end;
parameter_pca_proj.pca_mr_A_ = cell(n_study,1);
parameter_pca_proj.pca_mr_Z_ = cell(n_study,1);
for nstudy=0:n_study-1;
parameter_pca_proj.pca_mr_A_{1+nstudy} = mx__.mr_A__{1+nstudy};
parameter_pca_proj.pca_mr_Z_{1+nstudy} = mx__.mr_Z__{1+nstudy};
end;%for nstudy=0:n_study-1;
parameter_pca_proj.str_infix = ''; if ~isempty(pca_str_infix); parameter_pca_proj.str_infix = sprintf('pca_proj_D_%s',pca_str_infix); end;
if ~isempty(str_A_p); parameter_pca_proj.str_A_p = str_A_p; end;
if ~isempty(pca_rank); parameter_pca_proj.rank = pca_rank; end;
[ ...
 parameter ...
 parameter_pca_proj ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_ver16( ...
 parameter ...
,parameter_pca_proj ...
);
%%%%%%%%;
AnV_ = mda_read_r8(parameter_pca_proj.str_AnV);
ZnV_ = mda_read_r8(parameter_pca_proj.str_ZnV);
AZnV_ = AnV_ + ZnV_;
