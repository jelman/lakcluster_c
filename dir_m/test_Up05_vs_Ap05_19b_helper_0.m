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

dir_code = sprintf('/%s/rangan/dir_bcc/dir_lakcluster_c_dev',str_home);
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jelman',str_home);
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

dir_mat_replication = sprintf('%s/dir_mat_replication',dir_trunk);
if ~exist(dir_mat_replication,'dir'); disp(sprintf(' %% mkdir %s',dir_mat_replication')); mkdir(dir_mat_replication); end;
dir_jpg_replication = sprintf('%s/dir_jpg_replication',dir_trunk);
if ~exist(dir_jpg_replication,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg_replication')); mkdir(dir_jpg_replication); end;


%%%%%%%%;
% Consider first running test_Up05_vs_Ap05_19b.m ;
%%%%%%%%;

%%%%%%%%;
% Now load one of the data-files containing the projected values and labels. ;
%%%%%%%%;
%str_datafile = sprintf('trnUp05_tst_Ap05_ncontinent0');
%%%%%%%%;
% Alternatively, we can load the data from jeremy sent on 20230114. ;
%%%%%%%%;
%str_datafile = sprintf('Up05_vs_Ap05_replication_cont1');
str_datafile = sprintf('Up05_vs_Ap05_replication_cont2');

fname_mat = sprintf('%s/%s.mat',dir_mat_replication,str_datafile);
tmp_ = load(fname_mat);

%%%%%%%%;
% Now extract and plot the training and (transformed) testing data. ;
%%%%%%%%;
tmp_index_Up05_ = efind( (tmp_.mr_dvx_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==2) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==3) );
Y_dy__ = transpose(tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_Up05_,:));
label_Y_y_ = min(2,tmp_.mr_dvx_trnUp05_tstAp05_nix_(1+tmp_index_Up05_));
tmp_index_Ap05_ = efind( (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2) );
X_dx__ = transpose(tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_Ap05_,:));
label_X_x_ = tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_(1+tmp_index_Ap05_);
parameter_apm = struct('type','parameter_apm');
parameter_apm.k_use = 32; %<-- This is the initial number of nearest neighbors to use. ;
parameter_apm.k_gamma = 0.05; %<-- robust from 0.50 to 0.05 ;
parameter_apm.flag_disp = 0;
[parameter_apm,a_est_,A_est__,tmp_X_dx__] = ...
affine_point_match_1(parameter_apm,X_dx__,Y_dy__,[],[],label_X_x_,label_Y_y_);
%%%%;
% get limits from trnUp05. ;
%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,768]);
[~,AZnV_0_lim_,AZnV_0_tick_,AZnV_1_lim_,AZnV_1_tick_] = ...
test_scatter_and_heatmap_0( ...
 [] ...
, label_Y_y_ ...
, transpose(Y_dy__) ...
);
fname_fig_pre = sprintf('%s/%s_Up05_FIGA',dir_jpg_replication,str_datafile);
sgtitle(fname_fig_pre,'Interpreter','none');
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%;
% now use those limits for tstAp05. ;
%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,768]);
test_scatter_and_heatmap_0( ...
 [] ...
, label_X_x_ ...
, transpose(tmp_X_dx__) ...
, AZnV_0_lim_ ...
, AZnV_0_tick_ ...
, AZnV_1_lim_ ...
, AZnV_1_tick_ ...
);
%%%%;
fname_fig_pre = sprintf('%s/%s_apm_Ap05_FIGA',dir_jpg_replication,str_datafile);
sgtitle(fname_fig_pre,'Interpreter','none');
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;
