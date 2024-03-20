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

if (flag_verbose); disp(sprintf(' %% Comparing Up05 with Ap05 data. ;')); end;

dir_code = '/home/jelman/Github/lakcluster_c';
dir_trunk = '/home/jelman/Projects/AD_Biclustering/data/ADNI/ADNI_vs_UKB';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

dir_mat_replication = sprintf('%s/dir_mat_replication',dir_trunk);
if ~exist(dir_mat_replication,'dir'); disp(sprintf(' %% mkdir %s',dir_mat_replication')); mkdir(dir_mat_replication); end;
dir_jpg_replication = sprintf('%s/dir_jpg_replication',dir_trunk);
if ~exist(dir_jpg_replication,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg_replication')); mkdir(dir_jpg_replication); end;

%%%%%%%%;
% Consider first running test_Up05_vs_Ap05_18c.m ;
%%%%%%%%;

%%%%%%%%;
% Now load one of the data-files containing the projected values and labels. ;
%%%%%%%%;
ncontinent = 2; 
str_datafile = sprintf('trnUp05_tst_Ap05_ncontinent%d',ncontinent); 
k_gamma_use = 0.15; %<-- Continent1 converges at .50, continent2 must be lowered to 0.15 ;
%%%%%%%%;
% Alternatively, we can load the data from jeremey sent on 20230114. ;
%%%%%%%%;
%str_datafile = sprintf('Up05_vs_Ap05_replication_cont1'); k_gamma_use = 0.50; %<-- this is quite robust for this dataset, with similar results from [0.05,0.50]. ;
%str_datafile = sprintf('Up05_vs_Ap05_replication_cont2'); k_gamma_use = 0.05; %<-- this is less robust for this dataset, with similar results from [0.05,0.25]. ;

fname_mat = sprintf('%s/%s.mat',dir_mat_replication,str_datafile);
tmp_ = load(fname_mat);

parameter_apm = struct('type','parameter_apm');
parameter_apm.k_use = 32; %<-- This is the initial number of nearest neighbors to use. ;
parameter_apm.k_gamma = k_gamma_use ;
parameter_apm.flag_disp = 0;
parameter_apm.n_shuffle = 1; %<-- n_shuffle_interior. ;
parameter_apm.n_shuffle_exterior = 500;
n_shuffle_exterior = parameter_apm.n_shuffle_exterior;

tmp_dir = sprintf('%s/dir_%s',dir_mat_replication,str_datafile);
if ~exist(tmp_dir,'dir'); disp(sprintf(' %% mkdir %s',tmp_dir)); mkdir(tmp_dir); end;

%%%%%%%%;
% Now we try and find a mapping from loaded Ap05 to loaded Up05. ;
% For Ap05 we will use the projected values: ;
% tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_. ;
% For Up05 we will use the projected values: ;
% tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_. ;
% Note that affine_point_match_1 is designed to account for the labels in the data. ;
% For this point-cloud matching we will use only case vs ctrl labels (i.e., no bicluster labels). ;
%%%%%%%%;
nshuffle_exterior = 0;
tmp_fname_pre = sprintf('%s/%s_s%.4d',tmp_dir,str_datafile,nshuffle_exterior);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre);
if ~tmp_flag_skip;
if (flag_verbose); disp(sprintf(' %% nshuffle_exterior %d/%d',nshuffle_exterior,n_shuffle_exterior)); end;
tmp_index_Up05_ = efind( (tmp_.mr_dvx_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==2) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==3) );
tmp_prm_Up05_ = transpose(1:numel(tmp_index_Up05_));
tmp_index_Ap05_ = efind( (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2) );
tmp_prm_Ap05_ = transpose(1:numel(tmp_index_Ap05_));
if (nshuffle_exterior> 0);
rng(nshuffle_exterior); tmp_prm_Up05_ = randperm(numel(tmp_prm_Up05_)); tmp_prm_Ap05_ = randperm(numel(tmp_prm_Ap05_));
end;%if (nshuffle_exterior> 0);
tmp_mr_dvx_trnUp05_tstAp05_nix_ = tmp_.mr_dvx_trnUp05_tstAp05_nix_ ;
tmp_mr_dvx_trnUp05_tstAp05_nix_(1+tmp_index_Up05_) = tmp_mr_dvx_trnUp05_tstAp05_nix_(1+tmp_index_Up05_(tmp_prm_Up05_)) ;
tmp_mr_dvx_Ap05_from_trnUp05_tstAp05_nix_ = tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_;
tmp_mr_dvx_Ap05_from_trnUp05_tstAp05_nix_(1+tmp_index_Ap05_) = tmp_mr_dvx_Ap05_from_trnUp05_tstAp05_nix_(1+tmp_index_Ap05_(tmp_prm_Ap05_)) ;
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
,tmp_mr_dvx_trnUp05_tstAp05_nix_ ...
,tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_ ...
,tmp_mr_dvx_Ap05_from_trnUp05_tstAp05_nix_ ...
,tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_ ...
);
save(tmp_fname_mat ...
     ,'nshuffle_exterior' ...
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
end;%if ~tmp_flag_skip;
%%%%%%%%;
tmp_fname_mat = sprintf('%s.mat',tmp_fname_pre);
if  exist(tmp_fname_mat,'file');
load(tmp_fname_mat,'tmp_a_est_','tmp_A_est__');
a_est_ = tmp_a_est_;
A_est__ = tmp_A_est__;
end;%if  exist(tmp_fname_mat,'file');

%%%%%%%%;
% First we just plot the loaded data: ;
%%%%%%%%;
fname_fig_pre = sprintf('%s/%s_FIGA',dir_jpg_replication,str_datafile);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
figure(1+nf);nf=nf+1;clf;figbig;
markersize_big = 6;
markersize_sml = 4;
fontsize_use = 12;
subplot(1,2,1);
hold on;
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==1);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[0.0,1.0,1.0]);
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==2);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[0.5,0.0,0.5]);
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==3);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[1.0,0.0,1.0]);
legend({'Up05 ctrl','Up05 case','Up05 bicl'},'Location','NorthWest');
grid on;
set(gca,'Fontsize',fontsize_use);
hold off;
subplot(1,2,2);
hold on;
tmp_index_ = efind(tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1);
X_xd__ = tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
plot(X_xd__(:,1+0),X_xd__(:,1+1),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.0,0.8,0.8]);
tmp_index_ = efind(tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2);
X_xd__ = tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
plot(X_xd__(:,1+0),X_xd__(:,1+1),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.8,0.0,0.8]);
legend({'Ap05 ctrl','Ap05 case'},'Location','NorthWest');
grid on;
set(gca,'Fontsize',fontsize_use);
hold off;
sgtitle(sprintf('%s',str_datafile),'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');

%%%%%%%%;
% Now we can plot the results. ;
%%%%%%%%;
fname_fig_pre = sprintf('%s/%s_FIGB',dir_jpg_replication,str_datafile);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
figure(1+nf);nf=nf+1;clf;figbig;
markersize_big = 12;
markersize_sml = 8;
fontsize_use = 12;
subplot(1,1,1);
hold on;
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==1);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[0.0,1.0,1.0]);
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==2);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[0.5,0.0,0.5]);
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==3);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[1.0,0.0,1.0]);
tmp_index_ = efind(tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1);
X_xd__ = tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
plot(tmp_X_xd__(:,1+0),tmp_X_xd__(:,1+1),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.0,0.8,0.8]);
tmp_index_ = efind(tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2);
X_xd__ = tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
plot(tmp_X_xd__(:,1+0),tmp_X_xd__(:,1+1),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.8,0.0,0.8]);
legend({'Up05 ctrl','Up05 case','Up05 bicl','Ap05 ctrl','Ap05 case'},'Location','NorthWest');
set(gca,'Fontsize',fontsize_use);
hold off;
sgtitle(sprintf('%s',str_datafile),'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');

%%%%%%%%;
% Now estimate p-value. ;
%%%%%%%%;
tmp_dir = sprintf('%s/dir_%s',dir_mat_replication,str_datafile);
if ~exist(tmp_dir,'dir'); disp(sprintf(' %% mkdir %s',tmp_dir)); mkdir(tmp_dir); end;
n_shuffle_exterior = parameter_apm.n_shuffle_exterior;
for nshuffle_exterior=0:n_shuffle_exterior;
tmp_fname_pre = sprintf('%s/%s_s%.4d',tmp_dir,str_datafile,nshuffle_exterior);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre);
if ~tmp_flag_skip;
if (flag_verbose); disp(sprintf(' %% nshuffle_exterior %d/%d',nshuffle_exterior,n_shuffle_exterior)); end;
tmp_index_Up05_ = efind( (tmp_.mr_dvx_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==2) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==3) );
tmp_prm_Up05_ = transpose(1:numel(tmp_index_Up05_));
tmp_index_Ap05_ = efind( (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2) );
tmp_prm_Ap05_ = transpose(1:numel(tmp_index_Ap05_));
if (nshuffle_exterior> 0);
rng(nshuffle_exterior); tmp_prm_Up05_ = randperm(numel(tmp_prm_Up05_)); tmp_prm_Ap05_ = randperm(numel(tmp_prm_Ap05_));
end;%if (nshuffle_exterior> 0);
tmp_mr_dvx_trnUp05_tstAp05_nix_ = tmp_.mr_dvx_trnUp05_tstAp05_nix_ ;
tmp_mr_dvx_trnUp05_tstAp05_nix_(1+tmp_index_Up05_) = tmp_mr_dvx_trnUp05_tstAp05_nix_(1+tmp_index_Up05_(tmp_prm_Up05_)) ;
tmp_mr_dvx_Ap05_from_trnUp05_tstAp05_nix_ = tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_;
tmp_mr_dvx_Ap05_from_trnUp05_tstAp05_nix_(1+tmp_index_Ap05_) = tmp_mr_dvx_Ap05_from_trnUp05_tstAp05_nix_(1+tmp_index_Ap05_(tmp_prm_Ap05_)) ;
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
,tmp_mr_dvx_trnUp05_tstAp05_nix_ ...
,tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_ ...
,tmp_mr_dvx_Ap05_from_trnUp05_tstAp05_nix_ ...
,tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_ ...
);
save(tmp_fname_mat ...
     ,'nshuffle_exterior' ...
     ,'parameter_apm' ...
     ,'tmp_f_z_' ...
     ,'tmp_f_x_' ...
     ,'tmp_f_y_' ...
     ,'tmp_a_est_' ...
     ,'tmp_A_est__' ...
     );
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
end;%for nshuffle_exterior=0:n_shuffle_exterior;

%%%%%%%%;
% Now collate the results. ;
%%%%%%%%;
n_y = numel( efind( (tmp_.mr_dvx_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==2) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==3) ) );
n_x = numel( efind( (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2) ) );
n_z = max(n_x,n_y);
n_shuffle_exterior = parameter_apm.n_shuffle_exterior;
n_shuffle_interior = parameter_apm.n_shuffle;
f_ze__ = zeros(n_z,1+n_shuffle_exterior);
f_xe__ = zeros(n_x,1+n_shuffle_exterior);
f_ye__ = zeros(n_y,1+n_shuffle_exterior);
for nshuffle_exterior=0:n_shuffle_exterior;
tmp_fname_pre = sprintf('%s/%s_s%.4d',tmp_dir,str_datafile,nshuffle_exterior);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_fig_pre = sprintf('%s/%s_FIGC',dir_jpg_replication,str_datafile);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
tmp_x_ = linspace(0,1,n_x); tmp_y_ = linspace(0,1,n_y); tmp_z_ = linspace(0,1,n_z);
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig; p_row = 3; p_col = 3;
linewidth_sml = 0.5; linewidth_big = 2.0; ylim_ = 10*[-1,+1];
%%%%;
for ntab=0:p_col-1;
if ntab==0; tmp_w_ = tmp_z_; f_we__ = f_ze__; z_we__ = z_ze__; prctile_we__ = prctile_ze__; str_f_w = 'f_z_'; str_z_w = 'z_z_'; str_p_w = 'p_z_'; end;
if ntab==1; tmp_w_ = tmp_x_; f_we__ = f_xe__; z_we__ = z_xe__; prctile_we__ = prctile_xe__; str_f_w = 'f_x_'; str_z_w = 'z_x_'; str_p_w = 'p_x_'; end;
if ntab==2; tmp_w_ = tmp_y_; f_we__ = f_ye__; z_we__ = z_ye__; prctile_we__ = prctile_ye__; str_f_w = 'f_y_'; str_z_w = 'z_y_'; str_p_w = 'p_y_'; end;
%%%%;
subplot(p_row,p_col,1+ntab+0*p_col);
hold on;
plot(tmp_w_,f_we__(:,2:end),'k','LineWidth',linewidth_sml);
plot(tmp_w_,f_we__(:,1+0)  ,'r','LineWidth',linewidth_big);
hold off;
title(str_f_w,'Interpreter','none');
xlim([0,1]); xlabel('fraction nearest');
ylabel('f fraction');
set(gca,'XTick',0:0.05:1.0); grid on; xtickangle(90);
%%%%;
subplot(p_row,p_col,1+ntab+1*p_col);
hold on;
plot(tmp_w_,max(min(ylim_),min(max(ylim_),z_we__(:,2:end))),'k','LineWidth',linewidth_sml);
plot(tmp_w_,max(min(ylim_),min(max(ylim_),z_we__(:,1+0)  )),'r','LineWidth',linewidth_big);
hold off;
title(str_z_w,'Interpreter','none');
xlim([0,1]); xlabel('fraction nearest');
ylim(ylim_); ylabel('z score');
set(gca,'XTick',0:0.05:1.0); grid on; xtickangle(90);
%%%%;
subplot(p_row,p_col,1+ntab+2*p_col);
hold on;
plot(tmp_w_,prctile_we__(:,1+0)  ,'r','LineWidth',linewidth_big);
hold off;
title(str_p_w,'Interpreter','none');
xlim([0,1]); xlabel('fraction nearest');
ylim([0,1]); ylabel('p empirical');
set(gca,'YTick',0:0.05:1.0); set(gca,'XTick',0:0.05:1.0); grid on; xtickangle(90);
set(gca,'TickLength',[0,0]);
%%%%;
end;%for ntab=0:p_col-1;
%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Plots for manuscripts. Display only left column with altered style;
ntab = 0;
tmp_x_ = linspace(0,1,n_x); tmp_y_ = linspace(0,1,n_y); tmp_z_ = linspace(0,1,n_z);
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,768]);; p_row = 1; p_col = 3;
linewidth_sml = 0.5; linewidth_big = 2.0; ylim_ = 10*[-1,+1];
%%%%;
tmp_w_ = tmp_z_; f_we__ = f_ze__; z_we__ = z_ze__; prctile_we__ = prctile_ze__; str_f_w = 'f_z_'; str_z_w = 'z_z_'; str_p_w = 'p_z_'; 
%%%%;
subplot(p_row,p_col,1+ntab+0);
hold on;
plot(tmp_w_(tmp_w_<=0.5),f_we__(tmp_w_<=0.5,2:end),'k','LineWidth',linewidth_sml,'Color',[0 0 0 0.1]);
plot(tmp_w_(tmp_w_<=0.5),f_we__(tmp_w_<=0.5,1+0)  ,'r','LineWidth',linewidth_big);
hold off;
% title(str_f_w,'Interpreter','none');
xlim([0,.50]); 
set(gca,'XTick',0:0.05:.50,'FontSize',14); grid on; xtickangle(90);ytickformat('%.2f'),xtickformat('%.2f')
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
set(gca,'XTick',0:0.05:.50,'FontSize',14); grid on; xtickangle(90);ytickformat('%.2f'),xtickformat('%.2f')
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
set(gca,'YTick',0:0.05:1.0,'FontSize',14); set(gca,'XTick',0:0.05:.50,'FontSize',14); grid on; xtickangle(90); ytickformat('%.2f'),xtickformat('%.2f')
set(gca,'TickLength',[0,0]);
ylim([0,1]); 
xlabel('fraction nearest','FontSize',18);
ylabel('1-p empirical','FontSize',18);

%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save out max z of replication and index ;
[max_z, max_idx] = max(z_we__(tmp_w_<=0.5,1+0))
% Get fraction of training set to use as nearest neighbors that ;
% produces the max replication z ;
max_f = tmp_w_(max_idx);
tmp_fname_max_f = sprintf('%s/%s_max_f.txt',tmp_dir,str_datafile);
% write max_f to file
tmp_fid = fopen(tmp_fname_max_f,'w');
fprintf(tmp_fid,'%f',max_f);
fclose(tmp_fid);
