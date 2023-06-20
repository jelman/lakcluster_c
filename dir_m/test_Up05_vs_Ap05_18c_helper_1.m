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
% Consider first running test_Up05_vs_Ap05_18c.m ;
%%%%%%%%;

%%%%%%%%;
% Now load one of the data-files containing the projected values and labels. ;
%%%%%%%%;
%str_datafile = sprintf('trnUp05_tst_Ap05_ncontinent0');
%%%%%%%%;
% Alternatively, we can load the data from jeremy sent on 20230114. ;
%%%%%%%%;
str_datafile = sprintf('Up05_vs_Ap05_replication_cont1');

fname_mat = sprintf('%s/%s.mat',dir_mat_replication,str_datafile);
tmp_ = load(fname_mat);
  
a_est_ = zeros(2,1);
A_est__ = eye(2,2);
%%%%%%%%;
% First we just plot the loaded data: ;
%%%%%%%%;
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
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
plot(tmp_X_xd__(:,1+0),tmp_X_xd__(:,1+1),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.0,0.8,0.8]);
tmp_index_ = efind(tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2);
X_xd__ = tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
plot(tmp_X_xd__(:,1+0),tmp_X_xd__(:,1+1),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.8,0.0,0.8]);
legend({'Ap05 ctrl','Ap05 case'},'Location','NorthWest');
grid on;
set(gca,'Fontsize',fontsize_use);
hold off;
sgtitle(sprintf('%s',str_datafile),'Interpreter','none');
fname_fig_pre = sprintf('%s/%s_FIGA',dir_jpg_replication,str_datafile);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');

%%%%%%%%;
% Now we try and find a mapping from loaded Ap05 to loaded Up05. ;
% For Ap05 we will use the projected values: ;
% tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_. ;
% For Up05 we will use the projected values: ;
% tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_. ;
% Note that affine_point_match_1 is designed to account for the labels in the data. ;
% For this point-cloud matching we will use only case vs ctrl labels (i.e., no bicluster labels). ;
%%%%%%%%;
tmp_index_Up05_ = efind( (tmp_.mr_dvx_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==2) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==3) );
tmp_.Y_dy__ = transpose(tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_Up05_,:));
tmp_.label_Y_y_ = min(2,tmp_.mr_dvx_trnUp05_tstAp05_nix_(1+tmp_index_Up05_));
tmp_index_Ap05_ = efind( (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2) );
tmp_.X_dx__ = transpose(tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_Ap05_,:));
tmp_.label_X_x_ = tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_(1+tmp_index_Ap05_);
parameter_apm = struct('type','parameter_apm');
parameter_apm.k_use = 32; %<-- This is the initial number of nearest neighbors to use. ;
parameter_apm.k_gamma = 0.50; %<-- robust from 0.50 to 0.05 ;
parameter_apm.flag_disp = 0;
[parameter_apm,a_est_,A_est__,tmp_X_dx__] = ...
affine_point_match_1(parameter_apm,tmp_.X_dx__,tmp_.Y_dy__,[],[],tmp_.label_X_x_,tmp_.label_Y_y_);

%%%%%%%%;
% Now we can plot the results. ;
%%%%%%%%;
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
fname_fig_pre = sprintf('%s/%s_FIGB',dir_jpg_replication,str_datafile);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');

%%%%%%%%;
% Now we can use: ;
% 1. the labels (ctrls,case,bicl) from the training data (Up05) ;
% 2. the points (tmp_X_xd__) from the testing data (Ap05) ;
% to calculate approximate labels for the testing data. ;
% These approximate labels will be expressed in the form: ;
% vlXY_uxy___: a double-array of size numel(unique(Y_m_))-by-size(X_dx__,2)-by-size(Y_dy__,2) ;
% where vlXY_uxy___(1+nu,1+nx,1+ny) is the strength of association of: ;
% label nu to point nx using ny nearest neighbors from Y. ;
% More specifically: ;
% Assume we define k = 1+ny ; %<-- this is the number of nearest-neighbors we use to define the vector_label. ;
% and set: vlXY_ux__ = vlXY_ux___(:,:,1+ny), ;
% and then define: u_Y_m_u_ = unique(Y_m_) ; %<-- this is a listing of the unique labels from Y. ;
% then vlXY_ux__(:,1+nx) ;
% will be the vector label (of length numel(u_Y_m_u_)) associated with point nx defined by X_dx__(:,1+nx). ;
% The value vlXY_ux__(1+nu,1+nx) ;
% will be the value of this vector_label attributed to category u_Y_m_u_(1+nu). ;
% The array vcXY_uxy___ holds the associated counts. ;
% Note that the putative labels of X are not involved in this process. ;
%%%%%%%%;
Up05_stack_yd__ = zeros(0,2);
Up05_stack_m_ = zeros(0,1);
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==1);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
Up05_stack_yd__ = [Up05_stack_yd__;Y_yd__];
Up05_stack_m_ = [Up05_stack_m_;1*ones(size(Y_yd__,1),1)]; %<-- ctrl label. ;
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==2);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
Up05_stack_yd__ = [Up05_stack_yd__;Y_yd__];
Up05_stack_m_ = [Up05_stack_m_;2*ones(size(Y_yd__,1),1)]; %<-- case label. ;
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==3);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
Up05_stack_yd__ = [Up05_stack_yd__;Y_yd__];
Up05_stack_m_ = [Up05_stack_m_;3*ones(size(Y_yd__,1),1)]; %<-- bicl label. ;
Ap05_stack_xd__ = zeros(0,2);
Ap05_stack_m_ = zeros(0,1);
tmp_index_ = efind(tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1);
X_xd__ = tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
Ap05_stack_xd__ = [Ap05_stack_xd__;tmp_X_xd__];
Ap05_stack_m_ = [Ap05_stack_m_;1*ones(size(tmp_X_xd__,1),1)]; %<-- ctrl label. ;
tmp_index_ = efind(tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2);
X_xd__ = tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
Ap05_stack_xd__ = [Ap05_stack_xd__;tmp_X_xd__];
Ap05_stack_m_ = [Ap05_stack_m_;2*ones(size(tmp_X_xd__,1),1)]; %<-- case label. ;
%%%%%%%%;
[ ...
 parameter ...
,vlXY_uxy___ ...
,vcXY_uxy___ ...
,u_Y_m_u_ ...
] = ...
affine_point_match_vector_label_0( ...
 parameter ...
,transpose(Ap05_stack_xd__) ...
,transpose(Up05_stack_yd__) ...
,Up05_stack_m_ ...
);
%%%%%%%%;
% Now we can convert the vector-labels vlXY_uxy___ into conditional-labels. ;
% vbXY_xk__ is an array of size n_x by n_y, where: ;
% vbXY_xk__(1+nx,1+ny) holds the relative (i.e., conditional) bicluster-label for point nx, ;
% corresponding to the point at Ap05_stack_xd__(1+nx,:). ;
%%%%%%%%;
n_x = size(Ap05_stack_xd__,1);
n_y = size(Up05_stack_yd__,1);
vbXY_xk__ = zeros(n_x,n_y);
nu_bicl = 3-1; %<-- bicluster case-label index within u_Y_m_u_. ;
nu_case = 2-1; %<-- non-bicluster case-label index within u_Y_m_u_. ;
vbXY_xk__ = squeeze(vcXY_uxy___(1+nu_bicl,:,:)) ./ max(1,squeeze(sum(vcXY_uxy___(1+[nu_bicl,nu_case],:,:),1))) ;

%%%%%%%%;
% Now plot the results across a range of nearest-neighbor-fraction f. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
markersize_big = 12; markersize_med = 10; markersize_sml = 8; fontsize_use = 12;
s0_ctrl = '^'; c0_ctrl_ = 0.975*[1,1,1]; e0_ctrl_ = 0.975^2*[1,1,1];
s0_case = 's'; c0_case_ = 0.950*[1,1,1]; e0_case_ = 0.950^2*[1,1,1];
s0_bicl = 'o'; c0_bicl_ = 0.925*[1,1,1]; e0_bicl_ = 0.925^2*[1,1,1];
s1_ctrl = '^'; c1_ctrl_ = [1.0,0.9,0.8]; e1_ctrl_ = [1.0,0.9,0.8].^2;
s1_case = 'o'; 
c_use__ = colormap_80s; n_c_use = size(c_use__,1);
np=0;
for tmp_f = [0.01,0.05,0.10,0.20];
nk = max(0,min(n_y-1,floor(tmp_f*n_y)));
subplot(2,2,1+np);np=np+1;
hold on;
tmp_index_ = efind(Up05_stack_m_==0+1); %<-- ctrl label. ;
plot(Up05_stack_yd__(1+tmp_index_,1+0),Up05_stack_yd__(1+tmp_index_,1+1),s0_ctrl,'LineStyle','none','MarkerEdgeColor',e0_ctrl_,'MarkerFaceColor',c0_ctrl_,'MarkerSize',markersize_sml);
tmp_index_ = efind(Up05_stack_m_==1+1); %<-- case label. ;
plot(Up05_stack_yd__(1+tmp_index_,1+0),Up05_stack_yd__(1+tmp_index_,1+1),s0_case,'LineStyle','none','MarkerEdgeColor',e0_case_,'MarkerFaceColor',c0_case_,'MarkerSize',markersize_sml);
tmp_index_ = efind(Up05_stack_m_==2+1); %<-- bicl label. ;
plot(Up05_stack_yd__(1+tmp_index_,1+0),Up05_stack_yd__(1+tmp_index_,1+1),s0_bicl,'LineStyle','none','MarkerEdgeColor',e0_bicl_,'MarkerFaceColor',c0_bicl_,'MarkerSize',markersize_sml);
tmp_index_ = efind(Ap05_stack_m_==0+1); %<-- ctrl label. ;
plot(Ap05_stack_xd__(1+tmp_index_,1+0),Ap05_stack_xd__(1+tmp_index_,1+1),s1_ctrl,'LineStyle','none','MarkerEdgeColor',e1_ctrl_,'MarkerFaceColor',c1_ctrl_,'MarkerSize',markersize_med);
tmp_index_ = efind(Ap05_stack_m_==1+1); %<-- case label. ;
n_l = numel(tmp_index_);
for nl=0:n_l-1;
tmp_index = tmp_index_(1+nl);
nc_use = max(0,min(n_c_use-1,floor(n_c_use*vbXY_xk__(1+tmp_index,1+nk)/1.0)));
c_use_ = c_use__(1+nc_use,:);
plot(Ap05_stack_xd__(1+tmp_index,1+0),Ap05_stack_xd__(1+tmp_index,1+1),s1_case,'LineStyle','none','MarkerEdgeColor',c_use_.^2,'MarkerFaceColor',c_use_,'MarkerSize',markersize_med);
end;%for nl=0:n_l-1;
hold off;
set(gca,'FontSize',fontsize_use);
axisnotick();
title(sprintf('f %0.2f',tmp_f));
legend({'Up05 ctrl','Up05 case','Up05 bicl','Ap05 ctrl','Ap05 case'},'Location','NorthWest');
end;%for tmp_f = [0.05,0.10,0.15,0.20];
%%%%;
fname_fig_pre = sprintf('%s/%s_vector_label_FIGD',dir_jpg_replication,str_datafile);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;

%%%%%%%%;
% Now estimate p-value. ;
% Note that this p-value estimate assumes that the affine-point-transformation is fixed, ;
% and calculates the p-value using 'interior' shuffles only ;
% (i.e., permuting case-ctrl labels *after* the affine-point-transformation is determined). ;
%%%%%%%%;
Up05_stack_yd__ = zeros(0,2);
Up05_stack_m_ = zeros(0,1);
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==1);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
Up05_stack_yd__ = [Up05_stack_yd__;Y_yd__];
Up05_stack_m_ = [Up05_stack_m_;1*ones(size(Y_yd__,1),1)]; %<-- ctrl label. ;
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==2);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
Up05_stack_yd__ = [Up05_stack_yd__;Y_yd__];
Up05_stack_m_ = [Up05_stack_m_;2*ones(size(Y_yd__,1),1)]; %<-- case label. ;
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAp05_nix_==3);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
Up05_stack_yd__ = [Up05_stack_yd__;Y_yd__];
Up05_stack_m_ = [Up05_stack_m_;2*ones(size(Y_yd__,1),1)]; %<-- case label, and not bicl label just yet. ;
Ap05_stack_xd__ = zeros(0,2);
Ap05_stack_m_ = zeros(0,1);
tmp_index_ = efind(tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1);
X_xd__ = tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
Ap05_stack_xd__ = [Ap05_stack_xd__;tmp_X_xd__];
Ap05_stack_m_ = [Ap05_stack_m_;1*ones(size(tmp_X_xd__,1),1)]; %<-- ctrl label. ;
tmp_index_ = efind(tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2);
X_xd__ = tmp_.AZnV_D_Ap05_from_trnUp05_tstAp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
Ap05_stack_xd__ = [Ap05_stack_xd__;tmp_X_xd__];
Ap05_stack_m_ = [Ap05_stack_m_;2*ones(size(tmp_X_xd__,1),1)]; %<-- case label. ;
%%%%%%%%;
parameter_apm.n_shuffle = 32;
[ ...
 parameter_apm ...
,z_z_ ...
,z_zp__ ...
,z_x_ ...
,z_xp__ ...
,z_y_ ...
,z_yp__ ...
,f_z_ ...
,f_zp__ ...
,f_x_ ...
,f_xp__ ...
,f_y_ ...
,f_yp__ ...
] = ...
affine_point_match_p_0( ...
 parameter_apm ...
,transpose(Ap05_stack_xd__) ...
,Ap05_stack_m_ ...
,transpose(Up05_stack_yd__) ...
,Up05_stack_m_ ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_y = size(Up05_stack_yd__,1); n_x = size(Ap05_stack_xd__,1); n_z = max(n_x,n_y);
tmp_x_ = linspace(0,1,n_x); tmp_y_ = linspace(0,1,n_y); tmp_z_ = linspace(0,1,n_z);
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig; p_row = 3; p_col = 3;
linewidth_sml = 0.5; linewidth_big = 2.0; ylim_ = 10*[-1,+1];
%%%%;
for ntab=0:p_col-1;
if ntab==0; tmp_w_ = tmp_z_; f_w_ = f_z_; f_wp__ = f_zp__; z_w_ = z_z_; z_wp__ = z_zp__; str_f_w = 'f_z_'; str_z_w = 'z_z_'; str_p_w = 'p_z_'; end;
if ntab==1; tmp_w_ = tmp_x_; f_w_ = f_x_; f_wp__ = f_xp__; z_w_ = z_x_; z_wp__ = z_xp__; str_f_w = 'f_x_'; str_z_w = 'z_x_'; str_p_w = 'p_x_'; end;
if ntab==2; tmp_w_ = tmp_y_; f_w_ = f_y_; f_wp__ = f_yp__; z_w_ = z_y_; z_wp__ = z_yp__; str_f_w = 'f_y_'; str_z_w = 'z_y_'; str_p_w = 'p_y_'; end;
%%%%;
subplot(p_row,p_col,1+ntab+0*p_col);
hold on;
plot(tmp_w_,f_wp__,'k','LineWidth',linewidth_sml);
plot(tmp_w_,f_w_,'r','LineWidth',linewidth_big);
hold off;
title(str_f_w,'Interpreter','none');
xlim([0,1]); xlabel('fraction nearest');
ylabel('f fraction');
set(gca,'XTick',0:0.05:1.0); grid on;
%%%%;
subplot(p_row,p_col,1+ntab+1*p_col);
hold on;
plot(tmp_w_,max(min(ylim_),min(max(ylim_),z_wp__)),'k','LineWidth',linewidth_sml);
plot(tmp_w_,max(min(ylim_),min(max(ylim_),z_w_)),'r','LineWidth',linewidth_big);
hold off;
title(str_z_w,'Interpreter','none');
xlim([0,1]); xlabel('fraction nearest');
ylim(ylim_); ylabel('z score');
set(gca,'XTick',0:0.05:1.0); grid on;
%%%%;
subplot(p_row,p_col,1+ntab+2*p_col);
f_wp0__ = [f_w_ , f_wp__];
[~,tmp_ij_wp__] = sort(f_wp0__,2,'ascend');
[~,tmp_p_wp__] = sort(tmp_ij_wp__,2,'ascend');
tmp_p_wp__ = tmp_p_wp__/size(tmp_p_wp__,2);
hold on;
%plot(tmp_w_,tmp_p_wp__(:,2:end),'k','LineWidth',linewidth_sml);
plot(tmp_w_,tmp_p_wp__(:,1    ),'r','LineWidth',linewidth_big);
hold off;
title(str_p_w,'Interpreter','none');
xlim([0,1]); xlabel('fraction nearest');
ylim([0,1]); ylabel('p empirical');
set(gca,'XTick',0:0.05:1.0); grid on;
set(gca,'YTick',0:0.05:1.0); grid on;
set(gca,'TickLength',[0,0]);
%%%%;
end;%for ntab=0:p_col-1;
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGC',dir_jpg_replication,str_datafile);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now we estimate the p-value again, ;
% However, this time we calculate the p-value using a permutation-test ;
% where the label-shuffling is applied *before* the affine-transformation is determined. ;
% This time the number of 'exterior' shuffles denotes the number of trials, ;
% while the 'interior' shuffles are not used. ;
%%%%%%%%;
tmp_dir = sprintf('%s/dir_%s',dir_mat_replication,str_datafile);
if ~exist(tmp_dir,'dir'); disp(sprintf(' %% mkdir %s',tmp_dir)); mkdir(tmp_dir); end;
parameter_apm.n_shuffle_exterior = 128;
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
end;%for nshuffle_exterior=0:n_shuffle_exterior;

n_y = numel( efind( (tmp_.mr_dvx_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==2) | (tmp_.mr_dvx_trnUp05_tstAp05_nix_==3) ) );
n_x = numel( efind( (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==1) | (tmp_.mr_dvx_Ap05_from_trnUp05_tstAp05_nix_==2) ) );
n_z = max(n_x,n_y);
n_shuffle_exterior = parameter_apm.n_shuffle_exterior;
n_shuffle_interior = parameter_apm.n_shuffle;
f_zie___ = zeros(n_z,1+n_shuffle_interior,1+n_shuffle_exterior);
f_xie___ = zeros(n_x,1+n_shuffle_interior,1+n_shuffle_exterior);
f_yie___ = zeros(n_y,1+n_shuffle_interior,1+n_shuffle_exterior);
z_zie___ = zeros(n_z,1+n_shuffle_interior,1+n_shuffle_exterior);
z_xie___ = zeros(n_x,1+n_shuffle_interior,1+n_shuffle_exterior);
z_yie___ = zeros(n_y,1+n_shuffle_interior,1+n_shuffle_exterior);
for nshuffle_exterior=0:1+n_shuffle_exterior;
tmp_fname_pre = sprintf('%s/%s_s%.4d',tmp_dir,str_datafile,nshuffle_exterior);
tmp_fname_mat = sprintf('%s.mat',tmp_fname_pre);
if  exist(tmp_fname_mat,'file');
tmp_tmp_ = load(tmp_fname_mat);
f_zie___(:,1+0                     ,1+nshuffle_exterior) = tmp_tmp_.tmp_f_z_;
f_zie___(:,1+[1:n_shuffle_interior],1+nshuffle_exterior) = tmp_tmp_.tmp_f_zp__;
f_xie___(:,1+0                     ,1+nshuffle_exterior) = tmp_tmp_.tmp_f_x_;
f_xie___(:,1+[1:n_shuffle_interior],1+nshuffle_exterior) = tmp_tmp_.tmp_f_xp__;
f_yie___(:,1+0                     ,1+nshuffle_exterior) = tmp_tmp_.tmp_f_y_;
f_yie___(:,1+[1:n_shuffle_interior],1+nshuffle_exterior) = tmp_tmp_.tmp_f_yp__;
z_zie___(:,1+0                     ,1+nshuffle_exterior) = tmp_tmp_.tmp_z_z_;
z_zie___(:,1+[1:n_shuffle_interior],1+nshuffle_exterior) = tmp_tmp_.tmp_z_zp__;
z_xie___(:,1+0                     ,1+nshuffle_exterior) = tmp_tmp_.tmp_z_x_;
z_xie___(:,1+[1:n_shuffle_interior],1+nshuffle_exterior) = tmp_tmp_.tmp_z_xp__;
z_yie___(:,1+0                     ,1+nshuffle_exterior) = tmp_tmp_.tmp_z_y_;
z_yie___(:,1+[1:n_shuffle_interior],1+nshuffle_exterior) = tmp_tmp_.tmp_z_yp__;
clear tmp_tmp_;
end;%if  exist(tmp_fname_mat,'file');
end;%for nshuffle_exterior=0:n_shuffle_exterior;

%{
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 5; p_col = 7;
for nshuffle_exterior = 0:min(p_row*p_col-1,n_shuffle_exterior);
subplot(p_row,p_col,1+nshuffle_exterior);
hold on; plot(f_zie___(:,2:end,1+nshuffle_exterior),'k-'); plot(f_zie___(:,1+0,1+nshuffle_exterior),'r-'); hold off;
title(sprintf('%d',nshuffle_exterior));
end;%for nshuffle_exterior = 0:min(p_row*p_col-1,n_shuffle_exterior);
 %}

%%%%%%%%;
% Now we can plot the trace associated with this second p-value estimate. ;
% For now we just plot the raw trace (as calculated using z-scores). ;
% One can transform this into a plot based on empirical-p-values (i.e., ranks or percentiles). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
f_z0e__ = squeeze(f_zie___(:,1+0,:));
f_x0e__ = squeeze(f_xie___(:,1+0,:));
f_y0e__ = squeeze(f_yie___(:,1+0,:));
n_z = size(f_z0e__,1); tmp_z_ = linspace(0,1,n_z);
hold on; plot(tmp_z_,f_z0e__(:,2:end),'k-'); plot(tmp_z_,f_z0e__(:,1+0),'r-'); hold off;
