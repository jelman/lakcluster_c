clear;

%%%%%%%%;
run('/home/jelman/Github/lakcluster_c/dir_m/setup_0'); %<-- set up the paths. ;
flag_verbose = 1;
flag_disp = 1+flag_verbose; nf=0;
flag_replot = 0;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Assume path is set using dir_lakcluster_c/dir_m/setup_0.m. ;')); end;

if (flag_verbose); disp(sprintf(' %% Comparing Up05 with AFULLp05 data. ;')); end;

if (flag_verbose); disp(sprintf(' %% Comparing Up99 with Up05 data. ;')); end;
memory_GB = 128; %<-- maybe this should be increased? ;
if (flag_verbose); disp(sprintf(' %% trying with memory_GB %d ;',memory_GB)); end;

dir_code = '/home/jelman/Github/lakcluster_c';
dir_trunk = '/home/jelman/Projects/AD_Biclustering/data/ADNI/ADNI_vs_UKB';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

dir_mat_replication = sprintf('%s/dir_mat_replication',dir_trunk);
if ~exist(dir_mat_replication,'dir'); disp(sprintf(' %% mkdir %s',dir_mat_replication')); mkdir(dir_mat_replication); end;
dir_jpg_replication = sprintf('%s/dir_jpg_replication',dir_trunk);
if ~exist(dir_jpg_replication,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg_replication')); mkdir(dir_jpg_replication); end;


%%%%%%%%;
% First running test_Up05_vs_Ap05_18c_helper_0.m ;
%%%%%%%%;

%%%%%%%%;
% Now load one of the data-files containing the projected values and labels. ;
%%%%%%%%;
ncontinent = 1;
dataset = 'AFULLp05';
str_datafile = sprintf('trnUp05_tst_%s_ncontinent%d',dataset,ncontinent);

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
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==1);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[0.0,1.0,1.0]);
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==2);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[0.5,0.0,0.5]);
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==3);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[1.0,0.0,1.0]);
legend({'Up05 ctrl','Up05 case','Up05 bicl'},'Location','NorthWest');
grid on;
set(gca,'Fontsize',fontsize_use);
hold off;
subplot(1,2,2);
hold on;
tmp_index_ = efind(tmp_.mr_dvx_AFULLp05_from_trnUp05_tstAFULLp05_nix_==1);
X_xd__ = tmp_.AZnV_D_AFULLp05_from_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
plot(tmp_X_xd__(:,1+0),tmp_X_xd__(:,1+1),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.0,0.8,0.8]);
tmp_index_ = efind(tmp_.mr_dvx_AFULLp05_from_trnUp05_tstAFULLp05_nix_==2);
X_xd__ = tmp_.AZnV_D_AFULLp05_from_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
plot(tmp_X_xd__(:,1+0),tmp_X_xd__(:,1+1),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.8,0.0,0.8]);
legend({'AFULLp05 ctrl','AFULLp05 case'},'Location','NorthWest');
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
% Now we try and find a mapping from loaded AFULLp05 to loaded Up05. ;
% For AFULLp05 we will use the projected values: ;
% tmp_.AZnV_D_AFULLp05_from_trnUp05_tstAFULLp05_nix_p01_. ;
% For Up05 we will use the projected values: ;
% tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_. ;
% Note that affine_point_match_1 is designed to account for the labels in the data. ;
% For this point-cloud matching we will use only case vs ctrl labels (i.e., no bicluster labels). ;
%%%%%%%%;
tmp_index_Up05_ = efind( (tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==1) | (tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==2) | (tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==3) );
tmp_.Y_dy__ = transpose(tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_Up05_,:));
tmp_.label_Y_y_ = min(2,tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_(1+tmp_index_Up05_));
tmp_index_AFULLp05_ = efind( (tmp_.mr_dvx_AFULLp05_from_trnUp05_tstAFULLp05_nix_==1) | (tmp_.mr_dvx_AFULLp05_from_trnUp05_tstAFULLp05_nix_==2) );
tmp_.X_dx__ = transpose(tmp_.AZnV_D_AFULLp05_from_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_AFULLp05_,:));
tmp_.label_X_x_ = tmp_.mr_dvx_AFULLp05_from_trnUp05_tstAFULLp05_nix_(1+tmp_index_AFULLp05_);
parameter_apm = struct('type','parameter_apm');
parameter_apm.k_use = 32; %<-- This is the initial number of nearest neighbors to use. ;
parameter_apm.k_gamma = 0.50; %<-- Continent1 converges at .50, continent2 must be lowered to 0.15 ;
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
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==1);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[0.0,1.0,1.0]);
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==2);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[0.5,0.0,0.5]);
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==3);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
plot(Y_yd__(:,1+0),Y_yd__(:,1+1),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[1.0,0.0,1.0]);
tmp_index_ = efind(tmp_.mr_dvx_AFULLp05_from_trnUp05_tstAFULLp05_nix_==1);
X_xd__ = tmp_.AZnV_D_AFULLp05_from_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
plot(tmp_X_xd__(:,1+0),tmp_X_xd__(:,1+1),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.0,0.8,0.8]);
tmp_index_ = efind(tmp_.mr_dvx_AFULLp05_from_trnUp05_tstAFULLp05_nix_==2);
X_xd__ = tmp_.AZnV_D_AFULLp05_from_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
plot(tmp_X_xd__(:,1+0),tmp_X_xd__(:,1+1),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.8,0.0,0.8]);
legend({'Up05 ctrl','Up05 case','Up05 bicl','AFULLp05 ctrl','AFULLp05 case'},'Location','NorthWest');
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
% 2. the points (tmp_X_xd__) from the testing data (AFULLp05) ;
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
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==1);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
Up05_stack_yd__ = [Up05_stack_yd__;Y_yd__];
Up05_stack_m_ = [Up05_stack_m_;1*ones(size(Y_yd__,1),1)]; %<-- ctrl label. ;
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==2);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
Up05_stack_yd__ = [Up05_stack_yd__;Y_yd__];
Up05_stack_m_ = [Up05_stack_m_;2*ones(size(Y_yd__,1),1)]; %<-- case label. ;
tmp_index_ = efind(tmp_.mr_dvx_trnUp05_tstAFULLp05_nix_==3);
Y_yd__ = tmp_.AZnV_D_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
Up05_stack_yd__ = [Up05_stack_yd__;Y_yd__];
Up05_stack_m_ = [Up05_stack_m_;3*ones(size(Y_yd__,1),1)]; %<-- bicl label. ;
AFULLp05_stack_xd__ = zeros(0,2);
AFULLp05_stack_m_ = zeros(0,1);
tmp_index_ = efind(tmp_.mr_dvx_AFULLp05_from_trnUp05_tstAFULLp05_nix_==1);
X_xd__ = tmp_.AZnV_D_AFULLp05_from_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
AFULLp05_stack_xd__ = [AFULLp05_stack_xd__;tmp_X_xd__];
AFULLp05_stack_m_ = [AFULLp05_stack_m_;1*ones(size(tmp_X_xd__,1),1)]; %<-- ctrl label. ;
tmp_index_ = efind(tmp_.mr_dvx_AFULLp05_from_trnUp05_tstAFULLp05_nix_==2);
X_xd__ = tmp_.AZnV_D_AFULLp05_from_trnUp05_tstAFULLp05_nix_p01_(1+tmp_index_,:);
tmp_X_xd__ = transpose(a_est_ + A_est__*transpose(X_xd__));
AFULLp05_stack_xd__ = [AFULLp05_stack_xd__;tmp_X_xd__];
AFULLp05_stack_m_ = [AFULLp05_stack_m_;2*ones(size(tmp_X_xd__,1),1)]; %<-- case label. ;
%%%%%%%%;
[ ...
 parameter_apm ...
,vlXY_uxy___ ...
,vcXY_uxy___ ...
,u_Y_m_u_ ...
] = ...
affine_point_match_vector_label_0( ...
 parameter_apm ...
,transpose(AFULLp05_stack_xd__) ...
,transpose(Up05_stack_yd__) ...
,Up05_stack_m_ ...
);
%%%%%%%%;
% Now we can convert the vector-labels vlXY_uxy___ into conditional-labels. ;
% vbXY_xk__ is an array of size n_x by n_y, where: ;
% vbXY_xk__(1+nx,1+ny) holds the relative (i.e., conditional) bicluster-label for point nx, ;
% corresponding to the point at AFULLp05_stack_xd__(1+nx,:). ;
% as estimated using k=1+ny nearest-neighbors from Up05_stack_yd__. ;
%%%%%%%%;
n_x = size(AFULLp05_stack_xd__,1);
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
for tmp_f = [0.01,0.05,0.10,0.20,.30,.40];
nk = max(0,min(n_y-1,floor(tmp_f*n_y)));
subplot(2,3,1+np);np=np+1;
hold on;
tmp_index_ = efind(Up05_stack_m_==0+1); %<-- ctrl label. ;
plot(Up05_stack_yd__(1+tmp_index_,1+0),Up05_stack_yd__(1+tmp_index_,1+1),s0_ctrl,'LineStyle','none','MarkerEdgeColor',e0_ctrl_,'MarkerFaceColor',c0_ctrl_,'MarkerSize',markersize_sml);
tmp_index_ = efind(Up05_stack_m_==1+1); %<-- case label. ;
plot(Up05_stack_yd__(1+tmp_index_,1+0),Up05_stack_yd__(1+tmp_index_,1+1),s0_case,'LineStyle','none','MarkerEdgeColor',e0_case_,'MarkerFaceColor',c0_case_,'MarkerSize',markersize_sml);
tmp_index_ = efind(Up05_stack_m_==2+1); %<-- bicl label. ;
plot(Up05_stack_yd__(1+tmp_index_,1+0),Up05_stack_yd__(1+tmp_index_,1+1),s0_bicl,'LineStyle','none','MarkerEdgeColor',e0_bicl_,'MarkerFaceColor',c0_bicl_,'MarkerSize',markersize_sml);
tmp_index_ = efind(AFULLp05_stack_m_==0+1); %<-- ctrl label. ;
plot(AFULLp05_stack_xd__(1+tmp_index_,1+0),AFULLp05_stack_xd__(1+tmp_index_,1+1),s1_ctrl,'LineStyle','none','MarkerEdgeColor',e1_ctrl_,'MarkerFaceColor',c1_ctrl_,'MarkerSize',markersize_med);
tmp_index_ = efind(AFULLp05_stack_m_==1+1); %<-- case label. ;
n_l = numel(tmp_index_);
for nl=0:n_l-1;
tmp_index = tmp_index_(1+nl);
nc_use = max(0,min(n_c_use-1,floor(n_c_use*vbXY_xk__(1+tmp_index,1+nk)/1.0)));
c_use_ = c_use__(1+nc_use,:);
plot(AFULLp05_stack_xd__(1+tmp_index,1+0),AFULLp05_stack_xd__(1+tmp_index,1+1),s1_case,'LineStyle','none','MarkerEdgeColor',c_use_.^2,'MarkerFaceColor',c_use_,'MarkerSize',markersize_med);
end;%for nl=0:n_l-1;
hold off;
set(gca,'FontSize',fontsize_use);
axisnotick();
title(sprintf('f %0.2f',tmp_f));
legend({'Up05 ctrl','Up05 case','Up05 bicl','AFULLp05 ctrl','AFULLp05 case'},'Location','NorthWest');
end;%for tmp_f = [0.05,0.10,0.15,0.20,.30,.40];
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
% Based on the optimal nearest-neighbor-fraction f identified from the ;
% ADNI replication run in test_Up05_vs_Ap05_18c_helper_0.m, ;
% assign soft labels to all ADNI participants. ;
%%%%%%%%;
% Get name of training dataset;
str_trn = sprintf('trnUp05_tst_Ap05_ncontinent%d',ncontinent);
tmp_fname_max_f = sprintf('%s/dir_%s/%s_max_f.txt',dir_mat_replication,str_trn,str_trn);
% Read max_f from file
fid = fopen(tmp_fname_max_f,'r');
max_f = fscanf(fid,'%f');
fclose(fid);

% Find index matching optimal fraction
tmp_y_ = linspace(0,1,n_y);
[~, max_f_idx] = min(abs(tmp_y_-max_f));

% Extract conditional bicluster weights from optimal fraction of nearest ;
% neighbors ;
biclust_labels = vbXY_xk__(:,max_f_idx); 

% Load fam file 
fname_famext = '/home/jelman/Projects/AD_Biclustering/data/ADNI/ADNI_vs_UKB/dir_AFULLp05/dir_test2mds_maf01/test2mds_maf01_fam.ext';
[ ...
,n_patient ...
,fam_fid_ ...
,fam_iid_ ...
,fam_yid_ ...
,fam_xid_ ...
,fam_sex_ ...
,fam_dvx_ ...
,fam_dir_ ...
,fam_fidandiid_ ...
,~ ...
  ] = ...
load_famext_ver1( ...
 fname_famext ...
);
% Get subject ids belonging to continent ;
fam_iid_cont = fam_iid_(tmp_index_AFULLp05_ +1);
% Open the file for writing
tmp_fname_biclust_labels = sprintf('%s/dir_%s/AFULLp05_from_Up05_cap_biclust_labels_continent%d_.txt',dir_trunk,dataset,ncontinent);
fid = fopen(tmp_fname_biclust_labels, 'w');
% Write the data to the file
for i = 1:length(fam_iid_cont)
    fprintf(fid, '%s %d\n', fam_iid_cont{i}, biclust_labels(i));
end
% Close the file
fclose(fid); 