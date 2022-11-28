%%%%%%%%;
% renaming dir_Up05_maf01 as dir_Up05/dir_test2mds_maf01 ;
%%%%%%%%;

clear;
% platform = 'rusty';
% if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
% if (strcmp(platform,'access1')); str_home = 'data'; end;
% if (strcmp(platform,'OptiPlex')); str_home = 'home'; end;
% if (strcmp(platform,'eval1')); str_home = 'home'; end;
% if (strcmp(platform,'rusty')); str_home = 'mnt/home'; end;
%%%%%%%%;
run('/home/jelman/Github/lakcluster_c/dir_m/setup_0');
verbose = 1;
flag_disp = 1+verbose; nf=0;
flag_replot = 0;
if (verbose); disp(sprintf(' %% ;')); end;
if (verbose); disp(sprintf(' %% Assume path is set using dir_lakcluster_c/dir_m/setup_0.m. ;')); end;

dir_trunk = '/home/jelman/Projects/AD_Biclustering/data/UKB/ukb_geno_data';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

dir_0in = sprintf('%s/dir_Up05_maf01',dir_trunk);
if ~exist(dir_0in,'dir'); disp(sprintf(' %% Warning, %s not found',dir_0in)); return; end;
if  exist(dir_0in,'dir'); disp(sprintf(' %% %s found',dir_0in)); end;
dir_out = sprintf('%s/dir_Up05/dir_test2mds_maf01',dir_trunk);
if ~exist(dir_out,'dir'); disp(sprintf(' %% mkdir %s',dir_out)); mkdir(dir_out); end;
if  exist(dir_out,'dir'); disp(sprintf(' %% %s found',dir_out)); end;

filename_list_ = { ...
		   ,'Up05_maf01_A_full_n.b16' ...
		   ,'Up05_maf01_A_full_t.b16' ...
		   ,'Up05_maf01_A_01_n.b16' ...
		   ,'Up05_maf01_A_01_t.b16' ...
		   ,'Up05_maf01_bim.ext' ...
		   ,'Up05_maf01_fam.ext' ...
		   ,'Up05_maf01_mc_A.b16' ...
		   ,'Up05_maf01_mc_A_t05.b16' ...
		   ,'Up05_maf01_mc_A_t25.b16' ...
		   ,'Up05_maf01_mc_A_t50.b16' ...
		   ,'Up05_maf01_mc_T.b16' ...
		   ,'Up05_maf01_mc_T_m2r1.b16' ...
		   ,'Up05_maf01_mc_T_m2r2.b16' ...
		   ,'Up05_maf01_mc_T_m2r3.b16' ...
		   ,'Up05_maf01_mc_T_m2r4.b16' ...
		   ,'Up05_maf01_mr_A_01.b16' ...
		   ,'Up05_maf01_mr_A_continent0_01.b16' ...
		   ,'Up05_maf01_mr_A_continent0_full.b16' ...
		   ,'Up05_maf01_mr_A_continent1_01.b16' ...
		   ,'Up05_maf01_mr_A_continent1_full.b16' ...
		   ,'Up05_maf01_mr_A_continent2_01.b16' ...
		   ,'Up05_maf01_mr_A_continent2_full.b16' ...
		   ,'Up05_maf01_mr_A_full.b16' ...
		   ,'Up05_maf01_mr_Z_01.b16' ...
		   ,'Up05_maf01_mr_Z_continent0_01.b16' ...
		   ,'Up05_maf01_mr_Z_continent0_full.b16' ...
		   ,'Up05_maf01_mr_Z_continent1_01.b16' ...
		   ,'Up05_maf01_mr_Z_continent1_full.b16' ...
		   ,'Up05_maf01_mr_Z_continent2_01.b16' ...
		   ,'Up05_maf01_mr_Z_continent2_full.b16' ...
		   ,'Up05_maf01_mr_Z_full.b16' ...
		   ,'Up05_maf01_T_01_n.b16' ...
		   ,'Up05_maf01_T_01_t.b16' ...
		   ,'Up05_maf01_T_full_n.b16' ...
		   ,'Up05_maf01_T_full_t.b16' ...
		   ,'Up05_maf01_T_m2r1_01_n.b16' ...
		   ,'Up05_maf01_T_m2r1_01_t.b16' ...
		   ,'Up05_maf01_T_m2r1_full_n.b16' ...
		   ,'Up05_maf01_T_m2r1_full_n.mat' ...
		   ,'Up05_maf01_T_m2r1_full_t.b16' ...
		   ,'Up05_maf01_T_m2r1_kappa_squared.txt' ...
		   ,'Up05_maf01_T_m2r1_xx_n.mat' ...
		   ,'Up05_maf01_T_m2r2_01_n.b16' ...
		   ,'Up05_maf01_T_m2r2_01_t.b16' ...
		   ,'Up05_maf01_T_m2r2_full_n.b16' ...
		   ,'Up05_maf01_T_m2r2_full_n.mat' ...
		   ,'Up05_maf01_T_m2r2_full_t.b16' ...
		   ,'Up05_maf01_T_m2r2_kappa_squared.txt' ...
		   ,'Up05_maf01_T_m2r2_xx_n.mat' ...
		   ,'Up05_maf01_T_m2r3_01_n.b16' ...
		   ,'Up05_maf01_T_m2r3_01_t.b16' ...
		   ,'Up05_maf01_T_m2r3_full_n.b16' ...
		   ,'Up05_maf01_T_m2r3_full_n.mat' ...
		   ,'Up05_maf01_T_m2r3_full_t.b16' ...
		   ,'Up05_maf01_T_m2r3_kappa_squared.txt' ...
		   ,'Up05_maf01_T_m2r3_xx_n.mat' ...
		   ,'Up05_maf01_T_m2r4_01_n.b16' ...
		   ,'Up05_maf01_T_m2r4_01_t.b16' ...
		   ,'Up05_maf01_T_m2r4_full_n.b16' ...
		   ,'Up05_maf01_T_m2r4_full_n.mat' ...
		   ,'Up05_maf01_T_m2r4_full_t.b16' ...
		   ,'Up05_maf01_T_m2r4_kappa_squared.txt' ...
		   ,'Up05_maf01_T_m2r4_xx_n.mat' ...
};

n_l = numel(filename_list_);

for nl=0:n_l-1;
tmp_filename_0in = filename_list_{1+nl};
fname_0in = sprintf('%s/%s',dir_0in,tmp_filename_0in);
tmp_ij = strfind(tmp_filename_0in,'Up05');
tmp_filename_out = sprintf('%s%s','test2mds',tmp_filename_0in(tmp_ij+4:end));
fname_out = sprintf('%s/%s',dir_out,tmp_filename_out);
if  exist(fname_0in,'file') & ~exist(fname_out,'file');
str_command = sprintf('mv %s %s',fname_0in,fname_out);
disp(sprintf('%s',str_command));
system(str_command);
end;%if  exist(fname_0in,'file') & ~exist(fname_out,'file');
end;%for nl=0:n_l-1;

fname_del = sprintf('%s/test2mds_maf01_A_01_n.b16',dir_out),;
fname_sou = sprintf('%s/test2mds_maf01_A_full_n.b16',dir_out);
if  exist(fname_del,'file') &  exist(fname_sou,'file');
str_command = sprintf(' rm -f %s',fname_del);
disp(sprintf('%s',str_command));
system(str_command);
str_command = sprintf(' ln -s %s %s',fname_sou,fname_del);
disp(sprintf('%s',str_command));
system(str_command);
end;%if  exist(fname_del,'file') &  exist(fname_sou,'file');

fname_del = sprintf('%s/test2mds_maf01_A_01_t.b16',dir_out);
fname_sou = sprintf('%s/test2mds_maf01_A_full_t.b16',dir_out);
if  exist(fname_del,'file') &  exist(fname_sou,'file');
str_command = sprintf(' rm -f %s',fname_del);
disp(sprintf('%s',str_command));
system(str_command);
str_command = sprintf(' ln -s %s %s',fname_sou,fname_del);
disp(sprintf('%s',str_command));
system(str_command);
end;%if  exist(fname_del,'file') &  exist(fname_sou,'file');





