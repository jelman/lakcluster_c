verbose=1;
if (verbose); disp(sprintf(' %% [entering test_xxxcluster_fromdisk_ver16]')); end;
rng(0); nf=0;

if (verbose); disp(sprintf(' %% ')); end;
if (verbose); disp(sprintf(' %% First we generate a very simple example of dosage-data ;')); end;
if (verbose); disp(sprintf(' %% which exhibits a case-specific bicluster bcA, ;')); end;
if (verbose); disp(sprintf(' %% as well as a disjoint (larger) non-specific bicluster bcB. ;')); end;
str_lak_vs_dex = 'lak';
r_fraction_bcA = 0.125;
r_fraction_bcB = 0.250;
c_fraction_bcA = 0.0625;
c_fraction_bcB = 0.1250;
maf_lob = 0.25;
maf_upb = 0.75;

%%%%%%%%;
n_snp_cup = 1024*8;
bim_cup_khr_ = zeros(n_snp_cup,1);
bim_cup_vid_ = cell(n_snp_cup,1);
bim_cup_pos_ = zeros(n_snp_cup,1);
bim_cup_bpc_ = zeros(n_snp_cup,1);
bim_cup_al1_ = cell(n_snp_cup,1);
bim_cup_al2_ = cell(n_snp_cup,1);
tmp_al1_ = {'A','T','G','C'};
tmp_al2_ = {'G','C','A','T'};
for nsnp_cup=0:n_snp_cup-1;
bim_cup_khr_(1+nsnp_cup) = 1+mod(nsnp_cup,22);
bim_cup_vid_{1+nsnp_cup} = sprintf('vid%.4d',1+nsnp_cup);
bim_cup_pos_(1+nsnp_cup) = 1+nsnp_cup;
bim_cup_bpc_(1+nsnp_cup) = 1+nsnp_cup;
bim_cup_al1_{1+nsnp_cup} = sprintf('%s',tmp_al1_{1+mod(nsnp_cup,4)});
bim_cup_al2_{1+nsnp_cup} = sprintf('%s',tmp_al2_{1+mod(nsnp_cup,4)});
end;%for nsnp_cup=0:n_snp_cup-1;
cij_ = transpose(1:n_snp_cup);
cij_bcA_ = cij_(randperm(n_snp_cup,floor(n_snp_cup*c_fraction_bcA)));
n_snp_bcA = numel(cij_bcA_);
cij_tmp_ = setdiff(cij_,cij_bcA_);
cij_bcB_ = cij_tmp_(randperm(numel(cij_tmp_),floor(n_snp_cup*c_fraction_bcB)));
n_snp_bcB = numel(cij_bcB_);
%%%%%%%%;
n_patient_cup = 1024*4;
fam_cup_fid_ = cell(n_patient_cup,1);
fam_cup_iid_ = cell(n_patient_cup,1);
fam_cup_yid_ = cell(n_patient_cup,1);
fam_cup_xid_ = cell(n_patient_cup,1);
fam_cup_sex_ = zeros(n_patient_cup,1);
fam_cup_dvx_ = zeros(n_patient_cup,1);
fam_cup_fidandiid_ = cell(n_patient_cup,1);
tmp_sex_ = [1,2,1,2]; %<-- 1==male, 2==fema, ignored. ;
tmp_dvx_ = [1,1,2,2]; %<-- 1==ctrl, 2==case, not ignored. ;
for npatient_cup=0:n_patient_cup-1;
fam_cup_fid_{1+npatient_cup} = sprintf('fid%.4d',1+npatient_cup);
fam_cup_iid_{1+npatient_cup} = sprintf('iid%.4d',1+npatient_cup);
fam_cup_yid_{1+npatient_cup} = sprintf('yid%.4d',1+npatient_cup);
fam_cup_xid_{1+npatient_cup} = sprintf('xid%.4d',1+npatient_cup);
fam_cup_sex_(1+npatient_cup) = tmp_sex_(1+mod(npatient_cup,4));
fam_cup_dvx_(1+npatient_cup) = tmp_dvx_(1+floor(4*npatient_cup/n_patient_cup)); %<-- first half of patients are cases. ;
end;%for npatient_cup=0:n_patient_cup-1;
for npatient_cup=0:n_patient_cup-1;
fam_cup_fidandiid_{1+npatient_cup} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for npatient_cup=0:n_patient_cup-1;
rij_case_ = find(fam_cup_dvx_==2); rij_ctrl_ = find(fam_cup_dvx_==1); 
n_patient_case = numel(rij_case_); n_patient_ctrl = numel(rij_ctrl_);
rij_case_bcA_ = rij_case_(randperm(n_patient_case,floor(r_fraction_bcA*n_patient_case)));
rij_ctrl_bcA_ = [];
n_patient_case_bcA = numel(rij_case_bcA_);
n_patient_ctrl_bcA = numel(rij_ctrl_bcA_);
rij_case_tmp_ = setdiff(rij_case_,rij_case_bcA_);
rij_ctrl_tmp_ = setdiff(rij_ctrl_,rij_ctrl_bcA_);
rij_case_bcB_ = rij_case_tmp_(randperm(numel(rij_case_tmp_),floor(r_fraction_bcB*n_patient_case)));
rij_ctrl_bcB_ = rij_ctrl_tmp_(randperm(numel(rij_ctrl_tmp_),floor(r_fraction_bcB*n_patient_ctrl)));
n_patient_case_bcB = numel(rij_case_bcB_);
n_patient_ctrl_bcB = numel(rij_ctrl_bcB_);
rij_case_bcA_cmp_ = setdiff(rij_case_,rij_case_bcA_); n_patient_case_bcA_cmp = numel(rij_case_bcA_cmp_);
rij_case_bcB_cmp_ = setdiff(rij_case_,rij_case_bcB_); n_patient_case_bcB_cmp = numel(rij_case_bcB_cmp_);
rij_ctrl_bcB_cmp_ = setdiff(rij_ctrl_,rij_ctrl_bcB_); n_patient_ctrl_bcB_cmp = numel(rij_ctrl_bcB_cmp_);
%%%%%%%%;
frq_s_ = linspace(maf_lob,maf_upb,n_snp_cup);
ds1_cup_ps__ = zeros(n_patient_cup,n_snp_cup);
ds1_cup_ps__ = rand(n_patient_cup,n_snp_cup)<repmat(frq_s_,[n_patient_cup,1]);
ds2_cup_ps__ = zeros(n_patient_cup,n_snp_cup);
ds2_cup_ps__ = rand(n_patient_cup,n_snp_cup)<repmat(frq_s_,[n_patient_cup,1]);
if (strcmp(str_lak_vs_dex,'lak'));
%%%%;
% now correlate diplotypes within bcA. ;
%%%%;
ds1_cup_ps__(rij_case_bcA_,cij_bcA_) = repmat(transpose(linspace(0,1,n_patient_case_bcA)),[1,n_snp_bcA]) < repmat(frq_s_(cij_bcA_),[n_patient_case_bcA,1]);
ds1_cup_ps__(rij_case_bcA_,cij_bcA_) = ds1_cup_ps__(rij_case_bcA_(randperm(n_patient_case_bcA)),cij_bcA_);
ds2_cup_ps__(rij_case_bcA_,cij_bcA_) = repmat(transpose(linspace(0,1,n_patient_case_bcA)),[1,n_snp_bcA]) < repmat(frq_s_(cij_bcA_),[n_patient_case_bcA,1]);
ds2_cup_ps__(rij_case_bcA_,cij_bcA_) = ds2_cup_ps__(rij_case_bcA_(randperm(n_patient_case_bcA)),cij_bcA_);
%%%%;
% now correlate diplotypes within bcB. ;
%%%%;
ds1_cup_ps__(rij_case_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_case_bcB)),[1,n_snp_bcB]) < repmat(frq_s_(cij_bcB_),[n_patient_case_bcB,1]);
ds1_cup_ps__(rij_case_bcB_,cij_bcB_) = ds1_cup_ps__(rij_case_bcB_(randperm(n_patient_case_bcB)),cij_bcB_);
ds2_cup_ps__(rij_case_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_case_bcB)),[1,n_snp_bcB]) < repmat(frq_s_(cij_bcB_),[n_patient_case_bcB,1]);
ds2_cup_ps__(rij_case_bcB_,cij_bcB_) = ds2_cup_ps__(rij_case_bcB_(randperm(n_patient_case_bcB)),cij_bcB_);
ds1_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_ctrl_bcB)),[1,n_snp_bcB]) < repmat(frq_s_(cij_bcB_),[n_patient_ctrl_bcB,1]);
ds1_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = ds1_cup_ps__(rij_ctrl_bcB_(randperm(n_patient_ctrl_bcB)),cij_bcB_);
ds2_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_ctrl_bcB)),[1,n_snp_bcB]) < repmat(frq_s_(cij_bcB_),[n_patient_ctrl_bcB,1]);
ds2_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = ds2_cup_ps__(rij_ctrl_bcB_(randperm(n_patient_ctrl_bcB)),cij_bcB_);
%%%%;
end;%if (strcmp(str_lak_vs_dex,'lak'));
%%%%%%%%;
if (strcmp(str_lak_vs_dex,'dex'));
%%%%;
% now lift diplotypes within bcA. ;
%%%%;
frq_case_bcA_s_ = (frq_s_ >= 0.5) .* (frq_s_ + 0.5*(maf_upb-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + 0.5*(maf_lob-frq_s_)) ;
frq_case_bcA_cmp_s_ = (frq_s_*n_patient_case - frq_case_bcA_s_*n_patient_case_bcA)/(n_patient_case - n_patient_case_bcA);
ds1_cup_ps__(rij_case_bcA_,cij_bcA_) = repmat(transpose(linspace(0,1,n_patient_case_bcA)),[1,n_snp_bcA]) < repmat(frq_case_bcA_s_(cij_bcA_),[n_patient_case_bcA,1]);
ds1_cup_ps__(rij_case_bcA_,cij_bcA_) = ds1_cup_ps__(rij_case_bcA_(randperm(n_patient_case_bcA)),cij_bcA_);
ds2_cup_ps__(rij_case_bcA_,cij_bcA_) = repmat(transpose(linspace(0,1,n_patient_case_bcA)),[1,n_snp_bcA]) < repmat(frq_case_bcA_s_(cij_bcA_),[n_patient_case_bcA,1]);
ds2_cup_ps__(rij_case_bcA_,cij_bcA_) = ds2_cup_ps__(rij_case_bcA_(randperm(n_patient_case_bcA)),cij_bcA_);
ds1_cup_ps__(rij_case_bcA_cmp_,cij_bcA_) = rand(n_patient_case_bcA_cmp,n_snp_bcA) < repmat(frq_case_bcA_cmp_s_(cij_bcA_),[n_patient_case_bcA_cmp,1]);
ds1_cup_ps__(rij_case_bcA_cmp_,cij_bcA_) = ds1_cup_ps__(rij_case_bcA_cmp_(randperm(n_patient_case_bcA_cmp)),cij_bcA_);
ds2_cup_ps__(rij_case_bcA_cmp_,cij_bcA_) = rand(n_patient_case_bcA_cmp,n_snp_bcA) < repmat(frq_case_bcA_cmp_s_(cij_bcA_),[n_patient_case_bcA_cmp,1]);
ds2_cup_ps__(rij_case_bcA_cmp_,cij_bcA_) = ds2_cup_ps__(rij_case_bcA_cmp_(randperm(n_patient_case_bcA_cmp)),cij_bcA_);
%%%%;
% now lift diplotypes within bcB. ;
%%%%;
frq_case_bcB_s_ = (frq_s_ >= 0.5) .* (frq_s_ + 0.5*(maf_upb-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + 0.5*(maf_lob-frq_s_)) ;
frq_case_bcB_cmp_s_ = (frq_s_*n_patient_case - frq_case_bcB_s_*n_patient_case_bcB)/(n_patient_case - n_patient_case_bcB);
ds1_cup_ps__(rij_case_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_case_bcB)),[1,n_snp_bcB]) < repmat(frq_case_bcB_s_(cij_bcB_),[n_patient_case_bcB,1]);
ds1_cup_ps__(rij_case_bcB_,cij_bcB_) = ds1_cup_ps__(rij_case_bcB_(randperm(n_patient_case_bcB)),cij_bcB_);
ds2_cup_ps__(rij_case_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_case_bcB)),[1,n_snp_bcB]) < repmat(frq_case_bcB_s_(cij_bcB_),[n_patient_case_bcB,1]);
ds2_cup_ps__(rij_case_bcB_,cij_bcB_) = ds2_cup_ps__(rij_case_bcB_(randperm(n_patient_case_bcB)),cij_bcB_);
ds1_cup_ps__(rij_case_bcB_cmp_,cij_bcB_) = rand(n_patient_case_bcB_cmp,n_snp_bcB) < repmat(frq_case_bcB_cmp_s_(cij_bcB_),[n_patient_case_bcB_cmp,1]);
ds1_cup_ps__(rij_case_bcB_cmp_,cij_bcB_) = ds1_cup_ps__(rij_case_bcB_cmp_(randperm(n_patient_case_bcB_cmp)),cij_bcB_);
ds2_cup_ps__(rij_case_bcB_cmp_,cij_bcB_) = rand(n_patient_case_bcB_cmp,n_snp_bcB) < repmat(frq_case_bcB_cmp_s_(cij_bcB_),[n_patient_case_bcB_cmp,1]);
ds2_cup_ps__(rij_case_bcB_cmp_,cij_bcB_) = ds2_cup_ps__(rij_case_bcB_cmp_(randperm(n_patient_case_bcB_cmp)),cij_bcB_);
frq_ctrl_bcB_s_ = (frq_s_ >= 0.5) .* (frq_s_ + 0.5*(maf_upb-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + 0.5*(maf_lob-frq_s_)) ;
frq_ctrl_bcB_cmp_s_ = (frq_s_*n_patient_ctrl - frq_ctrl_bcB_s_*n_patient_ctrl_bcB)/(n_patient_ctrl - n_patient_ctrl_bcB);
ds1_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_ctrl_bcB)),[1,n_snp_bcB]) < repmat(frq_ctrl_bcB_s_(cij_bcB_),[n_patient_ctrl_bcB,1]);
ds1_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = ds1_cup_ps__(rij_ctrl_bcB_(randperm(n_patient_ctrl_bcB)),cij_bcB_);
ds2_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_ctrl_bcB)),[1,n_snp_bcB]) < repmat(frq_ctrl_bcB_s_(cij_bcB_),[n_patient_ctrl_bcB,1]);
ds2_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = ds2_cup_ps__(rij_ctrl_bcB_(randperm(n_patient_ctrl_bcB)),cij_bcB_);
ds1_cup_ps__(rij_ctrl_bcB_cmp_,cij_bcB_) = rand(n_patient_ctrl_bcB_cmp,n_snp_bcB) < repmat(frq_ctrl_bcB_cmp_s_(cij_bcB_),[n_patient_ctrl_bcB_cmp,1]);
ds1_cup_ps__(rij_ctrl_bcB_cmp_,cij_bcB_) = ds1_cup_ps__(rij_ctrl_bcB_cmp_(randperm(n_patient_ctrl_bcB_cmp)),cij_bcB_);
ds2_cup_ps__(rij_ctrl_bcB_cmp_,cij_bcB_) = rand(n_patient_ctrl_bcB_cmp,n_snp_bcB) < repmat(frq_ctrl_bcB_cmp_s_(cij_bcB_),[n_patient_ctrl_bcB_cmp,1]);
ds2_cup_ps__(rij_ctrl_bcB_cmp_,cij_bcB_) = ds2_cup_ps__(rij_ctrl_bcB_cmp_(randperm(n_patient_ctrl_bcB_cmp)),cij_bcB_);
%%%%;
end;%if (strcmp(str_lak_vs_dex,'dex'));
%%%%%%%%;
dsg_cup_ps__ = ds1_cup_ps__ + ds2_cup_ps__;
dsg_cup_frq_mss_ = mean(dsg_cup_ps__>=3,1); tmp_z_ = 1 - dsg_cup_frq_mss_;
dsg_cup_frq_nor_ = mean(dsg_cup_ps__==0,1)./tmp_z_;
dsg_cup_frq_xor_ = mean(dsg_cup_ps__==1,1)./tmp_z_;
dsg_cup_frq_and_ = mean(dsg_cup_ps__==2,1)./tmp_z_;
dsg_cup_p_opt_ = dsg_cup_frq_and_ + 0.5*dsg_cup_frq_xor_ ;
dsg_cup_q_opt_ = dsg_cup_frq_nor_ + 0.5*dsg_cup_frq_xor_ ;
tmp_and_ = dsg_cup_frq_and_.*log(dsg_cup_frq_and_./(dsg_cup_p_opt_.^2)); tmp_and_(find(~isfinite(tmp_and_)))=0;
tmp_xor_ = dsg_cup_frq_xor_.*log(dsg_cup_frq_xor_./(2*dsg_cup_p_opt_.*dsg_cup_q_opt_)); tmp_xor_(find(~isfinite(tmp_xor_)))=0;
tmp_nor_ = dsg_cup_frq_nor_.*log(dsg_cup_frq_nor_./(dsg_cup_q_opt_.^2)); tmp_nor_(find(~isfinite(tmp_nor_)))=0;
dsg_cup_I_opt_ = tmp_and_ + tmp_xor_ + tmp_nor_ ;
%%%%%%%%;

if (verbose); disp(sprintf(' %% ')); end;
if (verbose); disp(sprintf(' %% Now we visualize the dosage-data for the cases (top) and ctrls (bottom). ;')); end;
rij_case_prm_ = [ ... 
; rij_case_bcA_ ...
; rij_case_bcB_ ;
; setdiff(rij_case_,union(rij_case_bcA_,rij_case_bcB_)) ...
];
rij_ctrl_prm_ = [ ... 
; rij_ctrl_bcA_ ...
; rij_ctrl_bcB_ ;
; setdiff(rij_ctrl_,union(rij_ctrl_bcA_,rij_ctrl_bcB_)) ...
];
rij_prm_ = [ rij_case_prm_ ; rij_ctrl_prm_ ];
cij_prm_ = [ ...
; cij_bcA_ ...
; cij_bcB_ ...
; setdiff(cij_,union(cij_bcA_,cij_bcB_)) ...
];
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,512,1024]); fig80s;
subplot(2,1,1); imagesc(dsg_cup_ps__(rij_case_prm_,cij_prm_)); axisnotick; ylabel('case');
subplot(2,1,2); imagesc(dsg_cup_ps__(rij_ctrl_prm_,cij_prm_)); axisnotick; ylabel('ctrl');

if (verbose); disp(sprintf(' %% ')); end;
if (verbose); disp(sprintf(' %% Now we convert the dosage-data to b16 data. ;')); end;
if (verbose); disp(sprintf(' %% We put everything into a single study. ;')); end;
%%%%%%%%;
dir_trunk = sprintf('%s/dir_test_xxxcluster_fromdisk_ver16',pwd);
if (~exist(dir_trunk,'dir')); disp(sprintf(' %% mkdir %s',dir_trunk)); mkdir(dir_trunk); end;
n_study = 1; nstudy = 0;
study_name = sprintf('study%.2d',nstudy);
study_trunk = sprintf('dir_%s',study_name);
tmp_dir = sprintf('%s/%s',dir_trunk,study_trunk);
if (~exist(tmp_dir,'dir')); disp(sprintf(' %% mkdir %s',tmp_dir)); mkdir(tmp_dir); end;
n_snp = n_snp_cup;
index_nsnp_from_nstudy_ = 0:n_snp-1;
n_patient = n_patient_cup;
index_npatient_from_nstudy_ = 0:n_patient-1;
%%%%%%%%;
fn_fam = sprintf('%s/%s/%s.fam',dir_trunk,study_trunk,study_name);
fid = fopen(fn_fam,'w');
for npatient=0:n_patient-1;
np = index_npatient_from_nstudy_(1+npatient);
fprintf( ...
 fid ...
,'%s\t%s\t%s\t%s\t%d\t%d\n' ...
,fam_cup_fid_{1+np} ...
,fam_cup_iid_{1+np} ...
,fam_cup_yid_{1+np} ...
,fam_cup_xid_{1+np} ...
,fam_cup_sex_(1+np) ...
,fam_cup_dvx_(1+np) ...
);
end;%for npatient=0:n_patient-1;
fclose(fid);
%%%%%%%%;
fn_bim = sprintf('%s/%s/%s.bim',dir_trunk,study_trunk,study_name);
fid = fopen(fn_bim,'w');
for nsnp=0:n_snp-1;
ns = index_nsnp_from_nstudy_(1+nsnp);
fprintf( ...
 fid ...
,'%d\t%s\t%d\t%d\t%s\t%s\n' ...
,bim_cup_khr_(1+ns) ...
,bim_cup_vid_{1+ns} ...
,bim_cup_pos_(1+ns) ...
,bim_cup_bpc_(1+ns) ...
,bim_cup_al1_{1+ns} ...
,bim_cup_al2_{1+ns} ...
);
end;%for nsnp=0:n_snp-1;
fclose(fid);
%%%%%%%%;
fn_bed = sprintf('%s/%s/%s.bed',dir_trunk,study_trunk,study_name);
fid = fopen(fn_bed,'w');
tmp_dsg_ps__ = dsg_cup_ps__(1+index_npatient_from_nstudy_,1+index_nsnp_from_nstudy_);
[tmp_bed_ps__,n_patient_bed] = dsg_to_bed_0(n_patient,n_snp,tmp_dsg_ps__);
key1 = uint8(108); key2 = uint8(27); key3 = uint8(1);
fwrite(fid,key1,'uint8'); fwrite(fid,key2,'uint8'); fwrite(fid,key3,'uint8');
fwrite(fid,tmp_bed_ps__,'uint8');
fclose(fid);
%%%%%%%%;
n_mds_p = n_patient_cup;
n_mds_v = 2;
mds_pv__ = randn(n_mds_p,n_mds_v);
mds_fidandiid_p_ = cell(n_mds_p,1);
for npatient_cup=0:n_patient_cup-1;
mds_fidandiid_p_{1+npatient_cup} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for npatient_cup=0:n_patient_cup-1;
%%%%%%%%;
n_famex = 0;
index_npatient_cup_from_nfamex_ = randperm(n_patient_cup,n_famex)-1;
famex_fidandiid_p_ = cell(n_famex,1);
for nfamex=0:n_famex-1;
npatient_cup = index_npatient_cup_from_nfamex_(1+nfamex);
famex_fidandiid_p_{1+nfamex} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for nfamex=0:n_famex-1;
%%%%%%%%;
str_output_prefix = 'test';
parameter = struct('type','parameter');
parameter.verbose = 1;
parameter.str_output_prefix = str_output_prefix;
parameter.dir_trunk = dir_trunk;
parameter.ent_cutoff = 0.001;
parameter.maf_cutoff = 0.010;
[ ...
 parameter ...
] = ...
bed_to_b16_flip_ver7( ...
 parameter ...
,n_study ...
,{study_trunk} ...
,{study_name} ...
,n_mds_p ...
,n_mds_v ...
,mds_pv__ ...
,mds_fidandiid_p_ ...
,n_famex ...
,famex_fidandiid_p_ ...
);
%%%%%%%%;
str_maf = sprintf('maf%.2d',floor(100*parameter.maf_cutoff));
str_output_prefix_plus_maf = sprintf('%s_%s',str_output_prefix,str_maf);
str_output_prefix_local = sprintf('%s_',str_output_prefix_plus_maf);
dir_out = sprintf('%s/dir_%s',parameter.dir_trunk,str_output_prefix_plus_maf);
if (verbose>0); disp(sprintf(' %% str_output_prefix_local %s ;\n %% dir_out %s ;\n',str_output_prefix_local,dir_out)); end;

fname_b16 = sprintf('%s/%s_A_full_n.b16',dir_out,str_output_prefix_plus_maf);
A_full_n_ = binary_uncompress(fname_b16);
%%%%%%%%;
fname_famext = sprintf('%s/%s_fam.ext',dir_out,str_output_prefix_plus_maf);
[ ...
 n_patient_famext ...
,famext_fid_ ...
,famext_iid_ ...
,famext_yid_ ...
,famext_xid_ ...
,famext_sex_ ...
,famext_dvx_ ...
,famext_dir_ ...
,famext_fidandiid_ ...
,famext_ ...
  ] = ...
load_famext_ver1( ...
 fname_famext ...
);
[~,~,ij_fam_ext_from_cap_] = intersect(fam_cup_fid_,famext_fid_,'stable');
disp(sprintf(' %% famext_fid_ vs fam_cup_fid_: %0.16f',fnorm(ij_fam_ext_from_cap_-transpose(1:n_patient))));
%%%%%%%%;

fname_bimext = sprintf('%s/%s_bim.ext',dir_out,str_output_prefix_plus_maf);
[ ...
 n_snp_bimext ...
,bimext_khr_ ...
,bimext_vid_ ...
,bimext_gdi_ ...
,bimext_pdi_ ...
,bimext_al1_ ...
,bimext_al2_ ...
,bimext_alt_ ...
,bimext_ent_ ...
,bimext_frq_ ...
,bimext_mss_ ...
,bimext_maf_ ...
,bimext_name_ ...
,bimext_ ...
] = ...
load_bimext_ver1( ...
fname_bimext ...
);
[u_bimext_vid_,index_bim_ext_from_u_,index_bim_u_from_ext_] = unique(bimext_vid_,'stable');
index_bim_ext_from_u_ = index_bim_ext_from_u_ - 1;
index_bim_u_from_ext_ = index_bim_u_from_ext_ - 1;
[bim_vid_cap_,index_bim_cup_from_cap_,index_bim_u_from_cap_] = intersect(bim_cup_vid_,u_bimext_vid_,'stable');
index_bim_cup_from_cap_ = index_bim_cup_from_cap_ - 1;
index_bim_u_from_cap_ = index_bim_u_from_cap_ - 1;
[~,~,index_bim_cap_from_u_] = intersect(u_bimext_vid_,bim_vid_cap_,'stable');
index_bim_cap_from_u_ = index_bim_cap_from_u_ - 1;
index_bim_ext_from_cap_ = index_bim_ext_from_u_(1+index_bim_u_from_cap_);
index_bim_cup_from_ext_ = index_bim_cup_from_cap_(1+index_bim_cap_from_u_(1+index_bim_u_from_ext_));
ij_bim_cup_from_ext__ = sparse(1+index_bim_cup_from_ext_,1:n_snp_bimext,1,n_snp_cup,n_snp_bimext);
[~,cij_prm_inv_] = sort(cij_prm_,'ascend');
cij_prm_inv_from_ext_ = transpose(cij_prm_inv_) * ij_bim_cup_from_ext__;
[~,cij_ext_] = sort(cij_prm_inv_from_ext_,'ascend');
%%%%;
mr_bcx_ = zeros(n_patient,1);
mr_bcx_(rij_case_bcA_) = 2;
mr_bcx_(rij_case_bcB_) = 1;
mr_bcx_(rij_ctrl_bcB_) = 1;
%%%%;
mc_bcA_ = zeros(n_snp_bimext,1);
mc_bcB_ = zeros(n_snp_bimext,1);
for nsnp_bimext=0:n_snp_bimext-1;
mc_bcA_(1+nsnp_bimext) = sum(strcmp(bimext_vid_{1+nsnp_bimext},bim_cup_vid_(cij_bcA_)));
mc_bcB_(1+nsnp_bimext) = sum(strcmp(bimext_vid_{1+nsnp_bimext},bim_cup_vid_(cij_bcB_)));
end;%for nsnp_bimext=0:n_snp_bimext-1;
cij_ex2_ = [ find(mc_bcA_) ; find(mc_bcB_) ; setdiff(transpose(1:n_snp_bimext),find(mc_bcA_ + mc_bcB_)) ];
%%%%;
mc_bcx_ = zeros(n_snp_bimext,1);
mc_bcx_(find(mc_bcA_))=2;
mc_bcx_(find(mc_bcB_))=1;
n_snp_bcx = numel(find(mc_bcx_));
n_snp_ext_bcA = numel(find(mc_bcx_==2));
n_snp_ext_bcB = numel(find(mc_bcx_==1));
%%%%;
if (verbose); disp(sprintf(' %% ')); end;
if (verbose); disp(sprintf(' %% Now we visualize the b16-data for the cases (top) and ctrls (bottom). ;')); end;
if (verbose); disp(sprintf(' %% Column-1 orders the b16 data by snp-index, Column-2 by allele-frequency. ')); end;
if (verbose); disp(sprintf(' %% the column-indices corresponding to bcA are shown first, followed by bcB. ')); end;
if (verbose); disp(sprintf(' %% Column-3 shows the signal-strength of each bicluster in the cases alone. ')); end;
if (verbose); disp(sprintf(' %% Column-4 shows the case-vs-ctrl signal-strength of each bicluster. ')); end;
if (verbose); disp(sprintf(' %% Note that bcA is case-specific, while bcB extends across the cases and ctrls. ')); end;
%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,1024]); fig80s;
p_row = 2; p_col = 4;
subplot(p_row,p_col,1 + 0*p_col); imagesc(A_full_n_(rij_case_prm_,cij_ext_(1:n_snp_bcx))); axisnotick; ylabel('case'); title('ext');
subplot(p_row,p_col,1 + 1*p_col); imagesc(A_full_n_(rij_ctrl_prm_,cij_ext_(1:n_snp_bcx))); axisnotick; ylabel('ctrl'); title('ext');
subplot(p_row,p_col,2 + 0*p_col); imagesc(A_full_n_(rij_case_prm_,cij_ex2_(1:n_snp_bcx))); axisnotick; ylabel('case'); title('ex2');
subplot(p_row,p_col,2 + 1*p_col); imagesc(A_full_n_(rij_ctrl_prm_,cij_ex2_(1:n_snp_bcx))); axisnotick; ylabel('ctrl'); title('ex2');
subplot(p_row,p_col,3 + 0*p_col); hold on;
plot( 1:n_snp_bimext, mean(A_full_n_(rij_case_bcA_,cij_ext_),1) - mean(A_full_n_(rij_case_bcA_cmp_,cij_ext_),1) , 'k.' );
plot(n_snp_ext_bcA*[1,1],[-1,+1],'k-','LineWidth',3 );
xlim([1,n_snp_bcx]); ylim([-0.5,+0.5]); ylabel('signal strength'); title('bcA case');
subplot(p_row,p_col,3 + 1*p_col); hold on;
plot( 1:n_snp_bimext, mean(A_full_n_(rij_case_bcB_,cij_ext_),1) - mean(A_full_n_(rij_case_bcB_cmp_,cij_ext_),1) , 'k.' );
plot(n_snp_ext_bcA*[1,1],[-1,+1],'k-','LineWidth',3 );
xlim([1,n_snp_bcx]); ylim([-0.5,+0.5]); ylabel('signal strength'); title('bcB case');
subplot(p_row,p_col,4 + 0*p_col); hold on;
plot( 1:n_snp_bimext, mean(A_full_n_(rij_case_bcA_,cij_ext_),1) - mean(A_full_n_(rij_ctrl_,cij_ext_),1) , 'k.' );
plot(n_snp_ext_bcA*[1,1],[-1,+1],'k-','LineWidth',3 );
xlim([1,n_snp_bcx]); ylim([-0.5,+0.5]); ylabel('signal strength'); title('bcA case-vs-ctrl');
subplot(p_row,p_col,4 + 1*p_col); hold on;
plot( 1:n_snp_bimext, mean(A_full_n_(rij_case_bcB_,cij_ext_),1) - mean(A_full_n_(rij_ctrl_bcB_,cij_ext_),1) , 'k.' );
plot(n_snp_ext_bcA*[1,1],[-1,+1],'k-','LineWidth',3 );
xlim([1,n_snp_bcx]); ylim([-0.5,+0.5]); ylabel('signal strength'); title('bcB case-vs-ctrl');
%%%%;

n_mds_0in = 0; n_mds_repl = 0; ij_mds_use_ = [];
gamma = 0.05;
dir_code = sprintf('%s/..',pwd);
str_prefix = str_output_prefix_plus_maf;
parameter.str_lak_vs_dex = str_lak_vs_dex;
parameter.str_prefix = str_prefix;
parameter.dir_code = dir_code;
parameter.gamma = gamma;
parameter.n_mds = n_mds_0in;
parameter.n_mds_repl = n_mds_repl;
parameter.ij_mds_use_ = ij_mds_use_;
parameter.flag_rerun = 0;
[ ...
 parameter ...
] = ...
xxxcluster_fromdisk_uADZSZDA_ver16( ...
 parameter ...
);
parameter.str_trace_s0000 = sprintf('%s/out_trace_s0000.txt',parameter.dir_out_trace);
[ ...
 niter_ ...
,r_rem_ ...
,c_rem_ ...
,QR_ ...
,QC_ ...
,I_rem_ ...
,trace_ ...
] = ...
load_trace_ver0( ...
parameter.str_trace_s0000 ...
);
parameter.str_out_xdrop_a_s0000 = sprintf('%s/out_xdrop_a.txt',parameter.dir_out_s0000);
[ ...
 rindex_ ...
,cindex_ ...
,out_xdrop_ ...
] = ...
load_out_xdrop_ver0( ...
parameter.str_out_xdrop_a_s0000 ...
);
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,1024]); fig80s;
imagesc( [ A_full_n_(1+flip(rindex_),1+flip(cindex_)) ; A_full_n_(rij_ctrl_prm_,1+flip(cindex_)) ] );
