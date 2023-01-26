function ...
[ ...
 parameter ...
, AZnV_0_lim_ ...
, AZnV_0_tick_ ...
, AZnV_1_lim_ ...
, AZnV_1_tick_ ...
] = ...
test_scatter_and_heatmap_0( ...
 parameter ...
, mr_dvx_ ...
, AZnV_pr__ ...
, AZnV_0_lim_ ...
, AZnV_0_tick_ ...
, AZnV_1_lim_ ...
, AZnV_1_tick_ ...
);

str_thisfunction = 'test_scatter_and_heatmap_0';

na=0;
if nargin<1+na; parameter=[]; end; na=na+1;
if nargin<1+na; mr_dvx_=[]; end; na=na+1;
if nargin<1+na; AZnV_pr__=[]; end; na=na+1;
if nargin<1+na; AZnV_0_lim_=[]; end; na=na+1;
if nargin<1+na; AZnV_0_tick_=[]; end; na=na+1;
if nargin<1+na; AZnV_1_lim_=[]; end; na=na+1;
if nargin<1+na; AZnV_1_tick_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp = 0; end;
flag_disp = parameter.flag_disp;
if ~isfield(parameter,'markersize_sml'); parameter.markersize_sml = 12; end;
markersize_sml = parameter.markersize_sml;
if ~isfield(parameter,'markersize_big'); parameter.markersize_big = 16; end;
markersize_big = parameter.markersize_big;
if ~isfield(parameter,'fontsize_use'); parameter.fontsize_use = 12; end;
fontsize_use = parameter.fontsize_use;
if ~isfield(parameter,'n_pcol'); parameter.n_pcol = 8; end;
n_pcol = parameter.n_pcol;
if ~isfield(parameter,'n_pgap'); parameter.n_pgap = 2; end;
n_pgap = parameter.n_pgap;
if ~isfield(parameter,'n_h_0'); parameter.n_h_0 = 32+1; end;
n_h_0 = parameter.n_h_0;
if ~isfield(parameter,'n_h_1'); parameter.n_h_1 = 32+1; end;
n_h_1 = parameter.n_h_1;
if ~isfield(parameter,'n_tick_0'); parameter.n_tick_0 = 2; end;
n_tick_0 = parameter.n_tick_0;
if ~isfield(parameter,'n_tick_1'); parameter.n_tick_1 = 2; end;
n_tick_1 = parameter.n_tick_1;
if ~isfield(parameter,'h_DvX_scale'); parameter.h_DvX_scale = 0.025; end;
h_DvX_scale = parameter.h_DvX_scale;
if ~isfield(parameter,'box_expand'); parameter.box_expand = 1.25; end;
box_expand = parameter.box_expand;
if ~isfield(parameter,'c_use__'); parameter.c_use__ = colormap_pm; end;
c_use__ = parameter.c_use__;
if ~isfield(parameter,'legend_use_'); parameter.legend_use_ = {'ctrl','case','bicl'}; end;
legend_use_ = parameter.legend_use_;
if ~isfield(parameter,'legend_location'); parameter.legend_location = 'NorthWest'; end;
legend_location = parameter.legend_location;
%if ~isfield(parameter,'zzz'); parameter.zzz = ; end;
%zzz = parameter.zzz;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
tmp_index_ = efind( (mr_dvx_> 0) ); %<-- ctrl and case ;
if isempty(AZnV_0_lim_);
AZnV_0_min = min(AZnV_pr__(1+tmp_index_,1+0)); AZnV_0_max = max(AZnV_pr__(1+tmp_index_,1+0));
AZnV_0_lim_ = [AZnV_0_min,AZnV_0_max];
AZnV_0_lim_ = mean(AZnV_0_lim_) + box_expand*0.5*diff(AZnV_0_lim_)*[-1,+1];
end;%if isempty(AZnV_0_lim_);
if isempty(AZnV_0_tick_);
AZnV_0_tick_ = linspace(min(AZnV_0_lim_),max(AZnV_0_lim_),n_h_0);
end;%if isempty(AZnV_0_tick_);
if isempty(AZnV_1_lim_);
AZnV_1_min = min(AZnV_pr__(1+tmp_index_,1+1)); AZnV_1_max = max(AZnV_pr__(1+tmp_index_,1+1));
AZnV_1_lim_ = [AZnV_1_min,AZnV_1_max];
AZnV_1_lim_ = mean(AZnV_1_lim_) + box_expand*0.5*diff(AZnV_1_lim_)*[-1,+1];
end;%if isempty(AZnV_1_lim_);
if isempty(AZnV_1_tick_);
AZnV_1_tick_ = linspace(min(AZnV_1_lim_),max(AZnV_1_lim_),n_h_1);
end;%if isempty(AZnV_1_tick_);
%%%%%%%%;

%%%%%%%%;
tmp_index_ = efind( (mr_dvx_==1) ); %<-- ctrl ;
tmp_h_X_01__ = hist2d_0(AZnV_pr__(1+tmp_index_,1+0),AZnV_pr__(1+tmp_index_,1+1),n_h_0,n_h_1,AZnV_0_lim_,AZnV_1_lim_);
tmp_h_X_01__ = tmp_h_X_01__./max(1,sum(tmp_h_X_01__,'all'));
tmp_index_ = efind( (mr_dvx_==2) | (mr_dvx_==3) ); %<-- case ;
tmp_h_D_01__ = hist2d_0(AZnV_pr__(1+tmp_index_,1+0),AZnV_pr__(1+tmp_index_,1+1),n_h_0,n_h_1,AZnV_0_lim_,AZnV_1_lim_);
tmp_h_D_01__ = tmp_h_D_01__./max(1,sum(tmp_h_D_01__,'all'));
tmp_h_DvX_01__ = tmp_h_D_01__ - tmp_h_X_01__;
%%%%%%%%;

%figure(1); clf; set(gcf,'Position',1+[0,0,1024*2,768]);

%%%%%%%%;
subplot_{1} = subplot(1,n_pcol+n_pgap+n_pcol+1,0*n_pcol+0*n_pgap+[1:n_pcol]); cla;
hold on;
tmp_index_ = efind( (mr_dvx_==1) );
scatter(AZnV_pr__(1+tmp_index_,1+0),AZnV_pr__(1+tmp_index_,1+1),markersize_sml,0*ones(numel(tmp_index_),1),'filled','MarkerEdgeColor',0.75*[1,1,1]);
tmp_index_ = efind( (mr_dvx_==2) );
scatter(AZnV_pr__(1+tmp_index_,1+0),AZnV_pr__(1+tmp_index_,1+1),markersize_sml,1*ones(numel(tmp_index_),1),'filled','MarkerEdgeColor',0.75*[1,1,1]);
tmp_index_ = efind( (mr_dvx_==3) );
if numel(tmp_index_)> 0;
scatter(AZnV_pr__(1+tmp_index_,1+0),AZnV_pr__(1+tmp_index_,1+1),markersize_big,1*ones(numel(tmp_index_),1),'filled','MarkerEdgeColor',0.05*[1,1,1]);
end;%if numel(tmp_index_)> 0;
legend(legend_use_(1:max(mr_dvx_)),'Location',legend_location);
xlim(AZnV_0_lim_); ylim(AZnV_1_lim_); grid on; set(gca,'TickLength',[0,0]);
set(gca,'XTick',AZnV_0_tick_(1:n_tick_0:end),'XTickLabel',num2str(transpose(AZnV_0_tick_(1:n_tick_0:end)),'%+.1f'));
set(gca,'YTick',AZnV_1_tick_(1:n_tick_1:end),'YTickLabel',num2str(transpose(AZnV_1_tick_(1:n_tick_1:end)),'%+.1f'));
xtickangle(90); set(gca,'FontSize',fontsize_use);
title('scatter','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
%%%%%%%%;

%%%%%%%%;
subplot_{2} = subplot(1,n_pcol+n_pgap+n_pcol+1,1*n_pcol+1*n_pgap+[1:n_pcol]); cla;
imagesc(tmp_h_DvX_01__,h_DvX_scale*[-1,+1]);
set(gca,'ydir','normal'); set(gca,'TickLength',[0,0]);
set(gca,'XTick',1:n_tick_0:n_h_0,'XTickLabel',num2str(transpose(AZnV_0_tick_(1:n_tick_0:end)),'%+.1f'));
set(gca,'YTick',1:n_tick_0:n_h_1,'YTickLabel',num2str(transpose(AZnV_1_tick_(1:n_tick_1:end)),'%+.1f'));
xtickangle(90); set(gca,'FontSize',fontsize_use);
title('heatmap','Interpreter','none'); xlabel('PC1'); ylabel('PC2');
%%%%%%%%;

%%%%%%%%;
subplot_{3} = subplot(1,n_pcol+n_pgap+n_pcol+1,1+2*n_pcol+1*n_pgap); cla;
imagesc(transpose(1:64),[1,64]);
set(gca,'ydir','normal'); set(gca,'YAxisLocation','right'); set(gca,'FontSize',fontsize_use);
set(gca,'XTick',[]); set(gca,'YTick',[1,64],'YTickLabel',num2str(h_DvX_scale*[-1;+1],'%+.3f')); set(gca,'TickLength',[0,0]);
%%%%%%%%;

%%%%;
colormap(subplot_{1},c_use__);
colormap(subplot_{2},c_use__);
colormap(subplot_{3},c_use__);
%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
