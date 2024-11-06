function ctdplot(WWmeta,para)
% plot a raw ctd figure
% bz, jan 15, 2021

cbt   = cbrewer('div','RdYlBu',64,'cubic'); % for temperature
cbt   = flipud(cbt); cbt(cbt>1) = 1;
cbs   = cbrewer('div','RdYlGn',64,'cubic'); % for salinity
cbs   = flipud(cbs); cbs(cbs>1) = 1;
cbd   = cbrewer('div','PRGn',64,'cubic'); % for salinity
% cbd   = flipud(cbd); 
cbs(cbd>1) = 1;


load([WWmeta.gridpath,WWmeta.name_rbr,'_CTDgrid.mat']);

title_string = WWmeta.name_rbr;
title_string(strfind(title_string,'_'))=' ';

figure
ax(1) = subplot('position',[0.1,0.75,0.8,0.195]);
time = RBRgrid.std_profiles.time;
z = -RBRgrid.std_profiles.z;
h = pcolor(time,z,RBRgrid.std_profiles.T);
set(h, 'EdgeColor', 'none');
colormap(ax(1),cbt);
caxis([para.tscale(1) para.tscale(2)])
hclb = colorbar;
set(hclb,'location','eastoutside'); %colorbar location
set(hclb,'Position',[0.92 0.75 0.015 0.195]);
danum = min(time);
% set(gca,'xtick',danum:(max(time)-danum)/10:max(time),'tickdir','in');
set(gca,'xtick',danum:1:max(time),'tickdir','in');
% datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(time) max(time)])
title(['\fontsize{15}',title_string,' Temperature'])
set(gca,{'xticklabel'},{[]});
% hxlb = xlabel('Time (DD HH:MM)');
% set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)


ax(2) = subplot('position',[0.1,0.525,0.8,0.195]);
h = pcolor(time,z,RBRgrid.std_profiles.S);
set(h, 'EdgeColor', 'none');
colormap(ax(2),cbs);
caxis([para.sscale(1) para.sscale(2)])
hclb = colorbar;
set(hclb,'location','eastoutside'); %colorbar location
set(hclb,'Position',[0.92 0.525 0.015 0.195]);
danum = min(time);
% set(gca,'xtick',danum:(max(time)-danum)/10:max(time),'tickdir','in');
set(gca,'xtick',danum:1:max(time),'tickdir','in');
% datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(time) max(time)])
title(['\fontsize{15}',title_string,' Salinity'])
set(gca,{'xticklabel'},{[]});


ax(3) = subplot('position',[0.1,0.3,0.8,0.195]);
h = pcolor(time,z,RBRgrid.std_profiles.sig0);
set(h, 'EdgeColor', 'none');
colormap(ax(3),cbd);
caxis([para.dscale(1) para.dscale(2)])
hclb = colorbar;
set(hclb,'location','eastoutside'); %colorbar location
set(hclb,'Position',[0.92 0.3 0.015 0.195]);
danum = min(time);
% set(gca,'xtick',danum:(max(time)-danum)/10:max(time),'tickdir','in');
set(gca,'xtick',danum:1:max(time),'tickdir','in');
% datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(time) max(time)])
title(['\fontsize{15}',title_string,' Density'])
set(gca,{'xticklabel'},{[]});


ax(4) = subplot('position',[0.1,0.075,0.8,0.195]);
h = pcolorjw(time,z,log10(RBRgrid.std_profiles.chla));
set(h, 'EdgeColor', 'none');
colormap(ax(4),'jet');
caxis([para.cscale(1) para.cscale(2)])
hclb = colorbar;
set(hclb,'location','eastoutside'); %colorbar location
set(hclb,'Position',[0.92 0.075 0.015 0.195]);
danum = min(time);
% set(gca,'xtick',danum:(max(time)-danum)/10:max(time),'tickdir','in');
set(gca,'xtick',danum:1:max(time),'tickdir','in');
datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(time) max(time)])
title(['\fontsize{15}',title_string,' Chla (log10(counts))'])
hxlb = xlabel('Time (DD HH:MM)');
set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)

linkaxes(ax,'xy')
end