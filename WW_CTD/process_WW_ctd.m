% file to generate standard CTD data
% can be thought of as different sub-functions
% B zheng
% Dec. 21, 2020
% - % - % - % - % - % - % - % - % - % - % - % - % - % - %- %
%% convert .rsk file into matlab file
clear all
WWmeta.rbrpath = '/Users/Devon/Documents/GradSchool/DTS/Wirewalker/RBR/raw/';
WWmeta.rbrfile=dir(fullfile(WWmeta.rbrpath,'*.rsk'));
WWmeta.name_rbr = 'TLC23_dye_ThorpeScale';
WWmeta.matpath = '/Users/Devon/Documents/GradSchool/DTS/Wirewalker/RBR/mat/';  % path for converted mat file
WWmeta.propath = '/Users/Devon/Documents/GradSchool/DTS/Wirewalker/RBR/profile/';  % path for saving upcast profiles
WWmeta.gridpath = '/Users/Devon/Documents/GradSchool/DTS/Wirewalker/RBR/grid/';  % path for the gridded product
WWmeta.lat = 33;
WWmeta.lon = -117.5;

WWmeta.thorpescales.yn = 1; % Do thorpe scale analysis (1=yes, 0=no);

WWmeta.salspiking.yn = 1;  % Do salinity de-spiking? (1=yes, 0=no);
WWmeta.salspiking.win_len = 80; % window size for sal spiking calculation (I set to ~10s)
WWmeta.salspiking.Fs = 8; % CTD sampling frequency
WWmeta.salspiking.iPlot = 1; % Plot sal de-spiking stuff
WWmeta.salspiking.poly_order = 20; % Polynomial order for sal de-spiking fits
WWmeta.salspiking.dTorT = 1;  % 0 = use derivite of T and C, 1 = use raw data. Sould be equivelent?

WW_rskread(WWmeta);

%% read matlab files
WWmatread(WWmeta);

%% separate into profiles and then grid
WWprofile(WWmeta,15);  % 30 is subject to change

%% grid the CTD product
WWgrid(WWmeta,0.25);

%% sample plot - t/s/sig/chla
para.tscale = [32, 24];
para.sscale = [32, 35];
para.dscale = [19 24];
para.cscale = [1.7 2];

ctdplot(WWmeta,para)







