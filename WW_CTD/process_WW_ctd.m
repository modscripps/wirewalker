% file to generate standard CTD data
% can be thought of as different sub-functions
% B zheng
% Dec. 21, 2020
% - % - % - % - % - % - % - % - % - % - % - % - % - % - %- %
%% convert .rsk file into matlab file
clear all
WWmeta.rbrpath = '/Volumes/NorseTPADS/NORSE_ASTRAL_raw_data/NORSE23/DBASIS_upper_D2/';
WWmeta.rbrfile=dir(fullfile(WWmeta.rbrpath,'*.rsk'));
WWmeta.name_rbr = 'NORSE_23_WWH_D2';
WWmeta.matpath = '/Volumes/NorseTPADS/WW_NORSE_23/RBR/WW_H_D2/mat/';  % path for converted mat file
WWmeta.propath = '/Volumes/NorseTPADS/WW_NORSE_23/RBR/WW_H_D2/profile/';  % path for saving upcast profiles
WWmeta.gridpath = '/Volumes/NorseTPADS/WW_NORSE_23/RBR/WW_H_D2/grid/';  % path for the gridded product
WWmeta.rsktoolspath = '/Users/Devon/Documents/MATLAB/rbr-rsktoolsv353/';
addpath(genpath(WWmeta.rsktoolspath))
WWmeta.lat = 71;
WWmeta.lon = -6.5;

WWmeta.thorpescales.yn = 0; % Do thorpe scale analysis (1=yes, 0=no);

WWmeta.salspiking.yn = 1;  % Do salinity de-spiking? (1=yes, 0=no);
WWmeta.salspiking.win_len = 80; % window size for sal spiking calculation (I set to ~10s)
WWmeta.salspiking.Fs = 8; % CTD sampling frequency
WWmeta.salspiking.iPlot = 1; % Plot sal de-spiking stuff
WWmeta.salspiking.poly_order = 20; % Polynomial order for sal de-spiking fits
WWmeta.salspiking.dTorT = 1;  % 0 = use derivite of T and C, 1 = use raw data. Sould be equivelent?

WWmeta.LineNoiseSquasher = 0; % To fix a specific line noise issue in the ASTRAL DBASIS lower temperature record.

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







