% file to generate standard CTD data
% can be thought of as different sub-functions
% B zheng
% Dec. 21, 2020
% - % - % - % - % - % - % - % - % - % - % - % - % - % - %- %
%% convert .rsk file into matlab file
clear all
WWmeta.rbrpath = '/Users/Devon/Documents/GradSchool/ASTRAL/RBR_Raw/DBASIS/Upper/';
WWmeta.rbrfile=dir(fullfile(WWmeta.rbrpath,'*.rsk'));
WWmeta.name_rbr = 'ASTRAL_DBASIS_Upper_2';
WWmeta.matpath = '/Users/Devon/Documents/GradSchool/ASTRAL/WW_Processed/DBASIS/Upper/RBR/mat/';  % path for converted mat file
WWmeta.propath = '/Users/Devon/Documents/GradSchool/ASTRAL/WW_Processed/DBASIS/Upper/RBR/profile/';  % path for saving upcast profiles
WWmeta.gridpath = '/Users/Devon/Documents/GradSchool/ASTRAL/WW_Processed/DBASIS/Upper/RBR/grid/';  % path for the gridded product
WWmeta.rsktoolspath = '/Users/Devon/Documents/MATLAB/rbr-rsktoolsv353/';
addpath(genpath(WWmeta.rsktoolspath))
WWmeta.lat = 12;
WWmeta.lon = 85;

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







