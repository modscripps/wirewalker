function WWmatread(WWmeta)
% convert data from rskread to physical variables
% bz, june 15, 2021
%%

load([WWmeta.rbrpath,WWmeta.name_rbr])

%%
i=0;
for c=1:size(RSKread.data.values,2) % 
    switch RSKread.channels(c).longName
        case 'Pressure'
            out.P=RSKread.data.values(:,c)-10.13;  % pressure offset
        case 'Temperature'
            out.T=RSKread.data.values(:,c);
        case 'Conductivity'
            out.C=RSKread.data.values(:,c);
        case 'Backscatter'
            out.bs=RSKread.data.values(:,c);
        case 'Chlorophyll'
            out.chla=RSKread.data.values(:,c);
        case 'CDOM'
            out.cdom=RSKread.data.values(:,c);
        case 'PAR'
            out.par=RSKread.data.values(:,c);
        case 'Irradiance1'
            out.irr1=RSKread.data.values(:,c);
        case 'Irradiance2'
            out.irr2=RSKread.data.values(:,c);
        case 'Irradiance3'
            out.irr3=RSKread.data.values(:,c);
        case 'Temperature1'
            out.T=RSKread.data.values(:,c); 
        otherwise
            i=i+1;
            eval(sprintf('out.v%i=RSKread.data.values(:,c)',i));
            eval(sprintf('out.info.v%i=RSKread.channels(c).longName;',i));
    end  
end

if WWmeta.LineNoiseSquasher == 1
    out.T = SquashLineNoise(out.T,1000);
end

out.time=(RSKread.data.tstamp);
out.P(out.P<0)=0;
out.S   = gsw_SP_from_C(out.C,out.T,out.P);
SA      = gsw_SA_from_SP(out.S,out.P,WWmeta.lon,WWmeta.lat);  % lat/lon are subject to change
CT      = gsw_CT_from_t(SA,out.T,out.P);
out.rho = gsw_rho(SA,CT,0);
out.sig0 = out.rho-1000;

CTDall = out;
%% save ctd all
save([WWmeta.matpath,WWmeta.name_rbr,'_CTDall.mat'],'CTDall');
end