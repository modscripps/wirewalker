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
            if ~isfield(out,'T')
                out.T=RSKread.data.values(:,c);
            elseif isfield(out,'T')
                for n = 1:10
                    if isfield(out,['T' num2str(n)])
                        continue
                    else
                        eval(['out.T' num2str(n) '=RSKread.data.values(:,c);']);
                        break
                    end
                end
            end
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
            if ~isfield(out,'T')
                out.T=RSKread.data.values(:,c);
            elseif isfield(out,'T')
                for n = 1:10
                    if isfield(out,['T' num2str(n)])
                        continue
                    else
                        eval(['out.T' num2str(n) '=RSKread.data.values(:,c);']);
                        break
                    end
                end
            end
        case 'Temperature2'
            if ~isfield(out,'T')
                out.T=RSKread.data.values(:,c);
            elseif isfield(out,'T')
                for n = 1:10
                    if isfield(out,['T' num2str(n)])
                        continue
                    else
                        eval(['out.T' num2str(n) '=RSKread.data.values(:,c);']);
                        break
                    end
                end
            end
        case 'Dissolved O concentration'
            out.O2=RSKread.data.values(:,c);
        case 'Chlorophyll-a'
            out.chla=RSKread.data.values(:,c); 
        case 'FDOM'
            out.fdom=RSKread.data.values(:,c);
        case 'Sea Pressure'
            out.sea_pressure=RSKread.data.values(:,c);
        case 'Depth'
            out.depth=RSKread.data.values(:,c);
        case 'Speed of sound'
            out.sound_speed=RSKread.data.values(:,c);
        case 'Specific conductivity'
            out.specific_conductivity=RSKread.data.values(:,c);
        case 'Dissolved O saturation'
            out.O2_sat=RSKread.data.values(:,c);
        case 'Irradiance'
            if ~isfield(out,'irr1')
                out.irr1=RSKread.data.values(:,c);
            elseif isfield(out,'irr1')
                for n = 2:10
                    if isfield(out,['irr' num2str(n)])
                        continue
                    else
                        eval(['out.irr' num2str(n) '=RSKread.data.values(:,c);']);
                        break
                    end
                end
            end
        case 'Dissolved O21'
            out.O2=RSKread.data.values(:,c);
        case 'Dissolved O22'
            out.O2_sat=RSKread.data.values(:,c);
        case 'CT Cell Temperature'
            out.T_CTcell=RSKread.data.values(:,c);
        case 'Pressure Gauge Temperature'
            out.T_PressureGuage=RSKread.data.values(:,c);
        case 'Salinity'
            out.S_RBR=RSKread.data.values(:,c);
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
out.rho = gsw_pot_rho_t_exact(SA,out.T,out.P,0);
out.sig0 = out.rho-1000;

CTDall = out;
%% save ctd all
save([WWmeta.matpath,WWmeta.name_rbr,'_CTDall.mat'],'CTDall');
end