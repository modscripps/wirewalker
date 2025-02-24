function WWprofile(WWmeta,theshold)
% separate raw data into profiles
% bz, june 15, 2021
load([WWmeta.matpath,WWmeta.name_rbr,'_CTDall.mat']);

if WWmeta.salspiking.yn==1
    % Salinity de-spiking: Find smoothest part of upcasts and fit spectra
    dPdt = diff(CTDall.P)./(diff(CTDall.time)*86400);
    % remove weird pattern in TLC data
    pattern = abs(dPdt(1:end-2))<0.01 & abs(dPdt(2:end-1))<0.01 & dPdt(3:end)<-1;
    patt_s = find(pattern);
    patt_e = patt_s+2;
    pattern_mask = zeros(size(pattern));
    for n = 1:size(patt_s)
        pattern_mask(patt_s(n):patt_e(n)) = 1;
    end
    dPdt(pattern_mask==1) = NaN;

    % Get dPdt, a smoothed version, and the running variance
    dPdt_sm = smoothdata(dPdt,'gaussian',WWmeta.salspiking.win_len);
    dPdt_var = movstd(dPdt,WWmeta.salspiking.win_len,'omitnan');

    % isolate the smooth parts of the upcasts
    data_mask = [0;dPdt_sm]>nanmedian(dPdt(dPdt_sm<0))*0.8...
              | [0;dPdt]>nanmedian(dPdt(dPdt_sm<0))*0.8 ...
              | [0;dPdt_var]>nanmedian(dPdt_var(dPdt_sm<0))*1.5;

    [GainPoly,PhasePoly] = SalSpiking_Coeffs(WWmeta.salspiking.win_len,WWmeta.salspiking.Fs...
        ,CTDall.C,CTDall.T,data_mask,WWmeta.salspiking.iPlot,WWmeta.salspiking.dTorT,WWmeta.salspiking.poly_order);
end

[up,down,dataup,datadown] = get_ctd_2G(CTDall,theshold);  % find upcast/downcast, pay attention to this threshold value

dup=diff(up);
ind_prof=find(dup>1);
ind_prof = [0; ind_prof; length(up)];  % add the first index and the last index
RBRprofiles=struct([]);
fields=fieldnames(dataup);
tdata=dataup.time;
for i=1:length(ind_prof)-1
    for f=1:length(fields)
        wh_field=fields{f};
        if (length(tdata)==length(dataup.(wh_field)))
            RBRprofiles{i}.(wh_field)=dataup.(wh_field)(ind_prof(i)+1:ind_prof(i+1));
%             RBRprofiles{i}.info=CTDall.info;
        end
    end
end

index=find(abs(cellfun(@(x) x.P(end)-x.P(1),RBRprofiles))<5);
RBRprofiles(index) = [];

if WWmeta.salspiking.yn==1
    for n=1:length(RBRprofiles)
        if length(RBRprofiles{n}.P)>WWmeta.salspiking.poly_order*2
            [SP_Corr,SP_sm,SP] = SalSpiking_Profiles(WWmeta.salspiking.Fs,RBRprofiles{n}.C,...
                RBRprofiles{n}.T,RBRprofiles{n}.P,GainPoly,PhasePoly,0,WWmeta.salspiking.dTorT);
            RBRprofiles{n}.SP_corr = SP_Corr;
            RBRprofiles{n}.SP_sm = SP_sm;
            RBRprofiles{n}.SP = SP;
            RBRprofiles{n}.SA_corr = gsw_SA_from_SP(RBRprofiles{n}.SP_corr,RBRprofiles{n}.P,WWmeta.lon,WWmeta.lat);
            RBRprofiles{n}.rho_corr = gsw_rho(RBRprofiles{n}.SA_corr,gsw_CT_from_t(RBRprofiles{n}.SA_corr,RBRprofiles{n}.T,0),0);
        else
            RBRprofiles{n}.SP_corr = nan(size(RBRprofiles{n}.S));
            RBRprofiles{n}.SP_sm = nan(size(RBRprofiles{n}.S));
            RBRprofiles{n}.SP = nan(size(RBRprofiles{n}.S));
            RBRprofiles{n}.SA_corr = nan(size(RBRprofiles{n}.S));
            RBRprofiles{n}.rho_corr = nan(size(RBRprofiles{n}.S));            
        end
    end
end

%%
% calulate thorpe scales from full rate density data
if WWmeta.thorpescales.yn==1
    for n=1:length(RBRprofiles)
        if isfield(RBRprofiles{n},'rho_corr')
            [thorpescale,drdz_Smyth,drdz_Dillon] = f_ThorpeScales(RBRprofiles{n}.rho_corr,RBRprofiles{n}.P,RBRprofiles{n}.time);
        else
            [thorpescale,drdz_Smyth,drdz_Dillon] = f_ThorpeScales(RBRprofiles{n}.rho,RBRprofiles{n}.P,RBRprofiles{n}.time);
        end
        RBRprofiles{n}.ThorpeScale = thorpescale;
        RBRprofiles{n}.drdz_Dillon = drdz_Dillon;
        RBRprofiles{n}.drdz_Smyth = drdz_Smyth;
    end
end
%%
save([WWmeta.propath,WWmeta.name_rbr,'_CTDprofiles.mat'],'RBRprofiles');
end