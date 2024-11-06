function [SP_Corr,SP_sm,SP] = SalSpiking_Profiles(Fs,C,T,P,Gain_poly,Phase_poly,iPlot,dTorT)
    
    even_check = mod(length(T),2)==0;
    
    if even_check
        C = C(2:end);
        T = T(2:end);
        P = P(2:end);
    end
 
    % get frequency's for fft of a full profile
    switch dTorT
        case 0
            f = Fs*(0:((3*length(C)/2)))/(3*length(C));
            f2 = [f,-fliplr(f(2:end-1))];
        case 1
            f = Fs*(0:3*length(C)/2)/(3*length(C));
            f2 = [f,-fliplr(f(2:end))];
    end

    % Get values of the gain and phase polynomials for each frequency in 
    % the profile fft
    GainFit = polyval(Gain_poly,f2);
    PhaseFit = polyval(Phase_poly,f2);
    GainFit=GainFit./GainFit(1);
    
    nanmask = isnan(T) | isnan(C) | isnan(P);
    T = fillmissing(T,'linear');
    C = fillmissing(C,'linear');
    P = fillmissing(P,'linear');
    
    T_to_corr = [flipud(T);T;flipud(T)];
    
    switch dTorT
        case 0
            % get temperature profile fft
            ft_T = fft(diff(T_to_corr));
            % correct temperature profile fft
            ft_T_corr=ft_T.*GainFit'.*exp(1i.*PhaseFit');
            % get time domain profile back
            T_Corr = [T_to_corr(1);cumsum(real(ifft(ft_T_corr)))];
        case 1
            % get temperature profile fft
            ft_T = fft(T_to_corr);
            % correct temperature profile fft
            ft_T_corr=ft_T.*GainFit'.*exp(1i.*PhaseFit');
            % get time domain profile back
            T_Corr = real(ifft(ft_T_corr));
    end
    
    T_Corr = T_Corr(length(T)+1:2*length(T));
    
    T_Corr = T_Corr - mean(T_Corr,'omitnan');
    T_Corr = T_Corr + mean(T,'omitnan');
    
    % smooth profiles at 1/3 nyquist frequency
    C = smoothdata(C,'movmedian',3);
    [b,a] = butter(4,1/3);
    T_Corr_sm = filtfilt(b,a,T_Corr);
    T_sm = filtfilt(b,a,T);
    C_sm = filtfilt(b,a,C);
    P_sm = filtfilt(b,a,P);
        
    % calulate salinity
    SP_Corr = gsw_SP_from_C(C_sm,T_Corr_sm,P_sm);
    SP_sm = gsw_SP_from_C(C_sm,T_sm,P_sm);
    SP = gsw_SP_from_C(C,T,P);
    
    SP_Corr(nanmask) = NaN;
    SP_sm(nanmask) = NaN;
    SP(nanmask) = NaN;
    
    if even_check
        SP_Corr = [SP_Corr(1);SP_Corr];
        SP_sm = [SP_sm(1);SP_sm];
        SP = [SP(1);SP];
    end

    if iPlot
        figure;plot(SP)
        hold on;plot(SP_sm)
        hold on;plot(SP_Corr)
    end
end