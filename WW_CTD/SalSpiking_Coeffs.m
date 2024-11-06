function [Gain_poly,Phase_poly] = SalSpiking_Coeffs(win_len,Fs,C,T,data_mask,iPlot,dTorT,poly_order)
    
    if mod(length(T),2)==0
        C = C(2:end);
        T = T(2:end);
        data_mask = data_mask(2:end);
    end
    
    T_spec = T;C_spec = C;
    T_spec(data_mask)=NaN;C_spec(data_mask) = NaN;
    
    switch dTorT
        case 0
            % Take Derivitives of temperature and conductivity
            dT = diff(T_spec);
            dC = diff(C_spec);
        case 1
            dT = T_spec;
            dC = C_spec;
    end

    % Get conductivity derivitive-temperature derivitive cross spectra and 
    % temperature autospectra 
    [f_win,ft_cross_win]=CrossSpectra(dC,dT,win_len,0.5,'hamming');
    [~,ft_dT_win]=CrossSpectra(dT,dT,win_len,0.5,'hamming');

    % Calculate gain and phase between conductivy and temperature using
    % averaveraged spectra
    Gain = sqrt(ft_cross_win(1:win_len/2).^2./ft_dT_win(1:win_len/2).^2);
    Phase = angle(ft_cross_win(1:win_len/2));

    for p=2:length(Phase)
        while Phase(p)-Phase(p-1)>pi/32
            Phase(p)=Phase(p)-pi/32;
        end
        while Phase(p)-Phase(p-1)<-pi/32
            Phase(p)=Phase(p)+pi/32;
        end
    end

    % Make Gain an even function of frequency
    GainArray = [flip(Gain(2:end))' Gain'];
    % Make Phase an odd function of frequency
    PhaseArray = [fliplr(-Phase(2:end)') Phase'];
    % Make frequency function to match
    fArray = [flip(-f_win(2:end)) f_win]*Fs;

    % Fit 20th order polynomial to gain and phase arrays
    Phase_poly = polyfit(fArray,PhaseArray,poly_order);
    Gain_poly = polyfit(fArray,GainArray,poly_order);

    if iPlot
        figure;
        subplot(121)
        semilogx(f_win*Fs,Phase(1:win_len/2))
        hold on
        semilogx(f_win*Fs,polyval(Phase_poly,f_win*Fs))
        title('Phase');legend({'Cross spectra phase','Fit'},'location','northwest')
        ylabel('Phase (radians)');xlabel('Frequency (Hz)');

        subplot(122);
        semilogx(f_win*Fs,Gain(1:win_len/2))
        hold on
        semilogx(f_win*Fs,polyval(Gain_poly,f_win*Fs))
        title('Gain');legend({'Cross spectra gain','Fit'},'location','southwest')
        ylabel('Gain');xlabel('Frequency (Hz)');
        set(gcf,'Position',[22         362        1301         443])
    end
end