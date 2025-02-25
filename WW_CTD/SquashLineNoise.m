function [out]=SquashLineNoise(in,win_size)
%     win_size = 1000;
    
    % FFT input data
    [f,num_wins,fts] = fftWindowed(in,1000,0.5,'hamming');
    
    % identify narrow peaks
    norm_spec = log10(fts) - log10(smoothdata(fts,'movmedian',round(win_size/30)));
    [pks,locs] = findpeaks(norm_spec);
    
    fs_ind = locs(pks>0.25 & transpose(f(locs))*8>1);
    
    % select peaks
    fs = f(fs_ind)*8;
    
    % Plot raw spectra and selected peaks
    figure;
    p1=loglog(f*8,fts);
    hold on
    for n = 1:length(fs)
        [~,i]=min(abs(fs(n)-8*f));
        scatter(fs(n),fts(i),30,Color('c1'));
    end
    
    % generate frequency vetor for full timeseries
    f_all = 8*(0:(length(in)/2))/length(in);
    % specify filter bandwidth
    filt_bw = 0.03;
    
    % Generate filter
    filter = ones(size(f_all));
    for n = 1:length(fs)
        filt_mask = f_all>fs(n)-fs(n)*filt_bw/2 & f_all<fs(n)+fs(n)*filt_bw/2;
        filter(filt_mask)=0;
    end
    % Smooth filtert
    filter = smoothdata(filter,'gaussian',length(in)*filt_bw/5);
    
    filter_2 = transpose([filter,fliplr(filter)]);
    
    %Apply filter
    out = real(ifft(fft(in).*filter_2(1:length(in))));
    
    %% Plot filtered data
    [f,num_wins,fts_post] = fftWindowed(out,1000,0.5,'hamming');
    
    hold on;
    p2=loglog(f*8,nanmedian(fts_post,2));
    xlabel('Frequency (Hz)');
    ylabel('^{\circ}C^2/Hz');
    legend({'Pre-Filter','Post-Filter'})
    
    % 
    % figure;
    % idx(1) = 11140100;idx(2) = 11144100;
    % plot(CTDall.T(idx(1):idx(2)),CTDall.P(idx(1):idx(2)))
    % hold on
    % plot(T_use(idx(1):idx(2)),CTDall.P(idx(1):idx(2)))
end