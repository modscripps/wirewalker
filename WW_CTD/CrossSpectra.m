function [f,xy_list,num_wins]=CrossSpectra(data1,data2,win_length,overlap,window)
%%%
%     This function calculates cross spectra using an overlapping window method 
%     with an arbitrary overlap specified by "overlap". It takes 5 inputs
%     
%     data1 & 2: 2 same length vectors of data to be fft'd
%     win_length: the lenth of the overlaping windows into which the data ist to 
%     be subdivided
%     overlap: A number between zero and 1 specifying the fractional overlap 
%     between windows
%     window: the window function to be used, either 'hanning','hamming',or 
%     'boxcar'
%     
%     outputs are a sum of pwer spectra and the number of windows used. find the 
%     averaged power spectra as fft_list/num_wins
%%%

% set window
if strcmp(window,'hanning')
        window=hanning(win_length);
elseif strcmp(window,'hamming')
        window=hamming(win_length);
elseif strcmp(window,'boxcar')
        window=ones(win_length)/sqrt(win_length);
end

% initalize list of ffts
xy_list=zeros(win_length,1);
% set initail center-point
en=int32(ceil(win_length*(overlap)));
% initalize counter for number of fft's preformed
num_wins=0;
while int32(ceil(en+win_length*(1-overlap)))<length(data1)
    % get data in an overlaped window
    data1_use=squeeze(data1(int32(ceil(en-win_length*(overlap)))+1:int32(ceil(en+win_length*(1-overlap)))));
    data2_use=squeeze(data2(int32(ceil(en-win_length*(overlap)))+1:int32(ceil(en+win_length*(1-overlap)))));
    % check for completness of the dataset
    if ((sum(isnan(data1_use))/win_length < 0.1) && (sum(isnan(data2_use))/win_length < 0.1))
        % fill missing values with the mean
%         data1_use(isnan(data1_use)) = nanmean(data1_use);
%         data2_use(isnan(data2_use)) = nanmean(data2_use);
        % Detrend data
        data1_use=detrend(data1_use,'omitnan');
        data2_use=detrend(data2_use,'omitnan');
        data1_use=fillmissing(data1_use,'linear','EndValues',0);
        data2_use=fillmissing(data2_use,'linear','EndValues',0);
        % take fft, normalize, and add to list
        ft1=fft(data1_use.*window);
        ft2=fft(data2_use.*window);

        xy_list = xy_list+(ft1.*conj(ft2)*(win_length/sum(window.^2)));

        % increase counter by 1
        num_wins=num_wins+1;
    end
    
    % set next centerpoint value
    en=int32(en+ceil(win_length*(1-overlap)));

end

xy_list = (xy_list/num_wins)*2/(numel(window)^2);

f = (0:(win_length/2)-1)/win_length;
end