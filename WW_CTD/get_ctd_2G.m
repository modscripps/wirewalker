function [up,down,dataup,datadown] = get_ctd_2G(data,thhold)
% modified by bzheng
% thhold, threshold value to extract fake profiles that are too short
%
% rename variable to make it easy
pdata=double(data.P);
tdata=data.time;

% buid a filter 
dt=median(diff(tdata)); % sampling period
T=tdata(end)-tdata(1);  % length of the record
disp('check if time series is shorter than 3 hours')
if T<3/24  
    warning('time serie is less than 3 hours, very short for data processing, watch out the results')
end

disp('smooth the pressure to define up and down cast')
Nb  = 3; % filter order
fnb = 1/(2*dt); % Nyquist frequency
fc  = 1/200/dt; % 600 dt (give "large scale patern") 
[b,a]= butter(Nb,fc/fnb,'low');
filt_pdata=filtfilt(b,a,pdata);
dfilt_pdata=diff(filt_pdata);
pddfilt_pdata = dfilt_pdata(1:end-1).*dfilt_pdata(2:end);

% dfilt_pdata=diff(pdata);
% pddfilt_pdata = dfilt_pdata(1:end-1).*dfilt_pdata(2:end);

pind = find(pddfilt_pdata <= 0);
csind = [1; pind(1:end)+1];     % long column of cast start time
ceind = [pind(1:end); length(pdata)];   % cast endings
thind = [];
for i = 1:length(csind)
    if ceind(i)-csind(i) < thhold  % too short to be a profile
        thind = [thind i]; 
    end
end
csind(thind) = [];
ceind(thind) = [];
fprintf('identify %i profile \n',round(length(csind)/2))
%%
slope=nan*pdata;
for i=1:length(csind)
    [~,mI]=min(pdata(csind(i):ceind(i)));
    [~,MI]=max(pdata(csind(i):ceind(i)));
    if MI>mI
        slope(csind(i):ceind(i))=0;%downcast
    else
        slope(csind(i):ceind(i))=1;%upcast
    end
end

down = find(slope==0);
up = find(slope==1);


% Plot it
figure(1);clf;
plot(tdata(down),pdata(down),'b.')
hold on
plot(tdata(up),pdata(up),'r*')
set(gca,'ydir','reverse')
title('Verify rising/falling data separation');
ylabel('depth (m)')
set(gca,'xtick',tdata(1):(tdata(end)-tdata(1))/10:tdata(end),'tickdir','in');
datetick('x','dd HH:MM','keepticks');
xlabel('time (DD HH:MM)')
xlim([tdata(1) tdata(end)])

if nargout==4
    fields=fieldnames(data);
    for f=1:length(fields)
        wh_field=fields{f};
        if (length(tdata)==length(data.(wh_field))||(length(tdata)+1)==length(data.(wh_field)))
            switch length(size(data.(wh_field)))
                case 2
                    dataup.(wh_field)=data.(wh_field)(up,:);
                    datadown.(wh_field)=data.(wh_field)(down,:);
                case 3
                    dataup.(wh_field)=data.(wh_field)(up,:,:);
                    datadown.(wh_field)=data.(wh_field)(down,:,:);
                case 4
                    dataup.(wh_field)=data.(wh_field)(up,:,:,:);
                    datadown.(wh_field)=data.(wh_field)(down,:,:,:);
            end
        end
    end
end

end

