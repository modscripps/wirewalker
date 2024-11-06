function WWgrid(WWmeta,zgrid)
% grid ww profile data
% bz, june 15, 2021

load([WWmeta.propath,WWmeta.name_rbr,'_CTDprofiles.mat']);

%get the normal upcast (mean P of the upcast ~ median P of all the mean P)

Prbr=cellfun(@(x) mean(x.P),RBRprofiles);
critp= max(cellfun(@(x) max(x.P),RBRprofiles))-.5*std(Prbr);
critm= min(cellfun(@(x) min(x.P),RBRprofiles))+.5*std(Prbr);

timerbr=cellfun(@(x) mean(x.time),RBRprofiles);
timerbrOK=timerbr(Prbr>critm & Prbr<critp);
indOK=(Prbr>critm & Prbr<critp);
RBRprofiles=RBRprofiles(indOK);

% combine all upcasts
fields=fieldnames(RBRprofiles{1});
idx = size(length(RBRprofiles),2);
for f=1:length(fields)
    wh_field=fields{f};
    if ~strcmp(wh_field,'info')
        RBRgrid.(wh_field)=[];
        for t=1:length(RBRprofiles)
            F=RBRprofiles{t}.(wh_field);
            RBRgrid.(wh_field)=[RBRgrid.(wh_field); F];
            if f == length(fields)
                idx(t,2) = length(RBRgrid.(wh_field));
                idx(t,1) = length(RBRgrid.(wh_field)) - length(F) + 1;
            end
            
        end
    end
    
end
RBRgrid.idx = idx;

%% put gridded product into std_profiles struct
% if nargin==2
    zaxis=0:zgrid:max(cellfun(@(x) max(x.P),RBRprofiles));
% else
%     zaxis=0:.25:max(cellfun(@(x) max(x.P),RBRprofiles));
% end
dep_res = nanmedian(diff(zaxis));

Z=length(zaxis);
fields=fieldnames(RBRprofiles{1});
for f=1:length(fields)
    wh_field=fields{f};
    if ~strcmp(wh_field,'info')
        RBRgrid.std_profiles.(wh_field)=zeros([Z,sum(indOK),2]);
        for t=1:length(timerbrOK)
            F=RBRprofiles{t}.(wh_field);
            P_bin = ceil(RBRprofiles{t}.P/dep_res);
            P_bin(P_bin>Z | P_bin<1) = NaN;
            for n = 1:length(P_bin)
                if ~isnan(P_bin(n)) & ~isnan(F(n))
                    % Thorpe scale gridding keeps the maximum value in any
                    % given vertical bin
                    if strcmp(wh_field,'ThorpeScale') | strcmp(wh_field,'drdz_Dillon') | strcmp(wh_field,'drdz_Smyth')
                        if F(n)>RBRgrid.std_profiles.(wh_field)(P_bin(n),t,1)
                            RBRgrid.std_profiles.(wh_field)(P_bin(n),t,1) = F(n);
                        end
                        RBRgrid.std_profiles.(wh_field)(P_bin(n),t,2) = 1;
                    % Otherwise average value over the length of the
                    % vertical bin.
                    else
                        RBRgrid.std_profiles.(wh_field)(P_bin(n),t,1)=RBRgrid.std_profiles.(wh_field)(P_bin(n),t,1)+F(n);
                        RBRgrid.std_profiles.(wh_field)(P_bin(n),t,2)=RBRgrid.std_profiles.(wh_field)(P_bin(n),t,2)+1;
                    end
                end
            end
        end
        RBRgrid.std_profiles.(wh_field) = RBRgrid.std_profiles.(wh_field)(:,:,1)./RBRgrid.std_profiles.(wh_field)(:,:,2);
    end   
end
RBRgrid.std_profiles.z=zaxis';  % column array
RBRgrid.std_profiles.time=timerbrOK;
% RBRgrid.std_profiles.info=RBRprofiles{1}.info;

%%
save([WWmeta.gridpath,WWmeta.name_rbr,'_CTDgrid.mat'],'RBRgrid');
end