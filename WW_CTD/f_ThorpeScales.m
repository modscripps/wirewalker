% Calulates thorpe scales using wirewalker density profiles, as well as
% mean density gradient in each overturn following Dillon 1982 (Vertical
% Overturns' A Comparison of Thorpe and Ozmidov Length Scales) and Smyth
% et. al 2001 (The Efficiency of Mixing in Turbulent Patches: Inferences
% from Direct Simulations and Microstructure Observations).

% Density profiles are cleaned up following following Gargett and Garner
% 2008 (Determining Thorpe Scales from Ship-Lowered CTD Density Profiles)

function [thorpescale_out,drdz_Smyth_out,drdz_Dillon_out] = f_ThorpeScales(rho,pres,time)
%%
    % take center difference to find wirewalker direction reversals
    dp = diffdiff(pres,1);
    dt = diffdiff(time,1);
    mask1 = find(dp > 0);
    
    % use direction reversals to eliminate data between a direction reversal
    % and the point where the wirewalker returns to the same pressure
    % cooordinate as before the direction reversal

    % if no direction reversals, continue
    if isempty(mask1)
        next = length(pres);
    else
        next = mask1(1);
    end
    % create list of sections to cut from the data
    cut=[];
    while next < length(pres)
        % find pressures less than presure before direction reversal
        biggers = find(pres < pres(next));
        biggers(biggers <= next)=[];
        % If there are no pressures less (if reversal happens at the top of
        % a cast - common) eliminate data till end of record
        if isempty(biggers)
            biggers = length(pres);
        end
        % otherwise, eliminate data from pressure reversal to the point
        % were the pressure is less than before the reversal
        cut(end+1,:) = [next,biggers(1)];
        % find next pressure reversal
        next = min(mask1(mask1>biggers(1)));
    end            

    % remove pressure reversal data
    for i = 1:size(cut,1)
        rho(cut(i,1):cut(i,2)) = 9999;
        pres(cut(i,1):cut(i,2)) = 9999;
        time(cut(i,1):cut(i,2)) = 9999;
        dp(cut(i,1):cut(i,2)) = 9999;
        dt(cut(i,1):cut(i,2)) = 9999;
    end
    nanmask = rho==9999;
    rho(nanmask) = [];
    pres(nanmask) = [];
    time(nanmask) = [];
    dp(nanmask(2:end-1)) = [];
    dt(nanmask(2:end-1)) = [];

    % calculate intermediate density profiles to eliminate false overturs
    % due to sensor noise. Intermediate density keeps only density changes
    % above the noise threshold: otherwise the previous value of density it
    % used. This is calculated both upgoing and downgoing. (Gargett and Garner
    % 2008 (Determining Thorpe Scales from Ship-Lowered CTD Density Profiles))
    intup = zeros(size(rho));
    intdown = zeros(size(rho));
    for i = 1:length(rho)
       % set noise threshhold (in density units - kgm^-3)
       nl = 0.0005*2;
       % Calculate upgoing intermediate density
       if i~=1 && i~=length(rho) && abs(intup(i-1)-rho(i))<nl
           intup(i) = intup(i-1);
       elseif i == 1 || i == length(rho)
           intup(i) = rho(i);
       else
           intup(i)=rho(i);
       end
       % Calculate downgoing intermediate density
       if  i~=length(rho) && i~=1 && abs(intdown(end-i+2)-rho(end-i+1))<nl
           intdown(end-i+1) = intdown(end-i+2);
       elseif i == length(rho) || i==1
           intdown(end-i+1) = rho(end-i+1);
       else
           intdown(end-i+1)=rho(end-i+1);
       end
    end

    % Average the two to get processed (de-noised) density
    rho_proc = (intup+intdown)/2;

    % sort rho to get stable profile
    [sort_rho,is] = sort(rho_proc,'descend');
    % using time as a vertical variable, calculate how far a given
    % isopycnal is displaced from the stable state. This wil later be
    % converted to a vertial distance using the average vertical speed over
    % the overturn
    t_diff = time - time(is);

    % going from the bottom add density displacements from stability. Zero
    % crossings in this function indicate the edges of overturns
    rhodiff = rho_proc-sort_rho;
    rho_diff_sums = zeros(size(rhodiff));
    for i = 1:length(rhodiff)
        rho_diff_sums(i) = sum(rhodiff(1:i));
    end

    % number overturns in a given section
    turnnum = zeros(size(rhodiff));index = 1;reset = 1;i = 1;
    while i<=length(rhodiff)
        if rho_diff_sums(i) ~= 0 && reset == 0
            turnnum(i) = index;
        elseif rho_diff_sums(i) ~= 0 && reset == 1
            index = index + 1;
            reset = 0;
            turnnum(i) = index;
        elseif rho_diff_sums(i) == 0
            reset = 1;
        end
        i = i+1;
    end

    % calculate the thorpe scale in pressure units by multiplying the time
    % displacement by the average dp/dt over the overturn
    thorpescale = zeros(size(turnnum));drdz_Dillon = zeros(size(turnnum));
    drdz_Smyth = zeros(size(turnnum));
    for i = 1:max(turnnum)
        % keep only overturns with similar numbers of positive and negative
        % displacements, otherwise the overturn is likely the result of error.
        if sum(t_diff(turnnum==i)>0)/sum(t_diff(turnnum==i)<0) > 0.2 && ...
                sum(t_diff(turnnum==i)<0)/sum(t_diff(turnnum==i)>0) > 0.2 ...
                && sum(turnnum==i)>2

            spd = nanmean(dp(turnnum(2:end-1)==i)./dt(turnnum(2:end-1)==i));

            thorpescale(turnnum==i) = ones(size(turnnum(turnnum==i)))*...
                sqrt(sum(t_diff(turnnum==i).^2)/sum(turnnum==i)) * abs(spd);

            drdz_Dillon(turnnum==i) = ones(size(turnnum(turnnum==i)))*...
                mean(diff(sort_rho(turnnum==i))./diff(time(turnnum==i))) / spd;

            t_offset = t_diff(turnnum==i);
            drdz_Smyth(turnnum==i) = sqrt(mean((diffdiff(sort_rho(turnnum==i),1)./diffdiff(time(turnnum==i),1)/spd).^2.*...
                                     (t_offset(2:end-1)*spd).^2,'omitnan')/mean((t_offset(2:end-1)*spd).^2,'omitnan'));

        else
            thorpescale(turnnum==i) = zeros(size(turnnum(turnnum==i)));
            drdz_Smyth(turnnum==i) = zeros(size(turnnum(turnnum==i)));
            drdz_Dillon(turnnum==i) = zeros(size(turnnum(turnnum==i)));
        end
    end
    thorpescale_out = zeros(size(nanmask));
    drdz_Smyth_out = zeros(size(nanmask));
    drdz_Dillon_out = zeros(size(nanmask));

    thorpescale_out(~nanmask) = thorpescale;
    drdz_Smyth_out(~nanmask) = drdz_Smyth;
    drdz_Dillon_out(~nanmask) = drdz_Dillon;
    
end