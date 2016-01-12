function LL = LL_boxcox(Y,Xa,Xh,EstimOpt,b)

% save LL_boxcox_tmp
% return

eps_log = EstimOpt.EPS; % log-transformation limit
b_yt = b(EstimOpt.NVARA+1);
b_xat = b(EstimOpt.NVARA+2:end-EstimOpt.NVARH-1);

% transform variables
if EstimOpt.TransStructNew(1) > 0
    if abs(b_yt) > eps_log;
        tY = (Y.^b_yt-1)/b_yt;
    elseif abs(b_yt) < eps_log
        tY = log(Y);
    end
elseif EstimOpt.TransStruct(1) < 0
    tY = exp(Y*exp(b_yt));
    indx = find(EstimOpt.TransStruct(2:EstimOpt.NVARA+1) == EstimOpt.TransStruct(1));
    b_xat(indx) = exp(b_xat(indx));
else
	tY = Y;
end

ind_xbt = abs(b_xat) >= eps_log & EstimOpt.TransStruct(2:EstimOpt.NVARA+1)'>0; ...
ind_xbl = abs(b_xat) < eps_log & EstimOpt.TransStruct(2:EstimOpt.NVARA+1)'>0; ...
ind_xet = EstimOpt.TransStruct(2:EstimOpt.NVARA+1)' < 0; ...

Xa(:,ind_xbt) = (Xa(:,ind_xbt) .^ (b_xat(ind_xbt,ones(size(Xa,1),1))') - 1) ./b_xat(ind_xbt,ones(size(Xa,1),1))'; ...
Xa(:,ind_xbl) = log(Xa(:,ind_xbl)); ...
Xa(:,ind_xet) = exp((b_xat(ind_xet,ones(size(Xa,1),1))').*Xa(:,ind_xet));

% calculate LL
LL_eps = (tY - Xa*b(1:EstimOpt.NVARA)).^2; ...
if EstimOpt.NVARH == 0
	if isnan(b(EstimOpt.NVARA+1))
        LL = 0.5*log(2*pi) + b(end) + 0.5*LL_eps/exp(2*b(end));
    elseif EstimOpt.TransStruct(1) > 0
        LL = -(b(EstimOpt.NVARA+1)-1)*log(Y) + 0.5*log(2*pi) + b(end) + 0.5*LL_eps/exp(2*b(end));
    elseif EstimOpt.TransStruct(1) < 0
        LL = 0.5*log(2*pi) + b(end) + 0.5*LL_eps/exp(2*b(end))-b(EstimOpt.NVARA+1)-exp(b(EstimOpt.NVARA+1))*Y;
	end
else    
    b_xht = b(EstimOpt.NVARA*2+2:end-1);
    ind_xhbt = abs(b_xht) >= eps_log & EstimOpt.TransStruct(EstimOpt.NVARA+2:end)'>0; ...
    ind_xhbl = abs(b_xht) < eps_log & EstimOpt.TransStruct(EstimOpt.NVARA+2:end)'>0; ...
    ind_xhet = EstimOpt.TransStruct(EstimOpt.NVARA+2:end)' < 0; ...

    Xh(:,ind_xhbt) = ((Xh(:,ind_xhbt) .^ (b_xht(ind_xhbt,ones(size(Xh,1),1))') - 1) ./b_xht(ind_xhbt,ones(size(Xh,1),1))').^2; ...
    Xh(:,ind_xhbl) = log(Xh(:,ind_xhbl)).^2; ...
    Xh(:,ind_xhet) = exp(2*(b_xht(ind_xhet,ones(size(Xh,1),1))').*Xh(:,ind_xhet)); ...
    tXh = prod(Xh,2); ...    
    if isnan(b_yt)
        LL = 0.5*log(2*pi) + b(end) + 0.5*log(tXh) + 0.5*LL_eps./(exp(2*b(end))*tXh);
    elseif EstimOpt.TransStruct(1) > 0
        LL = -(b_yt-1)*log(Y) + 0.5*log(2*pi) + b(end) + 0.5*log(tXh) + 0.5*LL_eps./(exp(2*b(end))*tXh);
    elseif EstimOpt.TransStruct(1) < 0
        LL = 0.5*log(2*pi) + b(end) +0.5*log(tXh) + 0.5*LL_eps./(exp(2*b(end))*tXh)-b(EstimOpt.NVARA+1)-exp(b_yt)*Y;
    end
end




