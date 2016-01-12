function [f,g,h] = LL_boxcox_MATlike(Y,Xa,Xh,EstimOpt,OptimOpt,b0)

% save LL_box_cox_MATlike_tmp
% return
bt_tmp = zeros(sum(EstimOpt.TransActive0 == 1,2),1); ...
bt_tmp(EstimOpt.TransActive(EstimOpt.TransActive0 == 1) == 0) = EstimOpt.NotActive;
bt_tmp(EstimOpt.TransActive(EstimOpt.TransActive0 == 1) == 1) = b0(EstimOpt.NVARA+1:EstimOpt.NVARA+EstimOpt.NVARBCt+EstimOpt.NVARBCh);
b0 = [b0(1:EstimOpt.NVARA); bt_tmp; b0(end)];

bt_tmp = NaN(1 + EstimOpt.NVARA + EstimOpt.NVARH,1); ...
bt_tmp(~isnan(EstimOpt.TransStructNew)) = b0(EstimOpt.TransStructX + EstimOpt.NVARA); ...
b0 = [b0(1:EstimOpt.NVARA); bt_tmp; b0(end)];

LLfun = @(B) LL_boxcox(Y,Xa,Xh,EstimOpt,B);

LL = LLfun(b0);

f = sum(LL);

if isequal(OptimOpt.GradObj,'on')
    if isequal(OptimOpt.Hessian,'user-supplied')
        j = numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central',[ones(1,EstimOpt.NVARA),EstimOpt.TransActive,1]));...
        j(:,[ones(1,EstimOpt.NVARA),EstimOpt.TransActive,1] == 0) = []; ...
        g = sum(j,1)'; ...
        h = j'*j;
    else % OptimOpt.Hessian ~= 'user-supplied'
        j = numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),[ones(1,EstimOpt.NVARA),EstimOpt.TransActive,1]);...
        j(:,[ones(1,EstimOpt.NVARA),EstimOpt.TransActive,1] == 0) = []; ...
        g = sum(j,1)';    
    end 
end