function [Results, LL] = BOXCOX_loop(INPUT, EstimOpt, k)

if k < 1
    k = 1;
end

Results_loop = cell(k,1);
disp([num2str(1,'%1.0f'), '/', num2str(k, '%1.0f'), ' Box-Cox models'])

Results_loop{1} = BOXCOX(INPUT,EstimOpt);
LL = zeros(k,1);
LL(1) = Results_loop{1}.LL;
if k > 1
    b0 = Results_loop{1}.b0_old;
    for i = 2:k
        EstimOpt.b0 = b0.*unifrnd(0.3,1.8, length(b0),1);
        disp([num2str(i,'%1.0f'), '/', num2str(k, '%1.0f'), ' Box-Cox models'])
        Results_loop{i} = BOXCOX(INPUT,EstimOpt);
        LL(i) = Results_loop{i}.LL;
    end
    ind = find(LL == max(LL));
    Results = Results_loop{ind};
else
    Results = Results_loop{1};
end
    