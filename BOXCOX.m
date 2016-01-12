function Results = BOXCOX(INPUT, varargin) % INPUT,Results_old,EstimOpt

% save BOXCOX_tmp
% return

tic

Results.R_out = {};
Results.R = [];
Results.bhat = [];
Results.stats = [];

format shortG;
format compact;


%% check input:


% check and assign input arguments
if nargin == 1 
    EstimOpt = [];
elseif nargin == 2
    EstimOpt = varargin{1};
else
    error('Not enough input arguments')
end

disp(' ');
disp('__________________________________________________________________________________________________________________');
disp(' ');
disp('Box-Cox-Exp regression model ...')

% check heteroskedascity settings
if ~isfield(INPUT, 'Xh')
    INPUT.Xh = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVARH = size(INPUT.Xh,2);

if EstimOpt.NVARH == 0
    disp('with homoskedascity.')
else
    disp('with heteroskedascity.')
end


% check EstimOpt
if ~isfield(EstimOpt, 'EPS') || isempty(EstimOpt.EPS)
    EstimOpt.EPS = 1.e-12; % overall precision level
end
if ~isfield(EstimOpt, 'HessEstFix') || isempty(EstimOpt.HessEstFix)
    EstimOpt.HessEstFix = 1;
end
if ~isfield(EstimOpt, 'NP') || isempty(EstimOpt.NP)
    EstimOpt.NP = size(INPUT.Y,1);
end
if ~isfield(EstimOpt, 'NVARA') || isempty(EstimOpt.NVARA)
    EstimOpt.NVARA = size(INPUT.Xa,2);
end
if ~isfield(EstimOpt, 'spacing') || isempty(EstimOpt.spacing)
    EstimOpt.spacing = 4;
end
if ~isfield(EstimOpt, 'precision') || isempty(EstimOpt.precision)
    EstimOpt.precision = 4;
end

% set OptimOpt
OptimOpt = optimoptions('fminunc');
if isfield(EstimOpt, 'Algorithm') && ~isempty(EstimOpt.Algorithm)
    OptimOpt.Algorithm = EstimOpt.Algorithm;
else
    OptimOpt.Algorithm = 'quasi-newton'; %'trust-region';
    %OptimOpt.Algorithm = 'trust-region';
end
if isfield(EstimOpt, 'FinDiffType') && ~isempty(EstimOpt.FinDiffType)
    OptimOpt.FinDiffType = EstimOpt.FinDiffType;
else
    OptimOpt.FinDiffType = 'forward'; %central
end
if isfield(EstimOpt, 'GradObj') && ~isempty(EstimOpt.GradObj)
    OptimOpt.GradObj = EstimOpt.GradObj;
else
    OptimOpt.GradObj =  'off'; %'off';
end
if isfield(EstimOpt, 'Hessian') && ~isempty(EstimOpt.Hessian)
    if isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(EstimOpt.Hessian,'user-supplied')
        cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway.\n')
        OptimOpt.Hessian = 'off';
    elseif isequal(EstimOpt.Hessian,'off')
        OptimOpt.Hessian = EstimOpt.Hessian;
    elseif isequal(OptimOpt.Algorithm,'trust-region')
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect Hessian specification - switching to ''user-supplied''\n')
        OptimOpt.Hessian = 'user-supplied';
    else %if isequal(OptimOpt.Algorithm,'quasi-newton')
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect Hessian specification - switching to ''off''\n')
        OptimOpt.Hessian = 'off';
    end
elseif isequal(OptimOpt.Algorithm,'trust-region')
    OptimOpt.Hessian = 'user-supplied';
else % isequal(OptimOpt.Algorithm,'quasi-newton')
	OptimOpt.Hessian = 'off';
end
if isfield(EstimOpt, 'MaxIter') && ~isempty(EstimOpt.MaxIter)
    OptimOpt.MaxIter = EstimOpt.MaxIter;
else
    OptimOpt.MaxIter =  1e6; % Maximum number of iterations allowed
end
if isfield(EstimOpt, 'MaxFunEvals') && ~isempty(EstimOpt.MaxFunEvals)
    OptimOpt.MaxFunEvals = EstimOpt.MaxFunEvals;
else
    OptimOpt.MaxFunEvals =  1e8; % Maximum number of function evaluations allowed (1000)
end
if isfield(EstimOpt,'TolFun') && ~isempty(EstimOpt.TolFun)
    OptimOpt.TolFun = EstimOpt.TolFun; % df / gradient precision level
elseif isfield(EstimOpt,'EPS')
    OptimOpt.TolFun = EstimOpt.EPS;
end
if isfield(EstimOpt,'TolX') && ~isempty(EstimOpt.TolX)
    OptimOpt.TolX = EstimOpt.TolX; % step precision level
elseif isfield(EstimOpt,'EPS')
    OptimOpt.TolX = EstimOpt.EPS;
end
if isfield(EstimOpt,'FunValCheck') && ~isempty(EstimOpt.FunValCheck)
    OptimOpt.FunValCheck = EstimOpt.FunValCheck;
else
    OptimOpt.FunValCheck = 'on';
end
if isfield(EstimOpt,'Diagnostics') && ~isempty(EstimOpt.Diagnostics)
    OptimOpt.Diagnostics = EstimOpt.Diagnostics;
else
    OptimOpt.Diagnostics = 'on';
end
OptimOpt.OutputFcn = @outputf;



% check transformation settings

if ~isfield(EstimOpt, 'TransStruct') || isempty(EstimOpt.TransStruct)
    cprintf(rgb('DarkOrange'), 'WARNING: Transformation structure not specified - assuming no transformation\n');
	EstimOpt.TransStruct = NaN(1, EstimOpt.NVARA+1);
	if EstimOpt.NVARH > 0
%         cprintf(rgb('DarkOrange'), 'WARNING: Assuming exponential transformation for heteroskedascity\n');
        EstimOpt.TransStruct = [EstimOpt.TransStruct, -(1:EstimOpt.NVARH)];
	end
elseif length(EstimOpt.TransStruct) ~= (1 + EstimOpt.NVARA + EstimOpt.NVARH) && length(EstimOpt.TransStruct) ~= (1 + EstimOpt.NVARA)
    error('Incorrect length of transformation parameter indicators (EstimOpt.TransStruct)')
else
    EstimOpt.TransStruct(EstimOpt.TransStruct == 0) = NaN;
    if length(EstimOpt.TransStruct) == (1 + EstimOpt.NVARA) && EstimOpt.NVARH > 0
        if any(EstimOpt.TransStruct < 0)
            ind_tmp = min(EstimOpt.TransStruct(EstimOpt.TransStruct < 0)) - 1;
        else
            ind_tmp = -1;
        end
        EstimOpt.TransStruct = [EstimOpt.TransStruct, ind_tmp:-1:(ind_tmp - EstimOpt.NVARH + 1)];
        cprintf(rgb('DarkOrange'), 'WARNING: Assuming exponential transformation for heteroskedascity\n');
    elseif length(EstimOpt.TransStruct) == (1 + EstimOpt.NVARA + EstimOpt.NVARH) && EstimOpt.NVARH > 0 && any(isnan(EstimOpt.TransStruct(EstimOpt.NVARA+2:end)))
        error('Model with heteroskedascity requires specifying transformation parameters structure')
    end
    EstimOpt.TransStruct = EstimOpt.TransStruct(:)';
end

% replace transformation parameter names
EstimOpt.TransStructNew = NaN(size(EstimOpt.TransStruct));
kb = 1;
kb_used = [];
ke = -1;
ke_used = [];
for i = 1:length(EstimOpt.TransStruct)
	if ~isnan(EstimOpt.TransStruct(i))
        if EstimOpt.TransStruct(i) > 0 && ~any(EstimOpt.TransStruct(i) == kb_used)            
            if i == 1
                EstimOpt.TransStructNew(EstimOpt.TransStruct == EstimOpt.TransStruct(i)) = Inf; 
            else
                EstimOpt.TransStructNew(EstimOpt.TransStruct == EstimOpt.TransStruct(i)) = kb;
                kb = kb + 1;
            end
            kb_used = [kb_used, EstimOpt.TransStruct(i)];            
        elseif EstimOpt.TransStruct(i) < 0 && ~any(EstimOpt.TransStruct(i) == ke_used)
            if i == 1
                EstimOpt.TransStructNew(EstimOpt.TransStruct == EstimOpt.TransStruct(i)) = -Inf;
            else                
                EstimOpt.TransStructNew(EstimOpt.TransStruct == EstimOpt.TransStruct(i)) = ke;
                ke = ke - 1;
            end
            ke_used = [ke_used, EstimOpt.TransStruct(i)];
        end         
	end
end

k = 1;
kb_used = [];
ke_used = [];
EstimOpt.TransStructX = zeros(1, sum(~isnan(EstimOpt.TransStructNew),2));
for i = 1:length(EstimOpt.TransStruct)
	if ~isnan(EstimOpt.TransStruct(i))
        if EstimOpt.TransStruct(i) > 0 && ~any(EstimOpt.TransStruct(i) == kb_used)            
            EstimOpt.TransStructX(EstimOpt.TransStruct(~isnan(EstimOpt.TransStruct)) == EstimOpt.TransStruct(i)) = k;
            k = k + 1;
            kb_used = [kb_used, EstimOpt.TransStruct(i)];
        elseif EstimOpt.TransStruct(i) < 0 && ~any(EstimOpt.TransStruct(i) == ke_used)
            ke_used = [ke_used, EstimOpt.TransStruct(i)];
            EstimOpt.TransStructX(EstimOpt.TransStruct(~isnan(EstimOpt.TransStruct)) == EstimOpt.TransStruct(i)) = k;
            k = k + 1;
        end
	end
end

EstimOpt.TransUnique = unique(EstimOpt.TransStructNew(~isnan(EstimOpt.TransStructNew)));
EstimOpt.TransActive0 = ones(size(EstimOpt.TransStruct));
EstimOpt.TransActive0(isnan(EstimOpt.TransStruct)) = 0;
EstimOpt.TransActive0(EstimOpt.TransStruct == 0) = 0;

for i = 1:length(EstimOpt.TransUnique)
    idx = (EstimOpt.TransStructNew == EstimOpt.TransUnique(i)); ...
    idx(find(idx,1,'first')) = 0; ...
    EstimOpt.TransActive0(idx) = 0;  
end

if isfield(EstimOpt, 'TransActive') && ~isempty(EstimOpt.TransActive)
    EstimOpt.TransActive = EstimOpt.TransActive(:)';
    if length(EstimOpt.TransActive) ~= length(EstimOpt.TransStruct)
        error('Incorrect length of transformation parameters restrictions (EstimOpt.TransActive)')
    elseif any(EstimOpt.TransActive - EstimOpt.TransActive0 > 0)
        cprintf(rgb('DarkOrange'), 'WARNING: Transformation parameters restrictions incorrect - ignoring user input\n');
        EstimOpt.TransActive = EstimOpt.TransActive0;
    end
else
    EstimOpt.TransActive = EstimOpt.TransActive0;
end

if isfield(EstimOpt, 'NotActive') == 0 || numel(EstimOpt.NotActive) == 0 || length(EstimOpt.NotActive) ~= sum(EstimOpt.TransActive == 0 & EstimOpt.TransActive0 == 1,2)
    EstimOpt.NotActive = zeros(1,0);
else
    EstimOpt.NotActive = EstimOpt.NotActive(:)';
end

% check variable ranges ...
if sum([any(INPUT.Xa(:, EstimOpt.TransStruct(2:EstimOpt.NVARA+1) > 0) <= 0) , any(INPUT.Y) <= 0 && EstimOpt.TransStruct(1) > 0]) > 0
 	error('Box-Cox transformed variables must be strictly greater than zero')
end
if EstimOpt.NVARH > 0 && sum(any(min(INPUT.Xh(:,EstimOpt.TransStruct(EstimOpt.NVARA+2:end) > 0)) <= 1),2) > 0
    cprintf(rgb('DarkOrange'), 'WARNING: Box-Cox transformed heteroskedascity variables must be greater than one - adding a constant\n');
    Min_tmp = min(INPUT.Xh(:,EstimOpt.TransStruct(EstimOpt.NVARA+2:end) > 0));
    INPUT.Xh(:,EstimOpt.TransStruct(EstimOpt.NVARA+2:end) > 0) = INPUT.Xh(:,EstimOpt.TransStruct(EstimOpt.NVARA+2:end) > 0) + (1.01 - Min_tmp(ones(size(INPUT.Xh,1),1),:));
end
% if EstimOpt.NVARH > 0 && any(EstimOpt.TransStruct(EstimOpt.NVARA+2:end) > 0)
%     INPUT.Xh(:,EstimOpt.TransStruct(EstimOpt.NVARA+2:end) > 0) = INPUT.Xh(:,EstimOpt.TransStruct(EstimOpt.NVARA+2:end) > 0).^2;
% end

% check names 
if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= EstimOpt.NVARA
    EstimOpt.NamesA = (1:EstimOpt.NVARA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
elseif size(EstimOpt.NamesA,1) ~= EstimOpt.NVARA
    EstimOpt.NamesA = EstimOpt.NamesA';
end

if isfield(EstimOpt,'NamesY') == 0 || isempty(EstimOpt.NamesY) || length(EstimOpt.NamesY) ~= 1
    EstimOpt.NamesY = cellstr('Y');
end

if EstimOpt.NVARH > 0
    if isfield(EstimOpt,'NamesH') == 0 || isempty(EstimOpt.NamesH) || length(EstimOpt.NamesH) ~= EstimOpt.NVARH
        EstimOpt.NamesH = (1:EstimOpt.NVARH)';
        EstimOpt.NamesH = cellstr(num2str(EstimOpt.NamesH));
    elseif size(EstimOpt.NamesH,1) ~= (EstimOpt.NVARH)
        EstimOpt.NamesH = EstimOpt.NamesH';
    end
end

for n = 1:1+EstimOpt.NVARA+EstimOpt.NVARH
    if isnan(EstimOpt.TransStruct(n))
        EstimOpt.NamesTS{n,1} = NaN;
    elseif EstimOpt.TransStructNew(n) > 0
        if EstimOpt.TransStructNew(n) == Inf
            EstimOpt.NamesTS{n,1} = ['Lambda' num2str(0)];
        else
            EstimOpt.NamesTS{n,1} = ['Lambda' num2str(EstimOpt.TransStructNew(n))];
        end
    elseif EstimOpt.TransStructNew(n) < 0
        if EstimOpt.TransStructNew(n) == -Inf
            EstimOpt.NamesTS{n,1} = ['Theta' num2str(0)];
        else
            EstimOpt.NamesTS{n,1} = ['Theta' num2str(-EstimOpt.TransStructNew(n))];
        end
    end
end

EstimOpt.NVARBCt = sum(EstimOpt.TransActive(1:EstimOpt.NVARA+1) == 1);
EstimOpt.NVARBCh = sum(EstimOpt.TransActive(EstimOpt.NVARA+2:end) == 1);

warning off MATLAB:mir_warning_maybe_uninitialized_temporary


%% starting values 


% if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVARA + EstimOpt.NVARBCt + 1 + EstimOpt.NVARBCh
%     b0 = B_backup(:);
%     disp('Using the starting values from Backup')
% elseif ~isempty('EstimOpt') && isfield(EstimOpt,'b0') && ~isempty(EstimOpt.b0) % starting values provided
if ~isempty('EstimOpt') && isfield(EstimOpt,'b0') && ~isempty(EstimOpt.b0) % starting values provided
    if length(EstimOpt.b0) ~= EstimOpt.NVARA + EstimOpt.NVARBCt + 1 + EstimOpt.NVARBCh
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification - ignoring user starting values\n');
    else
        b0 = EstimOpt.b0(:);
    end
end
if  ~exist('b0','var')
    disp('Using linear regression estimates as starting values')
    [b0, ~, ~, ~, stats] = regress(INPUT.Y,INPUT.Xa);
    
%     b0 = [b0; ones(EstimOpt.NVARBCt,1); zeros(EstimOpt.NVARBCh,1); log((stats(4))^0.5)]; 
%     b0 = [b0; exp(1)*ones(EstimOpt.NVARBCt,1); ones(EstimOpt.NVARBCh,1); log((stats(4))^0.5)];
%     btmp1 = 0.01 .* ones(EstimOpt.NVARH,1); % linear coefficients for heteroskedascity parameters
    btmp2 = ones(EstimOpt.NVARBCt,1); % transformation parameters for NVARA
    btmp2(EstimOpt.TransStructNew(EstimOpt.TransActive(1:EstimOpt.NVARA+1) == 1) < 0) = 0; % exponentiation parameters    
    btmp3 = ones(EstimOpt.NVARBCh,1); % transformation parameters for NVARH
    btmp3(EstimOpt.TransStructNew(EstimOpt.NVARA + 1 + EstimOpt.TransActive(EstimOpt.NVARA+1:end) == 1) < 0) = 0; % exponentiation parameters
    b0 = [b0; btmp2; btmp3; log((stats(4))^0.5)];

%     ones(EstimOpt.NVARBCt,1); % transformation parameters for NVARA
%     btmp1(EstimOpt.TransStructNew(EstimOpt.TransActive(1:EstimOpt.NVARA+1) == 1) < 0) = 0; % exponentiation parameters
%     
%     .* EstimOpt.NVARBCt
%         ; zeros(EstimOpt.NVARBCh,1)
%     l = 1;
%     btmp = ones(EstimOpt.NVARBCt + EstimOpt.NVARBCh,1);
%     for i = 1:length(EstimOpt.TransStruct)
%         if EstimOpt.TransActive(i) == 1
%             if EstimOpt.TransStruct(i) < 0 && i == 1
%                 btmp(l) = log(0.1); 
%             elseif EstimOpt.TransStruct(i) < 0
%                 btmp(l) = 0.1;
%             end
%             l = l+1;
%         end
%     end
%     b0 = [b0; btmp; log((stats(4))^0.5)]; 
%     b0 = [b0; exp(1)*ones(EstimOpt.NVARBCt,1); ones(EstimOpt.NVARBCh,1); log((stats(4))^0.5)]; 
end
% EstimOpt.Active_ind = find(EstimOpt.TransActive0 == 1 && EstimOpt.TransActive == 0);
% EstimOpt.TransActive0(EstimOpt.Active_ind) = b0(EstimOpt.NVARA + EstimOpt.Active_ind);

%% Estimating LL0


OptimOpt_0 = optimoptions('fminunc');
OptimOpt_0.Algorithm = 'quasi-newton';
OptimOpt_0.GradObj = 'off';
OptimOpt_0.Hessian = 'off';
OptimOpt_0.Display = 'off';
OptimOpt_0.FunValCheck = 'off'; 
OptimOpt_0.Diagnostics = 'off';

LLfun0 = @(B) EstimOpt.NP*(0.5*log(2*pi) + B(2)) + 0.5*sum(((INPUT.Y - B(1)).^2)/exp(2*B(2)),1);
       
[b00, ~, ~, ~, stats0] = regress(INPUT.Y,ones(size(INPUT.Y)));
[Results.bhat0, LL0] = fminunc(LLfun0, [b00;log((stats0(4))^0.5)],OptimOpt_0);

Results.LL0 = -LL0;


%% estimation

disp(' ')
LLfun = @(B) LL_boxcox_MATlike(INPUT.Y,INPUT.Xa,INPUT.Xh,EstimOpt,OptimOpt,B);
if EstimOpt.HessEstFix == 0
    [Results.bhat,LL,Results.exitf,Results.output,Results.g,Results.hess] = fminunc(LLfun,b0,OptimOpt);
else
    [Results.bhat,LL,Results.exitf,Results.output,Results.g] = fminunc(LLfun,b0,OptimOpt);
end      

% save est_out


%% Output


Results.LL = -LL;
Results.b0_old = b0;

bt_tmp = zeros(sum(EstimOpt.TransActive0 == 1,2),1); ...
bt_tmp(EstimOpt.TransActive(EstimOpt.TransActive0 == 1) == 0) = EstimOpt.NotActive;
bt_tmp(EstimOpt.TransActive(EstimOpt.TransActive0 == 1) == 1) = Results.bhat(EstimOpt.NVARA+1:EstimOpt.NVARA+EstimOpt.NVARBCt+EstimOpt.NVARBCh);
Results.bhat = [Results.bhat(1:EstimOpt.NVARA); bt_tmp; Results.bhat(end)];

bt_tmp = NaN(1 + EstimOpt.NVARA + EstimOpt.NVARH,1); ...
bt_tmp(~isnan(EstimOpt.TransStructNew)) = Results.bhat(EstimOpt.TransStructX + EstimOpt.NVARA); ...
Results.bhat_full = [Results.bhat(1:EstimOpt.NVARA); bt_tmp; Results.bhat(end)];

if EstimOpt.HessEstFix == 1
	f = LL_boxcox(INPUT.Y,INPUT.Xa,INPUT.Xh,EstimOpt,Results.bhat_full);...
	Results.jacobian_full = numdiff(@(B) LL_boxcox(INPUT.Y,INPUT.Xa,INPUT.Xh,EstimOpt,B),f,Results.bhat_full,isequal(OptimOpt.FinDiffType,'central'),[ones(1,EstimOpt.NVARA),EstimOpt.TransActive,1]);
elseif EstimOpt.HessEstFix == 2
	Results.jacobian_full = jacobianest(@(B) LL_boxcox(INPUT.Y,INPUT.Xa,INPUT.Xh,EstimOpt,B),Results.bhat_full);
elseif EstimOpt.HessEstFix == 3
	Results.hess = hessian(@(B) sum(LL_boxcox(INPUT.Y,INPUT.Xa,INPUT.Xh,EstimOpt,B),1),Results.bhat_full);
end
if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
    Results.jacobian = Results.jacobian_full; ...
    Results.jacobian(:,[ones(1,EstimOpt.NVARA),EstimOpt.TransActive,1] == 0) = []; ...
    Results.hess = Results.jacobian'*Results.jacobian;
end

Results.ihess = inv(Results.hess);
Results.std = sqrt(diag(Results.ihess));

std_tmp = zeros(sum(EstimOpt.TransActive0 == 1,2),1); ...
std_tmp(EstimOpt.TransActive(EstimOpt.TransActive0 == 1) == 0) = NaN;
std_tmp(EstimOpt.TransActive(EstimOpt.TransActive0 == 1) == 1) = Results.std(EstimOpt.NVARA+1:EstimOpt.NVARA+EstimOpt.NVARBCt+EstimOpt.NVARBCh);
Results.std = [Results.std(1:EstimOpt.NVARA); std_tmp; Results.std(end)];

std_tmp = NaN(1 + EstimOpt.NVARA + EstimOpt.NVARH,1); ...
std_tmp(~isnan(EstimOpt.TransStruct)) = Results.std(EstimOpt.TransStructX + EstimOpt.NVARA); ...
Results.std_full = [Results.std(1:EstimOpt.NVARA); std_tmp; Results.std(end)];

Results.DetailsA = [Results.bhat(1:EstimOpt.NVARA) , Results.std(1:EstimOpt.NVARA) , pv(Results.bhat(1:EstimOpt.NVARA) , Results.std(1:EstimOpt.NVARA))];

Results.DetailsTS = [Results.bhat_full(EstimOpt.NVARA+1:EstimOpt.NVARA*2+1) , Results.std_full(EstimOpt.NVARA+1:EstimOpt.NVARA*2+1) , pv(Results.bhat_full(EstimOpt.NVARA+1:EstimOpt.NVARA*2+1) , Results.std_full(EstimOpt.NVARA+1:EstimOpt.NVARA*2+1))];
if any(EstimOpt.TransStruct(1:EstimOpt.NVARA+1) < 0)
    bx = bt_tmp(1:EstimOpt.NVARA+1);
    stdx = std_tmp(1:EstimOpt.NVARA+1);
    Results.DetailsTS(EstimOpt.TransStruct(1:EstimOpt.NVARA+1) < 0,:) = [exp(bx(EstimOpt.TransStruct(1:EstimOpt.NVARA+1) < 0)), exp(bx(EstimOpt.TransStruct(1:EstimOpt.NVARA+1) < 0)).*stdx(EstimOpt.TransStruct(1:EstimOpt.NVARA+1) < 0), pv(exp(bx(EstimOpt.TransStruct(1:EstimOpt.NVARA+1) < 0)), exp(bx(EstimOpt.TransStruct(1:EstimOpt.NVARA+1) < 0)).*stdx(EstimOpt.TransStruct(1:EstimOpt.NVARA+1) < 0))];
end
if EstimOpt.TransStruct(1) < 0 
    indx = find(EstimOpt.TransStruct(2:EstimOpt.NVARA+1) == EstimOpt.TransStruct(1));
    Results.DetailsTS([1; 1+indx],:) = [exp(exp(bt_tmp([1; 1+indx]))), exp(bt_tmp([1; 1+indx]) + exp(bt_tmp([1; 1+indx]))).*std_tmp([1; 1+indx]), pv(exp(exp(bt_tmp([1; 1+indx]))), exp(bt_tmp([1; 1+indx]) + exp(bt_tmp([1; 1+indx]))).*std_tmp([1; 1+indx]))];
end
Results.DetailsSIG = [exp(Results.bhat(end)) , exp(Results.bhat(end))*Results.std(end) , pv(exp(Results.bhat(end)) , exp(Results.bhat(end))*Results.std(end))];

if EstimOpt.NVARH > 0
	bx = bt_tmp(EstimOpt.NVARA+2:end);
	stdx = std_tmp(EstimOpt.NVARA+2:end);
	Results.DetailsHet = [Results.bhat_full(EstimOpt.NVARA*2+2:end-1) , Results.std_full(EstimOpt.NVARA*2+2:end-1) , pv(Results.bhat_full(EstimOpt.NVARA*2+2:end-1) , Results.std_full(EstimOpt.NVARA*2+2:end-1))];
	if any(EstimOpt.TransStruct(EstimOpt.NVARA+2:end) < 0)
        Results.DetailsHet(EstimOpt.TransStruct(EstimOpt.NVARA+2:end) < 0,:) = [exp(bx(EstimOpt.TransStruct(EstimOpt.NVARA+2:end) < 0)), exp(bx(EstimOpt.TransStruct(EstimOpt.NVARA+2:end) < 0)).*stdx(EstimOpt.TransStruct(EstimOpt.NVARA+2:end) < 0), pv(exp(bx(EstimOpt.TransStruct(EstimOpt.NVARA+2:end) < 0)), exp(bx(EstimOpt.TransStruct(EstimOpt.NVARA+2:end) < 0)).*stdx(EstimOpt.TransStruct(EstimOpt.NVARA+2:end) < 0))];
	end
else
    Results.DetailsHet = [];
end

Results.R = [Results.DetailsA; Results.DetailsTS; Results.DetailsHet; Results.DetailsSIG];

EstimOpt.params = length(Results.bhat) - length(EstimOpt.NotActive);
Results.stats = [Results.LL0; Results.LL; 1-Results.LL/Results.LL0; ((2*EstimOpt.params-2*Results.LL) + 2*EstimOpt.params*(EstimOpt.params+1)/(EstimOpt.NP-EstimOpt.params-1))/EstimOpt.NP; EstimOpt.NP; EstimOpt.params];
if EstimOpt.NVARH > 0
    Results.R_out = cell(4+EstimOpt.NVARA+3+2+EstimOpt.NVARH,11);
else
    Results.R_out = cell(4+EstimOpt.NVARA+3,11);
end
Results.R_out(1,1) = {'Box-Cox-exponential regression model'};
Results.R_out(2,2) = {'Linear coefficients'};
Results.R_out(2,6) = {'Transformation coefficients'};
Results.R_out(3,1:11) = [{'var'},{'coef.'},{'sig.'},{'st.err.'},{'p-value'},{'param.'},{'active'},{'coef.'},{'sig.'},{'st.err.'},{'p-value'}];
Results.R_out(4,1) = EstimOpt.NamesY;
Results.R_out(5:EstimOpt.NVARA+4,1) = EstimOpt.NamesA;
Results.R_out(5:EstimOpt.NVARA+4,2) = num2cell(Results.DetailsA(:,1));
Results.R_out(5:EstimOpt.NVARA+4,3) = cellstr(star_sig(Results.DetailsA(:,3)));
Results.R_out(5:EstimOpt.NVARA+4,4:5) = num2cell(Results.DetailsA(:,2:3));
Results.R_out(4:EstimOpt.NVARA+4,6) = EstimOpt.NamesTS(1:EstimOpt.NVARA+1);
Results.R_out(4:EstimOpt.NVARA+4,7) = num2cell(EstimOpt.TransActive(1:EstimOpt.NVARA+1)');
Results.R_out(4:EstimOpt.NVARA+4,8) = num2cell(Results.DetailsTS(:,1));
Results.R_out(4:EstimOpt.NVARA+4,9) = star_sig_cell(Results.DetailsTS(:,3));
Results.R_out(4:EstimOpt.NVARA+4,10:11) = num2cell(Results.DetailsTS(:,2:3));
Results.R_out(EstimOpt.NVARA+5,2) = {'Standard deviation of the error parameter'};
Results.R_out(EstimOpt.NVARA+6,1:5) = [{'var'},{'coef.'},{'sig.'},{'st.err.'},{'p-value'}];
Results.R_out(EstimOpt.NVARA+7,1) = {'Sigma'};
Results.R_out(EstimOpt.NVARA+7,2) = num2cell(Results.DetailsSIG(:,1));
Results.R_out(EstimOpt.NVARA+7,3) = star_sig_cell(Results.DetailsSIG(:,3));
Results.R_out(EstimOpt.NVARA+7,4:5) = num2cell(Results.DetailsSIG(:,2:3));
if EstimOpt.NVARH > 0
    Results.R_out(EstimOpt.NVARA+8,6) = {'Heteroskedascity parameters'};
    Results.R_out(EstimOpt.NVARA+9,1:11) = [{'var'},{''},{''},{''},{''},{'param.'},{'active'},{'coef.'},{'sig.'},{'st.err.'},{'p-value'}];
    Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,1) = EstimOpt.NamesH;
    Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,6) = EstimOpt.NamesTS(1+EstimOpt.NVARA+1:end);
    Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,7) = num2cell(EstimOpt.TransActive(EstimOpt.NVARA+2:end)');
%     num2cell((EstimOpt.TransActive(EstimOpt.NVARA+2:end)').*((EstimOpt.TransStruct(EstimOpt.NVARA+2:end)>0)'-(EstimOpt.TransStruct(EstimOpt.NVARA+2:end)<0)'));
    Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,9) = star_sig_cell(Results.DetailsHet(:,3));
    Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,[8,10,11]) = num2cell(Results.DetailsHet);
end

s = EstimOpt.spacing;
prec = EstimOpt.precision;
c1 = [Results.R_out(4:EstimOpt.NVARA+4,1);Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,1)];
s1 = max(cellfun(@(x) numel(x), c1));
c2 = [Results.R_out(5:EstimOpt.NVARA+4,2);Results.R_out(EstimOpt.NVARA+7,2);Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,2)];
c2 = c2(cellfun(@(x) ~isempty(x),c2));
c2 = c2(cellfun(@(x) isfinite(x),c2)); 
s2 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c2))+1,1])+prec+1;
c4 = [Results.R_out(5:EstimOpt.NVARA+4,4);Results.R_out(EstimOpt.NVARA+7,4);Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,4)];
c4 = c4(cellfun(@(x) ~isempty(x),c4));
c4 = c4(cellfun(@(x) isfinite(x),c4)); 
s4 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c4))+1,1])+prec+1;
c6 = [Results.R_out(4:EstimOpt.NVARA+4,6);Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,6)];
s6 = max(max(cellfun(@(x) numel(x), c6 ),6+1));
s7 = max(7,s);
c8 = [Results.R_out(5:EstimOpt.NVARA+4,8);Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,8)];
c8 = c8(cellfun(@(x) ~isempty(x),c8));
c8 = c8(cellfun(@(x) isfinite(x),c8)); 
s8 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c8))+1,1])+prec+1;
c10 = [Results.R_out(5:EstimOpt.NVARA+4,10); Results.R_out(EstimOpt.NVARA+10:EstimOpt.NVARA+9+EstimOpt.NVARH,10)];
c10 = c10(cellfun(@(x) ~isempty(x),c10));
c10 = c10(cellfun(@(x) isfinite(x),c10)); 
s10 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c10 ))+1,1])+prec+1;

fprintf('\n%*s\n\n', (s1+s2+3+s4+prec+2+s*5 + s6+1+s8+3+s10+prec+2+s*4 + size(Results.R_out{1,1},2))/2,Results.R_out{1,1})
fprintf('%-*s%-s\n', s1+s2+3+s4+prec+2+s*5, Results.R_out{2,2}, Results.R_out{2,6})
fprintf('%-*s%*s%*s%*s%s%-*s%*s%*s%*s%*s\n', s1,Results.R_out{3,1}, s2+s,Results.R_out{3,2}, s4+3+s,Results.R_out{3,4}, prec+2+s,Results.R_out{3,5},blanks(2*s), s6,Results.R_out{3,6}, s7,Results.R_out{3,7}, s8+s,Results.R_out{3,8}, s10+3+s,Results.R_out{3,10}, s+prec+2,Results.R_out{3,11})
fprintf('%-*s%-*s%*.0f%*.*f%-3s%*.*f%*.*f\n', s1+s2+3+s4+prec+2+s*5,Results.R_out{4,1}, s6,Results.R_out{4,6}, s7,Results.R_out{4,7}, s8+s,prec,Results.R_out{4,8}, Results.R_out{4,9}, s10+s,prec,Results.R_out{4,10}, 2+prec+s,prec,Results.R_out{4,11})
for i = 1:EstimOpt.NVARA
    fprintf('%-*s%*.*f%-3s%*.*f%*.*f%s%-*s%*.0f%*.*f%-3s%*.*f%*.*f\n', s1,Results.R_out{4+i,1}, s2+s,prec, Results.R_out{4+i,2}, Results.R_out{4+i,3}, s4+s,prec,Results.R_out{4+i,4}, 2+prec+s,prec,Results.R_out{4+i,5}, blanks(2*s),s6,Results.R_out{4+i,6}, s7,Results.R_out{4+i,7}, s8+s,prec,Results.R_out{4+i,8}, Results.R_out{4+i,9}, s10+s,prec,Results.R_out{4+i,10}, 2+prec+s,prec,Results.R_out{4+i,11})
end
fprintf('\n%-s\n', Results.R_out{EstimOpt.NVARA+5,2})
fprintf('%-*s%*s%*s%*s\n', s1,Results.R_out{EstimOpt.NVARA+6,1}, s2+s,Results.R_out{EstimOpt.NVARA+6,2}, s4+3+s,Results.R_out{EstimOpt.NVARA+6,4}, prec+2+s,Results.R_out{EstimOpt.NVARA+6,5})
fprintf('%-*s%*.*f%-3s%*.*f%*.*f\n', s1,Results.R_out{EstimOpt.NVARA+7,1}, s2+s,prec, Results.R_out{EstimOpt.NVARA+7,2}, Results.R_out{EstimOpt.NVARA+7,3}, s4+s,prec,Results.R_out{EstimOpt.NVARA+7,4}, 2+prec+s,prec,Results.R_out{EstimOpt.NVARA+7,5})
if EstimOpt.NVARH > 0
    fprintf('\n%-*s%-s\n', s1+s2+3+s4+prec+2+s*5, Results.R_out{EstimOpt.NVARA+8,2}, Results.R_out{EstimOpt.NVARA+8,6})
    fprintf('%-*s%*s%*s%*s%s%-*s%*s%*s%*s%*s\n', s1,Results.R_out{EstimOpt.NVARA+9,1}, s2+s,Results.R_out{EstimOpt.NVARA+9,2}, s4+3+s,Results.R_out{EstimOpt.NVARA+9,4}, prec+2+s,Results.R_out{EstimOpt.NVARA+9,5},blanks(2*s), s6,Results.R_out{EstimOpt.NVARA+9,6}, s7,Results.R_out{EstimOpt.NVARA+9,7}, s8+s,Results.R_out{EstimOpt.NVARA+9,8}, s10+3+s,Results.R_out{EstimOpt.NVARA+9,10}, s+prec+2,Results.R_out{EstimOpt.NVARA+9,11})
    for i = 1:EstimOpt.NVARH
        fprintf('%-*s%-*s%*.0f%*.*f%-3s%*.*f%*.*f\n', s1+s2+3+s4+prec+2+s*5,Results.R_out{EstimOpt.NVARA+9+i,1},s6,Results.R_out{EstimOpt.NVARA+9+i,6}, s7,Results.R_out{EstimOpt.NVARA+9+i,7}, s8+s,prec,Results.R_out{EstimOpt.NVARA+9+i,8}, Results.R_out{EstimOpt.NVARA+9+i,9}, s10+s,prec,Results.R_out{EstimOpt.NVARA+9+i,10}, 2+prec+s,prec,Results.R_out{EstimOpt.NVARA+9+i,11})
    end
end

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;


clocknote = clock;
tocnote = toc;
[~,DayName] = weekday(now,'long');

disp(' ');
disp(['LL at convergence: ',num2str(Results.LL,'%8.4f')])
disp(' ')
disp(['Estimation completed on ' DayName ', ' num2str(clocknote(1)) '-' sprintf('%02.0f',clocknote(2)) '-' sprintf('%02.0f',clocknote(3)) ' at ' sprintf('%02.0f',clocknote(4)) ':' sprintf('%02.0f',clocknote(5)) ':' sprintf('%02.0f',clocknote(6))])
disp(['Estimation took ' num2str(tocnote) ' seconds ('  num2str(floor(tocnote/(60*60))) ' hours ' num2str(floor(rem(tocnote,60*60)/60)) ' minutes ' num2str(rem(tocnote,60)) ' seconds).']);
disp(' ');