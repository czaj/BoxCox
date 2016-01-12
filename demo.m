clear
clc

%% Options
n = 10000; % no of respondents
rng(10000001);

%EstimOpt.Algorithm = 'trust-region';
EstimOpt.GradObj = 'off';
EstimOpt.Hessian = 'off';
EstimOpt.EPS = 1.e-12; % overall precision level
EstimOpt.HessEstFix = 1;

%% Generating Xa matrix
mu = [1,2,1];
sig = [3 4 3.5];
INPUT.Xa = [ones(n,1), unifrnd(mu(ones(n,1),:),sig(ones(n,1),:) )];


%% Models without heteroskedascity
% Only Y transformed with Box - Cox transformation
theta = 0.4;
eps = normrnd(0, 1, n,1);
b = [4;1 ;0.5; 1];

INPUT.Y = (1+theta*(INPUT.Xa*b + eps)).^(1/theta);
EstimOpt.TransStruct = [1,NaN,NaN,NaN, NaN]; % the same numbers for the same parameters, positive for box-cox transformation, negative for exponentiation, 0 for none
EstimOpt.b0 = [b; theta; 0];
Results.BOXCOX1 = BOXCOX(INPUT, EstimOpt); 

%Only Y transformed with exponential transformation
alpha = 1.2;
INPUT.Y = log(INPUT.Xa*b + eps)/alpha;
EstimOpt.TransStruct = [-1,NaN,NaN,NaN, NaN]; % the same numbers for the same parameters, positive for box-cox transformation, negative for exponentiation, 0 for none
EstimOpt.b0 = [b; log(alpha);0];
Results.BOXCOX2 = BOXCOX(INPUT, EstimOpt); 

%Y and one other variable transformed with Box - Cox transformation
Lambda = -0.5;
tmp = INPUT.Xa;
tmp(:,2) = (tmp(:,2).^Lambda-1)/Lambda;
INPUT.Y = (1+theta*(tmp*b + eps)).^(1/theta);
EstimOpt.TransStruct = [1,NaN,2,NaN, NaN]; % the same numbers for the same parameters, positive for box-cox transformation, negative for exponentiation, 0 for none
EstimOpt.b0 = [b; theta; Lambda;0];
Results.BOXCOX3 = BOXCOX(INPUT, EstimOpt);

%Y and two other variables transformed with exponential transformation
tmp = INPUT.Xa;
tmp(:,2) = exp(tmp(:,2)*Lambda); % Model wyrzuci exp(Lambda)
INPUT.Y = log(tmp*b + eps)/alpha;
EstimOpt.TransStruct = [-1,NaN,-2,NaN, NaN]; % the same numbers for the same parameters, positive for box-cox transformation, negative for exponentiation, 0 for none
EstimOpt.b0 = [b; log(alpha); Lambda;0];
Results.BOXCOX4 = BOXCOX(INPUT, EstimOpt);

%% Models with heteroskedascity
% Y transformed with Box - Cox transformation+one variable in variance also
lambdax = 0.8;
INPUT.Xh = unifrnd(1,2,n,1);
tmp = (INPUT.Xh.^lambdax - 1)/lambdax;
eps = normrnd(0, 0.5*tmp, n,1);
INPUT.Y = (1+theta*(INPUT.Xa*b + eps)).^(1/theta);
EstimOpt.b0 = [b; theta; lambdax;0];
EstimOpt.TransStruct = [1,NaN,NaN,NaN, NaN,2]; % the same numbers for the same parameters, positive for box-cox transformation, negative for exponentiation, 0 for none
Results.BOXCOX5 = BOXCOX(INPUT, EstimOpt);

% Y transformed with exponential transformation+one variable in variance also
eps = normrnd(0, 0.5*exp(-0.5*INPUT.Xh), n,1);
INPUT.Y = log(INPUT.Xa*b + eps)/alpha;
EstimOpt.TransStruct = [-1,NaN,NaN,NaN, NaN,-2]; % the same numbers for the same parameters, positive for box-cox transformation, negative for exponentiation, 0 for none
EstimOpt.b0 = [b; log(alpha); -0.5;0];
Results.BOXCOX6 = BOXCOX(INPUT, EstimOpt);

