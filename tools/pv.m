function result = pv(m,s,x)
if nargin < 3
   x = 0; 
end
if x == 0
    result=(1-normcdf(abs(m)./s,0,1))*2;
elseif x == 1
    result=(1-normcdf(abs(m-1)./s,0,1))*2;
end


