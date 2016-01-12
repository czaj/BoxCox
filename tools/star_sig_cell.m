function star_array = star_sig_cell(P)

for i = 1:length(P)
    if P(i) <= 0.01 
        star_array(i,1) = {'***'};
    elseif P(i) <= 0.05 
        star_array(i,1) = {'** '};
	elseif P(i) <= 0.1 
        star_array(i,1) = {'*  '};
    else 
        star_array(i,1) = {'   '};
    end
end
