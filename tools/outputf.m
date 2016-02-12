function stop = outputf(x,optimvalues,state)
global B_backup 
persistent LL_backup IterTime

stop = false;
if isequal(state,'init')
    disp('')
	fprintf('%6s %6s %8s %16s %17s %18s %17s %12s \n','Iter.','Eval.','dB','Step','f(x)','df(x)','Opt. Cond.','Iter. time');
elseif isequal(state,'iter')
    IterTocNote = toc(IterTime);
    if optimvalues.iteration == 0
        if isempty(optimvalues.stepsize)
            optimvalues.stepsize = 0;
        end
        fprintf('%4d %6d %15.10f %15.10f %19.10f %15.10f %15.10f %9.4f\n',optimvalues.iteration,optimvalues.funccount,0,optimvalues.stepsize,optimvalues.fval,0,optimvalues.firstorderopt,IterTocNote);
        B_backup = x;
        LL_backup = optimvalues.fval;
    else
        dB = max(abs(x - B_backup));
        fprintf('%4d %6d %15.10f %15.10f %19.10f %15.10f %15.10f %9.4f\n',optimvalues.iteration,optimvalues.funccount,dB,optimvalues.stepsize,optimvalues.fval,LL_backup - optimvalues.fval,optimvalues.firstorderopt,IterTocNote);
        B_backup = x;
        LL_backup = optimvalues.fval;
    end
end

IterTime = tic;


