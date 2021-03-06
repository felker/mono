function root = rtsafe(fn_handle,x1,x2,xacc,coef1,coef2,coef3,coef4)
    MAX_ITER = 400;
    %recall, fn_handle should return both the function value and derivative
    % ENSURE THE ROOT IS BRACKETED
    % ---------------------------
    [fl, df] = fn_handle(x1,coef1,coef2,coef3,coef4);
    [fh, df] = fn_handle(x2,coef1,coef2,coef3,coef4); %do we not care about derivative?
    if (fl == 0) root=x1; return; end
    if (fh == 0) root=x2; return; end
    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
        error('Root must be bracketed fl = %f fh = %f',fl,fh);
    end
    if fl < 0.0
        xl=x1; %low v high?
        xh=x2;
    else
        xh=x1;
        xl=x2;
    end
    rts = 0.5*(x1+x2);
    dxold = abs(x2-x1);
    dx = dxold;
    [f, df] = fn_handle(rts,coef1,coef2,coef3,coef4);
    
    for i=1:MAX_ITER
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (abs(2.0*f) > abs(dxold*df)))
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;
            if (xl == rts) root=rts; return; end
        else 
            dxold=dx;
            dx=f/df;
            temp=rts;
            rts = rts - dx;
            if (temp == rts) root=rts; return;  end
        end
        if (abs(dx) < xacc) root=rts; return;  end
        [f, df] = fn_handle(rts,coef1, coef2, coef3,coef4);
        if (f < 0.0)
            xl=rts;
        else
            xh=rts;
        end
    end
        error('Max iterations exceeded'); 
    end
        