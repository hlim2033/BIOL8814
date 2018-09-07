function dxdt = proteinproductino(t,x,pars)
m = x(1);
p = x(2);

dmdt = pars.a - pars.b*m;
dpdt = pars.r*m - pars.b*p;

dxdt = [dmdt; dpdt];
end
