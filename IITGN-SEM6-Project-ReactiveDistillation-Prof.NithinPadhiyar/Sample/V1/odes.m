function odematrix = odes(t,args)

X = args(1);
T = 273.15 + 57.1;
k = 9.732*10^8*exp(-6287.7/T);
CAo = 2;
MW = 74.08*10^-3;
rho = 932;
e = MW*CAo/rho;

odematrix(1) = k*CAo*((1-X)^2/(1+e*X));

end

