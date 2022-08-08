%args = [NA;NB;NC;ND;V]
function diffeqns = odes(t,args)

% A is Limiting Reagent
% T is TB(Methyl Aetate) = 57.1 C
% P is 1 atm
T = 273.15 + 57.1;
Keq = 2.32*exp(782.98/T);
k = 9.732*10^8*exp(-6287.7/T);

MW = 74.08*10^-3;
rho = 932;
alpha = MW/rho;


FI = 1;% Molar flow rate of Inert gas
l = 0.227632; %P(D)/Po = 0.227632/1
yD = l*args(4)/(args(1)+args(2)+args(3)+args(4));

diffeqns(1,1) = -(k/args(5))*(args(1)*args(2)-args(3)*args(4)/Keq);
diffeqns(2,1) = -(k/args(5))*(args(1)*args(2)-args(3)*args(4)/Keq);
diffeqns(3,1) = +(k/args(5))*(args(1)*args(2)-args(3)*args(4)/Keq);
diffeqns(4,1) = +(k/args(5))*(args(1)*args(2)-args(3)*args(4)/Keq) - FI*(yD/(1-yD));
diffeqns(5,1) = -alpha*args(4);
end
