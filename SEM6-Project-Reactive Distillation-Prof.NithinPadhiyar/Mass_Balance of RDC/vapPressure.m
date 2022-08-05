function value = vapPressure(n,T)

%% value of m
if(n == 1) % A
    value = 10^(4.68206 - (1642.54/(T-39.764)))/1;
end
if(n == 2) % B
    value = 10^(5.15853 - (1581.341/(T-33.50)))/1;
end
if(n == 3) % C
    value = 10^(5.0768 - (1659.793/(T-45.854)))/1;
end
 
if(n == 4) % D
    value = 10^(4.20364 - (1164.426/(T-52.69)))/1;
end
