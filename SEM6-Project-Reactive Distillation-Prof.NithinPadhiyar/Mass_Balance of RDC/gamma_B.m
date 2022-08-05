function value = gamma_B(x,T)


deltag = zeros(4,4);
alpha = zeros(4,4);
G = zeros(4,4);
deltag(1,1) = 0;
deltag(1,2) = 0.5767;
deltag(1,3) = -187.72;
deltag(1,4) = -635.89;
deltag(2,1) = -39.582;
deltag(3,1) = 981.10;
deltag(4,1) = 1218.87;

deltag(2,2) = 0;
deltag(4,2) = 566.15;
deltag(3,2) = 921.33;
deltag(2,4) = 456.94;
deltag(2,3) = -245.90;

deltag(3,3) = 0;
deltag(4,3) = 879.14;
deltag(3,4) = 1709.49;

deltag(4,4) = 0;
deltag = inv(deltag);
taw = deltag/(8.314*T);

alpha(4,2) = 1.0293;
alpha(2,4) = 1.0293;

alpha(4,3) = 0.3830;
alpha(3,4) = 0.3830;

alpha(4,1) = 0.3600;
alpha(1,4) = 0.3600;

alpha(2,3) = 0.2989;
alpha(3,2) = 0.2989;

alpha(2,1) = 0.3055;
alpha(1,2) = 0.3055;

alpha(3,1) = 0.2960;
alpha(1,3) = 0.2960;

for i = 1:4
    for j = 1:4
        G(i,j) = exp(-alpha(i,j)*taw(i,j));
    end
end

sum = 0;

for i = 1:4
    sum = sum + G(i,2)*x(i);
end

%% Formula
value1 = 0;
for j = 1:4
    value1 = value1 + taw(j,2)*G(j,2)*x(j)/sum;
end
temp = 0;
for k = 1:4
    temp = temp + G(k,2)*x(k);
end
value1 = value1/temp;
value2 = 0;
for j = 1:4
    sum2 = 0;
    sum3 = 0;
    for m = 1:4
        sum2 = sum2 + G(m,j)*x(m);
    end
    for n = 1:4
        sum3 = sum3 + x(n)*taw(n,j)*G(n,j);
    end
    
    value2 = value2 + (x(j)*G(2,j)/sum2)*(taw(2,j) - sum3/sum2);
end
value = exp(value1 + value2);
% value = vapPressure(2,T);
    
    
    
    
    
    
    
end