clc;
clear all;

interval = [0 5];
initial = [0];

% Assumptions: 
% 1. Immediate Evaporation
% 2. Concentration of Liquid Evaporate is zero always

%args=[X]

[t,Vec] = ode45(@(t,args) odes(t,args), interval, initial);

figure('Name','Conversion Profile:Case 1','NumberTitle','off')
plot(t,Vec);
xlabel('time(hrs)') 
ylabel('Conversion') 
