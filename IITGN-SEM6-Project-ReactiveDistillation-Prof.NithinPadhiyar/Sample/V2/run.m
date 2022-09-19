clc;
clear all;


interval = [0 5];
initial = [2 2 0 0 1];%Mols initial

% Assumptions: 
% 1. Inert gas also present in the system and has constant concentration
% 2. Liquid Vapour Equilibrium by Raoult's Law
% 3. Only Methyl Acetate Evaporates
% 4. Total Pressure always remain constant
% 5. Constant Temperature
%args=[FA FB FC FD V]
[t,Vec] = ode45(@(t,args) odes(t,args), interval, initial);

figure('Name','Overall Composition Profile','NumberTitle','off')%Overall
plot(t,Vec(:,1));
hold on
plot(t,Vec(:,2));
hold on
plot(t,Vec(:,3));
hold on
plot(t,Vec(:,4))
xlabel('time(hrs)')
ylabel('Composition in Reactor(mols)') 
legend({'FA-Acid','FB-Methanol','FC-Water','FD-Acetate'},'Location','northeast')


figure('Name','Acid and Acetate Profile','NumberTitle','off')%Acid and Acetate
plot(t,Vec(:,1));
hold on
plot(t,Vec(:,4));
xlabel('time(hrs)') 
ylabel('Composition in Reactor(mols)') 
legend({'FA-Acid','FD-Acetate'},'Location','northeast')

figure('Name','Conversion Profile: Case 2','NumberTitle','off')%Conversion
plot(t,1-Vec(:,1)/2);
xlabel('time(hrs)') 
ylabel('Conversion') 

