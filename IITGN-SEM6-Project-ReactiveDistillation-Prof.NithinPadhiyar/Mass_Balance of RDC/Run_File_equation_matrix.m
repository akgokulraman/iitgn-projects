clear all;
clc;

%% Initiation
interval = [0.001 1.8]; % Time in Hours
initial = ones(1,30)*0.2; % A B C
% for i = 4:4:40
% initial(i) = 1 - 3*0.21;
% end

% initial(13) = 0.4962;
% initial(14) = 0.4808;
% initial(15) = 0.023;
% initial(16) = 0.01;
% 
% initial(25) = 0.4962;
% initial(26) = 0.4808;
% initial(27) = 0.023;
% initial(28) = 0.01;
initial(20) = 0.8;
initial(19) = 0.05;
initial(21) = 0.05;

initial(10) = 0.8;
initial(11) = 0.05;
initial(12) = 0.05;

% initial(28) = 0.8;
% initial(29) = 0.1;
% initial(30) = 0.1;
%% Assumptions 

% P is constant
% Ideal nature of methyl acetate
% Methyl acetate forms vapour immediately after forming.

%% ODE SOlving
% args = [x of different constituents at different plates]
% yA + yB + yC + yD = 1
[time,XConc]=ode45(@(t,args) equation_matrix(t,args), interval, initial);
XConc = abs(XConc);
a = size(XConc);
count_rows = int16(a(1));
count_col = int16(a(2));
y = zeros(count_rows,count_col);
for r = 1:count_rows
    for c = 1:3:count_col
        y(r,c) = XConc(r,c)*gamma_A([XConc(r,c);XConc(r,c+1);XConc(r,c+2);1-(XConc(r,c)+XConc(r,c+1)+XConc(r,c+2))],T(time(r),(c-1)/3)+273.15)*vapPressure(1,T(time(r),(c-1)/3)+273.15);
    end
    for c = 2:3:count_col
        y(r,c) = XConc(r,c)*gamma_B([XConc(r,c-1);XConc(r,c);XConc(r,c+1);1-(XConc(r,c-1)+XConc(r,c)+XConc(r,c+1))],T(time(r),(c-2)/3)+273.15)*vapPressure(2,T(time(r),(c-2)/3)+273.15);
    end
    for c = 3:3:count_col
        y(r,c) = XConc(r,c)*gamma_C([XConc(r,c-2);XConc(r,c-1);XConc(r,c);1-(XConc(r,c-2)+XConc(r,c-1)+XConc(r,c))],T(time(r),(c-3)/3)+273.15)*vapPressure(3,T(time(r),(c-3)/3)+273.15);
    end
end

yD = zeros(count_rows,count_col/3);
xD = zeros(count_rows,count_col/3);
for r = 1:count_rows
    for c = 3:3:count_col
        yD(r,c/3) = 0.98-(y(r,c)+y(r,c-1)+y(r,c-2));
    end
end
for r = 1:count_rows
    for c = 3:3:count_col
        xD(r,c/3) = 0.98-(XConc(r,c)+XConc(r,c-1)+XConc(r,c-2));
    end
end
% for r = 1:count_rows
%     for c = 3:3:count_col
%         xD(r,c/3) = yD(r,c/3)/(gamma_D([XConc(r,c-2);XConc(r,c-1);XConc(r,c);xD(r,c/3)],T(time(r),(c-3)/3)+273.15) * vapPressure(4,T(time(r),(c-3)/3)+273.15));
%     end
% end
% for r = 1:count_rows
%     for c = 3:3:count_col
%         yD(r,c/3) = xD(r,c/3)*gamma_D([XConc(r,c-2);XConc(r,c-1);XConc(r,c);xD(r,c/3)],T(time(r),(c-3)/3)+273.15) * vapPressure(4,T(time(r),(c-3)/3)+273.15);
%     end
% end
xD = abs(xD);
 
%% Plot 1

figure('Name','Liquid Mol-Fraction of D vs time','NumberTitle','off')
plot(time,xD(:,1));
hold on
plot(time,xD(:,2));
hold on
plot(time,xD(:,3));
hold on
plot(time,xD(:,4));
hold on
plot(time,xD(:,5));
hold on
plot(time,xD(:,6));
hold on
plot(time,xD(:,7));
hold on
plot(time,xD(:,8));
hold on
plot(time,xD(:,9));
hold on
plot(time,xD(:,10));
hold on
xlabel('time(hrs)')
ylabel('Liquid MolFraction of Methyl Acetate') 
legend('D','8','7','6','5','4','3','2','1','B')

%% Plot 2
xD_Final = xD(count_rows,:);
XD_Final = transpose(xD_Final);
N = [1:1:10];
figure('Name','Liquid Mol-Fraction vs Number','NumberTitle','off')
plot(N,xD_Final);
xticklabels({'D','8','7','6','5','4','3','2','1','B'})
xlabel('Tray No')
ylabel('Liquid MolFraction of Methyl Acetate at the end') 

%% Plot 3
% xA = zeros(count_rows,10);
% xA(:,1) = XConc(:,1);
% xA(:,2) = XConc(:,4);
% xA(:,3) = XConc(:,7);
% xA(:,4) = XConc(:,10);
% xA(:,5) = XConc(:,13);
% xA(:,6) = XConc(:,16);
% xA(:,7) = XConc(:,19);
% xA(:,8) = XConc(:,22);
% xA(:,9) = XConc(:,25);
% xA(:,10) = XConc(:,28);
% figure('Name','Liquid Mol-Fraction vs Time','NumberTitle','off')
% plot(time,xA(:,1));
% hold on
% plot(time,xA(:,2));
% hold on
% plot(time,xA(:,3));
% hold on
% plot(time,xA(:,4));
% hold on
% plot(time,xA(:,5));
% hold on
% plot(time,xA(:,6));
% hold on
% plot(time,xA(:,7));
% hold on
% plot(time,xA(:,8));
% hold on
% plot(time,xA(:,9));
% hold on
% plot(time,xA(:,10));
% hold on
% xlabel('time(hrs)')
% ylabel('Liquid MolFraction of Acetic Acid') 
% legend('D','8','7','6','5','4','3','2','1','B')



% figure('Name','Liquid Mol-Fraction vs Time of Acid and Acetate','NumberTitle','off')
% plot(time,xA(:,4));
% hold on
% plot(time,xD(:,4));
% xlabel('time(hrs)')
% ylabel('Liquid MolFraction of Acetic Acid') 
% legend('Acid','Acetate')