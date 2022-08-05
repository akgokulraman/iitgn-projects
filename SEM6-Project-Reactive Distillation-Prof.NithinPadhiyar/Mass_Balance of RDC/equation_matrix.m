function odematrix = equation_matrix(t,args)


% A - Acid
% B - Methanol
% C - Water
% D - Methyl Acetate

%% Constants
F3 = 30;
F6 = 30;
M = ones(1,10);
DR = 9; 
BR = 51;
RR = 5;
% t = 2;
%% Differential Variables
% xDA = args(1);
% xDB = args(2);
% xDC = args(3);
% xDD = args(4);
% x8A = args(5);
% x8B = args(6);
% x8C = args(7);
% x8D = args(8);
% x7A = args(9);
% x7B = args(10);
% x7C = args(11);
% x7D = args(12);
% x6A = args(13);
% x6B = args(14);
% x6C = args(15);
% x6D = args(16);
% x5A = args(17);
% x5B = args(18);
% x5C = args(19);
% x5D = args(20);
% x4A = args(21);
% x4B = args(22);
% x4C = args(23);
% x4D = args(24);
% x3A = args(25);
% x3B = args(26);
% x3C = args(27);
% x3D = args(28);
% x2A = args(29);
% x2B = args(30);
% x2C = args(31);
% x2D = args(32);
% x1A = args(33);
% x1B = args(24);
% x1C = args(35);
% x1D = args(36);
% xBA = args(37);
% xBB = args(38);
% xBC = args(39);
% xBD = args(40);
xDA = args(1);
xDB = args(2);
xDC = args(3); 
x8A = args(4);
x8B = args(5);
x8C = args(6);
x7A = args(7);
x7B = args(8);
x7C = args(9);
x6A = args(10);
x6B = args(11);
x6C = args(12);
x5A = args(13);
x5B = args(14);
x5C = args(15);
x4A = args(16);
x4B = args(17);
x4C = args(18);
x3A = args(19);
x3B = args(20);
x3C = args(21);
x2A = args(22);
x2B = args(23);
x2C = args(24);
x1A = args(25);
x1B = args(26);
x1C = args(27);
xBA = args(28);
xBB = args(29);
xBC = args(30);
xDD = 1 - (xDA + xDB + xDC);
x8D = 1 - (x8A + x8B + x8C);
x7D = 1 - (x7A + x7B + x7C);
x6D = 1 - (x6A + x6B + x6C);
x5D = 1 - (x5A + x5B + x5C);
x4D = 1 - (x4A + x4B + x4C);
x3D = 1 - (x3A + x3B + x3C);
x2D = 1 - (x2A + x2B + x2C);
x1D = 1 - (x1A + x1B + x1C);
xBD = 1 - (xBA + xBB + xBC);

%% Finding gammas of Methyl Acetate



%% Reaction Variables

R8A = -1* M(1,2) *k(273.15+T(t,8)) *(x8A*x8B - x8C*x8D/Keq(273.15+T(t,8)));
R8B = -1* M(1,2) *k(273.15+T(t,8)) *(x8A*x8B - x8C*x8D/Keq(273.15+T(t,8)));
R8C = +1* M(1,2) *k(273.15+T(t,8)) *(x8A*x8B - x8C*x8D/Keq(273.15+T(t,8)));
% R8D = +1* M(1,2) *k(273.15+T(t,8)) *(x8A*x8B - x8C*x8D/Keq(273.15+T(t,8)));

R7A = -1* M(1,3) *k(273.15+T(t,7)) *(x7A*x7B - x7C*x7D/Keq(273.15+T(t,7)));
R7B = -1* M(1,3) *k(273.15+T(t,7)) *(x7A*x7B - x7C*x7D/Keq(273.15+T(t,7)));
R7C = +1* M(1,3) *k(273.15+T(t,7)) *(x7A*x7B - x7C*x7D/Keq(273.15+T(t,7)));
% R7D = +1* M(1,3) *k(273.15+T(t,7)) *(x7A*x7B - x7C*x7D/Keq(273.15+T(t,7)));

R6A = -1* M(1,4) *k(273.15+T(t,6)) *(x6A*x6B - x6C*x6D/Keq(273.15+T(t,6)));
R6B = -1* M(1,4) *k(273.15+T(t,6)) *(x6A*x6B - x6C*x6D/Keq(273.15+T(t,6)));
R6C = +1* M(1,4) *k(273.15+T(t,6)) *(x6A*x6B - x6C*x6D/Keq(273.15+T(t,6)));
% R6D = +1* M(1,4) *k(273.15+T(t,6)) *(x6A*x6B - x6C*x6D/Keq(273.15+T(t,6)));

R5A = -1* M(1,5) *k(273.15+T(t,5)) *(x5A*x5B - x5C*x5D/Keq(273.15+T(t,5)));
R5B = -1* M(1,5) *k(273.15+T(t,5)) *(x5A*x5B - x5C*x5D/Keq(273.15+T(t,5)));
R5C = +1* M(1,5) *k(273.15+T(t,5)) *(x5A*x5B - x5C*x5D/Keq(273.15+T(t,5)));
% R5D = +1* M(1,5) *k(273.15+T(t,5)) *(x5A*x5B - x5C*x5D/Keq(273.15+T(t,5)));

R4A = -1* M(1,6) *k(273.15+T(t,4)) *(x4A*x4B - x4C*x4D/Keq(273.15+T(t,4)));
R4B = -1* M(1,6) *k(273.15+T(t,4)) *(x4A*x4B - x4C*x4D/Keq(273.15+T(t,4)));
R4C = +1* M(1,6) *k(273.15+T(t,4)) *(x4A*x4B - x4C*x4D/Keq(273.15+T(t,4)));
% R4D = +1* M(1,6) *k(273.15+T(t,4)) *(x4A*x4B - x4C*x4D/Keq(273.15+T(t,4)));

R3A = -1* M(1,7) *k(273.15+T(t,3)) *(x3A*x3B - x3C*x3D/Keq(273.15+T(t,3)));
R3B = -1* M(1,7) *k(273.15+T(t,3)) *(x3A*x3B - x3C*x3D/Keq(273.15+T(t,3)));
R3C = +1* M(1,7) *k(273.15+T(t,3)) *(x3A*x3B - x3C*x3D/Keq(273.15+T(t,3)));
% R3D = +1* M(1,7) *k(273.15+T(t,3)) *(x3A*x3B - x3C*x3D/Keq(273.15+T(t,3)));

R2A = -1* M(1,8) *k(273.15+T(t,2)) *(x2A*x2B - x2C*x2D/Keq(273.15+T(t,2)));
R2B = -1* M(1,8) *k(273.15+T(t,2)) *(x2A*x2B - x2C*x2D/Keq(273.15+T(t,2)));
R2C = +1* M(1,8) *k(273.15+T(t,2)) *(x2A*x2B - x2C*x2D/Keq(273.15+T(t,2)));
% R2D = +1* M(1,8) *k(273.15+T(t,2)) *(x2A*x2B - x2C*x2D/Keq(273.15+T(t,2)));

R1A = -1* M(1,9) *k(273.15+T(t,1)) *(x1A*x1B - x1C*x1D/Keq(273.15+T(t,1)));
R1B = -1* M(1,9) *k(273.15+T(t,1)) *(x1A*x1B - x1C*x1D/Keq(273.15+T(t,1)));
R1C = +1* M(1,9) *k(273.15+T(t,1)) *(x1A*x1B - x1C*x1D/Keq(273.15+T(t,1)));
% R1D = +1* M(1,9) *k(273.15+T(t,1)) *(x1A*x1B - x1C*x1D/Keq(273.15+T(t,1)));

RDC = +1* M(1,1) *k(273.15+T(t,9)) *(xDA*xDB - xDC*xDD/Keq(273.15+T(t,9)));
RBC = +1* M(1,10) *k(273.15+T(t,0)) * (xBA*xBB - xBC*xBD/Keq(273.15+T(t,0)));

RC = [RDC;R8C;R7C;R6C;R5C;R4C;R3C;R2C;R1C;RBC];

%% V and L

L = zeros(1,9);
V = zeros(1,8);
Liq = DR*RR;

L(9) = Liq;
L(8) = Liq;
L(7) = Liq;
L(6) = Liq;
L(5) = L(6) + H(273.15+T(t,5))*RC(5);
L(4) = L(5) + H(273.15+T(t,4))*RC(4);
L(3) = L(4);
L(2) = L(3);
L(1) = L(2);
% for i = 8:-1:1
%     L(i) = L(i+1) + (H(273.15+T(t,i)) * RC(10-i)); %change lamda/deltaHv
% end
% 
% 
V0 = L(1) - BR; % Not Sure
V(1) = V0;
V(2) = V(1);
V(3) = V(2);
V(4) = V(3) - (H(273.15+T(t,4))) * RC(4);
V(5) = V(4) - (H(273.15+T(t,5))) * RC(5);
V(6) = V(5);
V(7) = V(6);
V(8) = V(7);

% V(1) = V0 - (H(273.15+T(t,1))) * RC(10-1); %change lamda/deltaHv
% for i = 2:8
%     V(i) = V(i-1) - (H(273.15+T(t,i))) * RC(10-i); %change lamda/deltaHv 
% end


%% Material Balance
%% D - Condesner(10)

% xDD = 1 - (xDA + xDB + xDC);
% x8D = 1 - (x8A + x8B + x8C);

y8A = vapPressure(1,273.15+T(t,8)) * gamma_A([x8A;x8B;x8C;x8D],273.15+T(t,8)) * x8A;
y8B = vapPressure(2,273.15+T(t,8)) * gamma_B([x8A;x8B;x8C;x8D],273.15+T(t,8)) * x8B;
y8C = vapPressure(3,273.15+T(t,8)) * gamma_C([x8A;x8B;x8C;x8D],273.15+T(t,8)) * x8C;
% y8D = vapPressure(4,273.15+T(t,8)) * gamma_D([x8A;x8B;x8C;x8D],273.15+T(t,8)) * x8D;
y8D = 1 - (y8A+y8B+y8C);


odematrix(1,1) = (V(8)*y8A - DR*(1+RR)*xDA)/M(1,1);
odematrix(2,1) = (V(8)*y8B - DR*(1+RR)*xDB)/M(1,1);
odematrix(3,1) = (V(8)*y8C - DR*(1+RR)*xDC)/M(1,1);
% odematrix(4,1) = (V(8)*y8D - DR*(1+RR)*xDD)/M(1,1);

%% 273.15+Tray 8 - Rectifying

% x7D = 1 - (x7A + x7B + x7C);

y7A = vapPressure(1,273.15+T(t,7)) * gamma_A([x7A;x7B;x7C;x7D],273.15+T(t,7)) * x7A;
y7B = vapPressure(2,273.15+T(t,7)) * gamma_B([x7A;x7B;x7C;x7D],273.15+T(t,7)) * x7B;
y7C = vapPressure(3,273.15+T(t,7)) * gamma_C([x7A;x7B;x7C;x7D],273.15+T(t,7)) * x7C;
% y7D = vapPressure(3,273.15+T(t,7)) * gamma_D([x7A;x7B;x7C;x7D],273.15+T(t,7)) * x7D;
y7D = 1 - (y7A+y7B+y7C);
odematrix(4,1) = (L(9)*xDA + V(7)*y7A - (L(8)*x8A + V(8)*y8A) )/M(1,2);
odematrix(5,1) = (L(9)*xDB + V(7)*y7B - (L(8)*x8B + V(8)*y8B) )/M(1,2);
odematrix(6,1) = (L(9)*xDC + V(7)*y7C - (L(8)*x8C + V(8)*y8C) )/M(1,2);

% odematrix(5,1) = (L(9)*xDA + V(7)*y7A - (L(8)*x8A + V(8)*y8A) )/M(1,2);
% odematrix(6,1) = (L(9)*xDB + V(7)*y7B - (L(8)*x8B + V(8)*y8B) )/M(1,2);
% odematrix(7,1) = (L(9)*xDC + V(7)*y7C - (L(8)*x8C + V(8)*y8C) )/M(1,2);
% odematrix(8,1) = (L(9)*xDD + V(7)*y7D - (L(8)*x8D + V(8)*y8D) )/M(1,2);

%% 273.15+Tray 7 - Rectifying

% x6D = 1 - (x6A + x6B + x6C);

y6A = vapPressure(1,273.15+T(t,6)) * gamma_A([x6A;x6B;x6C;x6D],273.15+T(t,6)) * x6A;
y6B = vapPressure(2,273.15+T(t,6)) * gamma_B([x6A;x6B;x6C;x6D],273.15+T(t,6)) * x6B;
y6C = vapPressure(3,273.15+T(t,6)) * gamma_C([x6A;x6B;x6C;x6D],273.15+T(t,6)) * x6C;
% y6D = vapPressure(4,273.15+T(t,6)) * gamma_D([x6A;x6B;x6C;x6D],273.15+T(t,6)) * x6D;

y6D = 1 - (y6A+y6B+y6C);
odematrix(7,1) = (L(8)*x8A + V(6)*y6A - (L(7)*x7A + V(7)*y7A) )/M(1,3);
odematrix(8,1) = (L(8)*x8B + V(6)*y6B - (L(7)*x7B + V(7)*y7B) )/M(1,3);
odematrix(9,1) = (L(8)*x8C + V(6)*y6C - (L(7)*x7C + V(7)*y7C) )/M(1,3);
% odematrix(9,1) = (L(8)*x8A + V(6)*y6A - (L(7)*x7A + V(7)*y7A) )/M(1,3);
% odematrix(10,1) = (L(8)*x8B + V(6)*y6B - (L(7)*x7B + V(7)*y7B) )/M(1,3);
% odematrix(11,1) = (L(8)*x8C + V(6)*y6C - (L(7)*x7C + V(7)*y7C) )/M(1,3);
% odematrix(12,1) = (L(8)*x8D + V(6)*y6D - (L(7)*x7D + V(7)*y7D) )/M(1,3);

%% 273.15+Tray 6 - Feed of Acetic Acid

z6A = 0.9;
z6B = 0.05;
z6C = 0.05;
z6D = 0.0001;

% x5D = 1 - (x5A + x5B + x5C);
% 
y5A = vapPressure(1,273.15+T(t,5)) * gamma_A([x5A;x5B;x5C;x5D],273.15+T(t,5)) * x5A;
y5B = vapPressure(2,273.15+T(t,5)) * gamma_B([x5A;x5B;x5C;x5D],273.15+T(t,5)) * x5B;
y5C = vapPressure(3,273.15+T(t,5)) * gamma_C([x5A;x5B;x5C;x5D],273.15+T(t,5)) * x5C;
% y5D = vapPressure(4,273.15+T(t,5)) * gamma_C([x5A;x5B;x5C;x5D],273.15+T(t,5)) * x5D;

y5D = 1 - (y5A+y5B+y5C);
odematrix(10,1) = (L(7)*x7A + V(5)*y5A - (L(6)*x6A + V(6)*y6A) + R6A + F6*z6A)/M(1,4);
odematrix(11,1) = (L(7)*x7B + V(5)*y5B - (L(6)*x6B + V(6)*y6B) + R6B + F6*z6B)/M(1,4);
odematrix(12,1) = (L(7)*x7C + V(5)*y5C - (L(6)*x6C + V(6)*y6C) + R6C + F6*z6C)/M(1,4);
% odematrix(13,1) = (L(7)*x7A + V(5)*y6A - (L(6)*x6A + V(6)*y6A) + R6A + F6*z6A)/M(1,4);
% odematrix(14,1) = (L(7)*x7B + V(5)*y6B - (L(6)*x6B + V(6)*y6B) + R6B + F6*z6B)/M(1,4);
% odematrix(15,1) = (L(7)*x7C + V(5)*y6C - (L(6)*x6C + V(6)*y6C) + R6C + F6*z6C)/M(1,4);
% odematrix(16,1) = (L(7)*x7D + V(5)*y6D - (L(6)*x6D + V(6)*y6D) + R6C + F6*z6D)/M(1,4);

%% 273.15+Tray 5 - Reaction

% x4D = 1 - (x4A + x4B + x4C);

y4A = vapPressure(1,273.15+T(t,4)) * gamma_A([x4A;x4B;x4C;x4D],273.15+T(t,4)) * x4A;
y4B = vapPressure(2,273.15+T(t,4)) * gamma_B([x4A;x4B;x4C;x4D],273.15+T(t,4)) * x4B;
y4C = vapPressure(3,273.15+T(t,4)) * gamma_C([x4A;x4B;x4C;x4D],273.15+T(t,4)) * x4C;
% y4D = vapPressure(4,273.15+T(t,4)) * gamma_C([x4A;x4B;x4C;x4D],273.15+T(t,4)) * x4D;

y4D = 1 - (y4A+y4B+y4C);
odematrix(13,1) = (L(6)*x6A + V(4)*y4A - (L(5)*x5A + V(5)*y5A) + R5A)/M(1,5);
odematrix(14,1) = (L(6)*x6B + V(4)*y4B - (L(5)*x5B + V(5)*y5B) + R5B)/M(1,5);
odematrix(15,1) = (L(6)*x6C + V(4)*y4C - (L(5)*x5C + V(5)*y5C) + R5C)/M(1,5);
% odematrix(17,1) = (L(6)*x6A + V(4)*y4A - (L(5)*x5A + V(5)*y5A) + R5A)/M(1,5);
% odematrix(18,1) = (L(6)*x6B + V(4)*y4B - (L(5)*x5B + V(5)*y5B) + R5B)/M(1,5);
% odematrix(19,1) = (L(6)*x6C + V(4)*y4C - (L(5)*x5C + V(5)*y5C) + R5C)/M(1,5);
% odematrix(20,1) = (L(6)*x6D + V(4)*y4D - (L(5)*x5D + V(5)*y5D) + R5D)/M(1,5);


%% 273.15+Tray 4 - Reaction
% 
% x3D = 1 - (x3A + x3B + x3C);

y3A = vapPressure(1,273.15+T(t,3)) * gamma_A([x3A;x3B;x3C;x3D],273.15+T(t,3)) * x3A;
y3B = vapPressure(2,273.15+T(t,3)) * gamma_B([x3A;x3B;x3C;x3D],273.15+T(t,3)) * x3B;
y3C = vapPressure(3,273.15+T(t,3)) * gamma_C([x3A;x3B;x3C;x3D],273.15+T(t,3)) * x3C;
% y3D = vapPressure(4,273.15+T(t,3)) * gamma_C([x3A;x3B;x3C;x3D],273.15+T(t,3)) * x3D;

y3D = 1 - (y3A+y3B+y3C);
odematrix(16,1) = (L(5)*x5A + V(3)*y3A - (L(4)*x4A + V(4)*y4A) + R4A)/M(1,6);
odematrix(17,1) = (L(5)*x5B + V(3)*y3B - (L(4)*x4B + V(4)*y4B) + R4B)/M(1,6);
odematrix(18,1) = (L(5)*x5C + V(3)*y3C - (L(4)*x4C + V(4)*y4C) + R4C)/M(1,6);
% odematrix(21,1) = (L(5)*x5A + V(3)*y3A - (L(4)*x4A + V(4)*y4A) + R4A)/M(1,6);
% odematrix(22,1) = (L(5)*x5B + V(3)*y3B - (L(4)*x4B + V(4)*y4B) + R4B)/M(1,6);
% odematrix(23,1) = (L(5)*x5C + V(3)*y3C - (L(4)*x4C + V(4)*y4C) + R4C)/M(1,6);
% odematrix(24,1) = (L(5)*x5D + V(3)*y3D - (L(4)*x4D + V(4)*y4D) + R4D)/M(1,6);

%% 273.15+Tray 3 - Feed of methanol

z3B = 0.9;
z3A = 0.05;
z3C = 0.05;
z3D = 0.0001;
% 
% x2D = 1 - (x2A + x2B + x2C);

y2A = vapPressure(1,273.15+T(t,2)) * gamma_A([x2A;x2B;x2C;x2D],273.15+T(t,2)) * x2A;
y2B = vapPressure(2,273.15+T(t,2)) * gamma_B([x2A;x2B;x2C;x2D],273.15+T(t,2)) * x2B;
y2C = vapPressure(3,273.15+T(t,2)) * gamma_C([x2A;x2B;x2C;x2D],273.15+T(t,2)) * x2C;
% y2D = vapPressure(4,273.15+T(t,2)) * gamma_C([x2A;x2B;x2C;x2D],273.15+T(t,2)) * x2D;

y2D = 1 - (y2A+y2B+y2C);
odematrix(19,1) = (L(4)*x4A + V(2)*y2A - (L(3)*x3A + V(3)*y3A) + R3A + F3*z3A)/M(1,7);
odematrix(20,1) = (L(4)*x4B + V(2)*y2B - (L(3)*x3B + V(3)*y3B) + R3B + F3*z3B)/M(1,7);
odematrix(21,1) = (L(4)*x4C + V(2)*y2C - (L(3)*x3C + V(3)*y3C) + R3C + F3*z3C)/M(1,7);
% odematrix(25,1) = (L(4)*x4A + V(2)*y2A - (L(3)*x3A + V(3)*y3A) + R3A + F3*z3A)/M(1,7);
% odematrix(26,1) = (L(4)*x4B + V(2)*y2B - (L(3)*x3B + V(3)*y3B) + R3B + F3*z3B)/M(1,7);
% odematrix(27,1) = (L(4)*x4C + V(2)*y2C - (L(3)*x3C + V(3)*y3C) + R3C + F3*z3C)/M(1,7);
% odematrix(28,1) = (L(4)*x4D + V(2)*y2D - (L(3)*x3D + V(3)*y3D) + R3C + F3*z3D)/M(1,7);

%% 273.15+Tray 2 - Strip

% x1D = 1 - (x1A + x1B + x1C);

y1A = vapPressure(1,273.15+T(t,1)) * gamma_A([x1A;x1B;x1C;x1D],273.15+T(t,1)) * x1A;
y1B = vapPressure(2,273.15+T(t,1)) * gamma_B([x1A;x1B;x1C;x1D],273.15+T(t,1)) * x1B;
y1C = vapPressure(3,273.15+T(t,1)) * gamma_C([x1A;x1B;x1C;x1D],273.15+T(t,1)) * x1C;
% y1D = vapPressure(4,273.15+T(t,1)) * gamma_C([x1A;x1B;x1C;x1D],273.15+T(t,1)) * x1D;

y1D = 1 - (y1A+y1B+y1C);
odematrix(22,1) = (L(3)*x3A + V(1)*y1A - (L(2)*x2A + V(2)*y2A) )/M(1,8);
odematrix(23,1) = (L(3)*x3B + V(1)*y1B - (L(2)*x2B + V(2)*y2B) )/M(1,8);
odematrix(24,1) = (L(3)*x3C + V(1)*y1C - (L(2)*x2C + V(2)*y2C) )/M(1,8);
% odematrix(29,1) = (L(3)*x3A + V(1)*y1A - (L(2)*x2A + V(2)*y2A) )/M(1,8);
% odematrix(30,1) = (L(3)*x3B + V(1)*y1B - (L(2)*x2B + V(2)*y2B) )/M(1,8);
% odematrix(31,1) = (L(3)*x3C + V(1)*y1C - (L(2)*x2C + V(2)*y2C) )/M(1,8);
% odematrix(32,1) = (L(3)*x3D + V(1)*y1D - (L(2)*x2D + V(2)*y2D) )/M(1,8);

%% 273.15+Tray 1 - Strip

% xBD = 1 - (xBA + xBB + xBC);

yBA = vapPressure(1,273.15+T(t,0)) * gamma_A([xBA;xBB;xBC;xBD],273.15+T(t,0)) * xBA;
yBB = vapPressure(2,273.15+T(t,0)) * gamma_B([xBA;xBB;xBC;xBD],273.15+T(t,0)) * xBB;
yBC = vapPressure(3,273.15+T(t,0)) * gamma_C([xBA;xBB;xBC;xBD],273.15+T(t,0)) * xBC;
% yBD = vapPressure(4,273.15+T(t,0)) * gamma_C([xBA;xBB;xBC;xBD],273.15+T(t,0)) * xBD;

yBD = 1 - (yBA+yBB+yBC);
odematrix(25,1) = (L(2)*x2A + V0*yBA - (L(1)*x1A + V(1)*y1A) )/M(1,9);
odematrix(26,1) = (L(2)*x2B + V0*yBB - (L(1)*x1B + V(1)*y1B) )/M(1,9);
odematrix(27,1) = (L(2)*x2C + V0*yBC - (L(1)*x1C + V(1)*y1C) )/M(1,9);
% odematrix(33,1) = (L(2)*x2A + V0*yBA - (L(1)*x1A + V(1)*y1A) )/M(1,9);
% odematrix(34,1) = (L(2)*x2B + V0*yBB - (L(1)*x1B + V(1)*y1B) )/M(1,9);
% odematrix(35,1) = (L(2)*x2C + V0*yBC - (L(1)*x1C + V(1)*y1C) )/M(1,9);
% odematrix(36,1) = (L(2)*x2D + V0*yBD - (L(1)*x1D + V(1)*y1D) )/M(1,9);

%% B - Boiler(9)
VB = V0;
odematrix(28,1) = (L(1)*x1A - BR*xBA - VB*yBA)/M(1,10);
odematrix(29,1) = (L(1)*x1B - BR*xBB - VB*yBB)/M(1,10);
odematrix(30,1) = (L(1)*x1C - BR*xBC - VB*yBC)/M(1,10);
% odematrix(37,1) = (L(1)*x1A - BR*xBA - VB*yBA)/M(1,10);
% odematrix(38,1) = (L(1)*x1B - BR*xBB - VB*yBB)/M(1,10);
% odematrix(39,1) = (L(1)*x1C - BR*xBC - VB*yBC)/M(1,10);
% odematrix(40,1) = (L(1)*x1D - BR*xBD - VB*yBD)/M(1,10);





