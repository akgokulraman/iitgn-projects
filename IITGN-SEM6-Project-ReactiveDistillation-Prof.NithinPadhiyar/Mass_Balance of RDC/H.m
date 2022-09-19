function frac = H(T)

%% delta H data
HFDo = -410; % Methyl Acetate kJ/mol
CpD = 92.58; % J/molK @ 25C

HFCo = -285.83; % Water kJ/mol
t = T/1000;
A = -203.6060; B = 1523.290; C = -3196.413; D = 2474.455;
CpC = A + B*t + C*t^2 + D*t^3; % J/molK @25C

HFBo = -238.42; % Methanol kJ/mol
CpB = 79.5; % J/mol K

HFAo = -483.52; % Acetic Acid kJ/mol
CpA = 123.1; %J/molK

HFD = HFDo + CpD*(T-(25+273.15))/1000; % kJ/mol
HFC = HFCo + CpC*(T-(25+273.15))/1000; % kJ/mol
HFB = HFBo + CpB*(T-(25+273.15))/1000; % kJ/mol
HFA = HFAo + CpA*(T-(25+273.15))/1000; % kJ/mol

lamda = HFD +HFC -(HFB +HFA);
deltaHv = 	30.32+ 42.6+51.6+ 35.21; % kJ/mol
frac = lamda/deltaHv;