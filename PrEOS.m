%% Peng Robinson Equation of State Calculator
% Input in necessary variables for substance
% Units are T = K; P = bar; V = m^3/mol


function Output = PrEOS(Tc,Pc,om,Cp,T,P)
Tr = 298.15; %reference temp at room temp
Pr = 1; %reference pressure at 1 bar

% Calculating necessary variables for dep. functions
kap = 0.37464+1.54226*om-0.26992*om^2;
R = 0.00008314;
a = (1 + kap*(1 - sqrt(T/Tc)))^2; %alf in mathCAD
b = 0.0788*R*Tc/Pc;
ac = (0.45724*(R^2)*(Tc^2))/Pc;
alpha = (1 + kap*(1- sqrt(T/Tc)))^2;
A = ac*alpha;
AA = A*P/((R*T)^2); %substitute for CA
BB = (P*b)/(R*T); %substitute for CB 
Da = -0.45724*((R^2)*(Tc^2)/Pc)*kap*sqrt(a/(T*Tc));

% Solving for the roots of the cubic equation to get Z
Z = roots([1 -(1-BB) (AA - 3*BB^2 - 2*BB) -(AA*BB - BB^2 -BB^3)]);
% Get only the real roots from equation
for i = 1:3
   if isreal(Z(i)) == 1
       Z(i) = Z(i);
   else
       Z(i) = 0;
   end
Z = sort(Z,'ascend');

% Calculate fugacity of liquid and gas phase
fl = (Z(1)-1) - log(Z(1) - BB) - (AA/(2*sqrt(2)*BB))*log((Z(1) +(1+sqrt(2))*BB)/(Z(1) + (1-sqrt(2))*BB));
fv = (Z(3)-1) - log(Z(3) - BB) - (AA/(2*sqrt(2)*BB))*log((Z(3) +(1+sqrt(2))*BB)/(Z(3) + (1-sqrt(2))*BB));
% Fugacity values
fugl = P*exp(fl);
fugv = P*exp(fv);
% Fugacity coefficients phi
phil = fugl/P;
phiv = fugv/P;

% Calculate entropy of liquid and gas phase
DelSL = (R*log(Z(1) - BB) + (Da/(2*sqrt(2)))*log((Z(1) + (1+sqrt(2))*BB)/(Z(1) + (1-sqrt(2))*BB)))*10^5;
DelSV = (R*log(Z(3) - BB) + (Da/(2*sqrt(2)))*log((Z(3) + (1+sqrt(2))*BB)/(Z(3) + (1-sqrt(2))*BB)))*10^5;

%Calculate enthalpy of liquid and gas phase
DelHL = (R*T*(Z(1) -1) + ((T*Da-A)/(2*sqrt(2)*b))*log((Z(1) + (1 + sqrt(2))*BB)/(Z(1) + (1-sqrt(2))*BB)))*10^5;
DelHV = (R*T*(Z(3) -1) + ((T*Da-A)/(2*sqrt(2)*b))*log((Z(3) + (1 + sqrt(2))*BB)/(Z(3) + (1-sqrt(2))*BB)))*10^5;

%Calcualte ideal roperty as reference states

DelHIG = Cp(1)*(T-Tr) + Cp(2)*(T^2-Tr^2)/2 + Cp(3)*(T^3-Tr^3)/3 + Cp(4)*(T^4-Tr^4)/4;
DelSIG = Cp(1)*log(T/Tr) + Cp(2)*(T-Tr) + Cp(3)*(T^2-Tr^2)/2 + Cp(4)*(T^3-Tr^3)/3 - R*(10^5)*log(P/Pr);

%Actually calculating the values for enthalpy and entropy based off ideal
%gas

SL = DelSL + DelSIG;
SV = DelSIG + DelSV;
HL = (DelHIG + DelHL)/1000;
HV = (DelHIG + DelHV)/1000;
VL = Z(1)*R*T/P;
VV = Z(3)*R*T/P;

%Output as table of values with units
VarNames = {'Variable', 'Liquid', 'Vapor'}; %Column Headers
Var = {'Temperature K', 'Pressure (bar)', 'Volume (m^3/mol)', 'Compressibility', 'Fugacity (bar)', 'Phi', 'Enthalpy (kJ/mol)' , 'Entropy (J/mol*K)'}; %Row Headers
Liq = [T, P, VL, Z(1), fugl, phil, HL, SL];
Vap = [T, P, VV, Z(3), fugv, phiv, HV, SV];

Output = table(Var', Liq', Vap', 'VariableNames', VarNames);
end