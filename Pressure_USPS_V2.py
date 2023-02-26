# Calculate Pressure in USPS Box with Soda

from math import sqrt
from sympy import *

# Enter the dimensions of box
w = float(input("Enter Box Width: "))
l = float(input("Enter Box Length: "))
h = float(input("Enter Box Height: "))

# Enter Temperature
T = float(input("Enter Temperature: "))
# Convert Temp based on input
Tunit = input("Temperature Units (use one letter abb.): ")
if Tunit == "F" or Tunit == 'f': 
    Tact = ((T-32)*5/9) + 273.15
elif Tunit == "C" or Tunit == 'c': 
    Tact = T + 273.15
else: Tunit = T

# Enter Equilibrium Constant
EquiC = float(input("Enter CO2/liquid equilibrium constant: "))
EquiAir = float(input("Enter air/liquid equilibrium constant: "))

# Constants
R = 8.314 #J/mol*K

# Enter how much of box is filled with soda
Fill = float(input("How full is the box? Enter Fraction: "))

# Enter predicted concentraion of CO2 in soda
Conc = float(input("Enter CO2 Concentration: "))

#Caculate values based on inputs
V = w*l*h
VL = Fill*V
Vv = (1-Fill)*V
nco2L = VL*Conc
n = EquiC*nco2L
nair = Vv*1000/(R*T)

# Use ideal gas law to estimate pressure
Ptot = (n+nair)*R*T/Vv

print("Ideal Pressure is ",Ptot," atm")

# Now to calculate without assuming ideal gas or neglecting equilibrium

# Below uses Soave-Redlich-Kwong equation, mass balance and equlibrium
# Equations to determine mols of total gas (air+CO2)

# Define Constants. Values from Engineering Toolbox or Google
# For simplicity, CO2 and air are assumed to have similar critical properties
Pc = 72.8 #atm
Tc = 304.13 #K
omega = 0.224
ac = 0.42748*(R**2)*(Tc**2)/Pc
b = 0.08664*R*Tc/Pc
kappa = 0.480 + 1.574*omega - 0.176*omega**2
a = ac*(1 + kappa*(1-sqrt(T/Tc))**2)

H = 1/0.034 # [=] kg*bar/mol

# Initial Conditions
Patm = 101325 # Pascals - Assume air in box initially at atmospheric
Pco2i = float(input("Enter Initial CO2 Partial Pressure: ")) 
# Can assume no CO2 has left the soda at start
Ptotini = Patm + Pco2i

#Define SRK function
rhomair = symbols('rhomair', real=True)
SRK = [(R*T*rhomair/(1- b*rhomair)) - ((a*rhomair**2)/(1+b*rhomair)) - Ptotini]


# At start, no VLE - Gas phase is pure air, liquid is just water and CO2
# Use SRK law to determine air initial conditions
rhomair = nonlinsolve(SRK,rhomair) # Molar Density [=] moles/m^3
print(rhomair)

from operator import itemgetter

rhomair_use = max(rhomair)[0]
print(rhomair_use)
# Calculate moles of CO2 in soda.
# Source: https://teacherscollegesj.org/how-much-co2-is-in-a-soda/
# 0.001 to convert L to m^3 and 44 g/mol is MW of CO2
nCO2i = 7*VL*0.001/44
nairi = rhomair_use*Vv
Hair = 0.1 # Vol air to Vol water
nair_liqf = Hair*VL/29

rhom = (nairi - nair_liqf)/Vv
Pair_final = (R*T*rhom/(1- b*rhom)) - (a*rhom**2)/(1+b*rhom)
PCO2_final = H*VL


Ptotsrk = PCO2_final + Pair_final
print("Pressure from SRK Method is ", Ptotsrk," atm")
