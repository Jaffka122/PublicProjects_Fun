# Calculate Pressure in USPS Box with Soda

from math import sqrt

# Enter the dimensions of box
w = float(input("Enter Box Width: "))
l = float(input("Enter Box Length: "))
h = float(input("Enter Box Height: "))

# Enter Temperature
T = float(input("Enter Temperature: "))
# Convert Temp based on input
Tunit = input("Temperature Units (use one letter abb.): ")
if Tunit == "F": 
    Tact = ((T-32)*5/9) + 273.15
elif Tunit == "C": 
    Tact = T + 273.15
else: Tunit = T

# Enter Equilibrium Constant
Equi = float(input("Enter gas/liquid equilibrium constant: "))

# Enter how much of box is filled with soda
Fill = float(input("How full is the box? Enter Fraction: "))

# Enter predicted concentraion of CO2 in soda
Conc = float(input("Enter CO2 Concentration: "))

#Caculate values based on inputs
V = w*l*h
VL = Fill*V
Vv = (1-Fill)*V
nco2L = VL*Conc
n = Equi*nco2L

# Constants
R = 8.314 #J/mol*K

# Use ideal gas law to estimate pressure
Ptot = n*R*T/Vv

print("Ideal Pressure is ",Ptot," atm")

# Now to calculate without assuming ideal gas

# Below uses Soave-Redlich-Kwong equation

# Define Constants. Values from Engineering Toolbox or Google

Pc = 72.8 #atm
Tc = 304.13 #K
omega = 0.224
ac = 0.42748*(R**2)*(Tc**2)/Pc
b = 0.08664*R*Tc/Pc
kappa = 0.480 + 1.574*omega - 0.176*omega**2
a = ac*(1 + kappa*(1-sqrt(T/Tc))**2)
rhom = n/Vv # Molar Density [=] moles/m^3

Ptotni = (R*T*rhom/(1- b*rhom)) - (a*rhom**2)/(1+b*rhom)

print("Pressure from SRK Method is ", Ptotni," atm")