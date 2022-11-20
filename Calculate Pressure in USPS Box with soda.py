# Calculate Pressure in USPS Box with Soda

# Enter the dimensions of box
w = input("Enter Box Width")
l = input("Enter Box Length")
h = input("Enter Box Height")

# Enter Temperature
T = input("Enter Temperature")
# Convert Temp based on input
Tunit = input("Temperature Units (use one letter abb.)")
if Tunit == "F": 
    Tact = ((T-32)*5/9) + 273.15
elif Tunit == "C": 
    Tact = T + 273.15
else: Tunit = T

# Enter Equilibrium Constant
Equi = input("Enter gas/liquid equilibrium constant")

# Enter how much of box is filled with soda
Fill = input("How full is the box?")

# Enter predicted concentraion of CO2 in soda
Conc = input("Enter CO2 Concentration")

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