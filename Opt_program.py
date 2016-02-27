"""
This is and attempt to port over an existing solar optimization program into python.
eventually, it will be more complex than the originial program, and have a GUI.
The goal being to eventually have a way to draw a layout for fingers on a solar cell, and 
immediatly calculate the resistance.
"""

#Variable Unit Conversions
Ang = 1e-10 # Angstrom unit in meters
nm = 1e-9 #Nanometer unit in meters
um = 1e-6 #Micron unit in meters
mm = 1e-3 #Milimeter unit in meters
cm = 1e-2 #centimeter unit in meters



# Instantiate CONSTANTS AND STANDARD PARAMATERS FOR CALCULATIONS
rho_au = 2.44e-8 #Resistivity of gold in [ohm-m] @ 20 C
q = 1.602e-19 # Coulomg unit of charge
N_n = 2e17 / cm^3 #doping of n type base layer
N_p = 5e18 / cm^3 #doping of p type emitter


# Initialize cell thickness detail all should be prompt's laters

#Emitter layer thickness 
thick_emit = 500*nm 
"""type_emit = input('N or P type')"""
#Base layer thickness 
thick_base = 500*um
"""type_emit = input("N or P type')"""

#Metal thickness's
back_metal = 2000 * Ang 
front_metal = 3000*Ang

#Initializing the cell size parameters
L_cell = input('Enter the Cell side (mm): ')*mm
A_cell = L_cell^2

#All parameters for optimization are based on coverage area
cvrg = input('Input Desired coverage area (%): ')/100

#Area of the front metal contact
A_metal = A_cell*cvrg

# Margin is the distance the fingers must be from the edge, this should be a set input
Margin = 25* um 

#BUSBAR parameters

#minimum
W_bus_min = input('minimum busbar width (um): ')*um
#busbar widths....think about this more on the algorithim front
"""W_bus = W_bus_min + (10*um*(array*L_cell/mm))"""
#Length of busbar
L_bus = L_cell - (2*Margin)
#Area of the busbar in m^2
A_bus = L_bus * W_bus 

#end of BUSBAR Parameters

#FINGER parameters.... These should mostly be prompts

#Minimum finger width
W_fin_min = 10*um
#range of finger widths
"""W_fin = W_fin_min + 2*um*array"""
#Finger length max
L_fin = (L_bus - W_bus)/2
#Potential Finger Area
A_fin = A_metal - A_bus
#Number of finger's calculations
"""N_fin = use floor function to find out arrays number of finger"""

#End of Finger parameters

#Characteristics of the cell. Some of these (FF) will have to be prompts 

#Vmp= 0.14
#Voc= 0.20
#Isc = 0.082
#FF = 0.531

#End of Characteristics of the cell

#Resistance Calculatoins
R_fc = 1.66 
R_emit = (1/
