import aerotbx

h=0.5   #outlet height [m]
p0_pa=2.0
gamma=1.4
M0=2.0
phi0=0.0 # Degrees
nu0=aerotbx.flowprandtlmeyer(gamma=gamma,M=M0)[1]   # Output: [m,nu,mu]

# Know from common knowledge
phi4=phi7=phi9=0.0

# First use 3 characteristcs solving the problem
# See drawn figure
# First calculate the Mach number in region 3 using the pressure ratio since the expension is isentropic
p_0_p_t=aerotbx.flowisentropic(gamma=gamma, M=M0)[2] # Output: [mach,T,P,rho,area]
p_3_p_t=p_0_p_t*(1/p0_pa) # Calculate ratio total pressure vs pressure in 3
p_a_p_t=p_3_p_t

# Using this total pressure ratio the mach number in section 3 can be determined
M3=aerotbx.flowisentropic(gamma=gamma, p=p_3_p_t)[0][0]  # Returns an array of mach numbers
nu3=aerotbx.flowprandtlmeyer(gamma=gamma,M=M3)[1][0]

# Use MOC to obtain mu3 (gamma + from 0 to 3)
phi3=nu3-nu0+phi0

# No go back through expansion fan
dphi_I=(phi3-phi0)/3.0  # Phi step over fan I

phi1=phi0+dphi_I
phi2=phi0+2*dphi_I

nu1=nu0-phi0+phi1
nu2=nu0-phi0+phi2

M1=aerotbx.flowprandtlmeyer(gamma=gamma,nu=nu1)[0][0]
M2=aerotbx.flowprandtlmeyer(gamma=gamma,nu=nu2)[0][0]

# Now the first expansion fan is fully known, go one with the interference region
# Theta 5 is known to be 0
dphi_II=phi4-phi1
# Use MOC (gamma- 1 to 4)
nu4=nu1+phi1-phi4
# Use two times the MOC (gamma-25   gamma+45)
nu5=(nu2+nu4+phi2-phi4)*0.5
phi5=(nu2+phi2-nu4+phi4)*0.5

# Solve for region 6 using region 3 and 5
nu6=0.5*(nu5-phi5+nu3+phi3)
phi6=0.5*(phi5-nu5+nu3+phi3)

# phi7 is already known to be zero, using gamma-57
nu7=nu5+phi5-phi7

# Region 8
nu8=0.5*(nu7-phi7+nu6+phi6)
phi8=0.5*(phi7-nu7+nu6+phi6)

nu9=nu8+phi8-phi9
M9=aerotbx.flowprandtlmeyer(nu=nu9)[0][0]

# the pressure in 10 is the ambiant pressure
p_10_p_t=p_a_p_t
M10=aerotbx.flowisentropic(gamma=gamma, p=p_10_p_t)[0][0]
nu10=aerotbx.flowprandtlmeyer(gamma=gamma,M=M10)[1][0]

phi10=nu10-nu6+phi6

# Region 11 again use gamma+811 and gamma-1011
nu11=0.5*(nu8-phi8+nu10+phi10)
phi11=0.5*(phi8-nu8+nu10+phi10)

# Region 13
p_13_p_t=p_a_p_t
M13=aerotbx.flowisentropic(gamma=gamma, p=p_13_p_t)[0][0]
nu13=aerotbx.flowprandtlmeyer(gamma=gamma,M=M13)[1][0]

phi13=nu13-nu11+phi11

# Region 12
nu12=0.5*(nu9-phi9+nu11+phi11)
phi12=0.5*(phi9-nu9+nu11+phi11)

# Region 14
nu14=0.5*(nu12-phi12+nu13+phi13)
phi14=0.5*(phi12-nu12+nu13+phi13)

# Region 15
p_15_p_t=p_a_p_t
M15=aerotbx.flowisentropic(gamma=gamma, p=p_15_p_t)[0][0]
nu15=aerotbx.flowprandtlmeyer(gamma=gamma,M=M15)[1][0]

phi15=nu15-nu14+phi14
print nu4









