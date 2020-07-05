import aerotbx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

gamma = 1.4
M_e = 2.0
P_eP_a = 2.0
phi_e = 0.0 # deg
n = 20 # number of characteristics in the centered expansion fan
h = 1.0 # So total tunnel height is two times this value

# Plotting
plot_pres_c = True
plot_regions = True
plot_chars = True
plot_s = True
plot_pres_s = True
import matplotlib # Set graph font size
matplotlib.rcParams.update({'font.size': 16})

# Initialise the flow variables
phi = np.zeros(n+1 + 3*int(round(n*(n+1)/2.)))
nu = np.zeros(n+1 + 3*int(round(n*(n+1)/2.)))
M = np.zeros(n+1 + 3*int(round(n*(n+1)/2.)))
mu = np.zeros(n+1 + 3*int(round(n*(n+1)/2.)))

# Compute PM angle and Mach angle in the exit
M[0],nu[0],mu[0] = aerotbx.flowprandtlmeyer(gamma=gamma,M=M_e)

# Compute the total pressure ratio in the exit
P_eP_t = aerotbx.flowisentropic(gamma=gamma,M=M_e)[2] # Outputs [mach, T, P, rho, area]

# Compute total pressure ratio in the ambient after the expansion
P_aP_t = P_eP_t/P_eP_a

# With the total pressure ratio, compute the Mach number behind the expansion fan (region 3)
M[n] = aerotbx.flowisentropic(gamma=gamma,p=P_aP_t)[0][0]

# Compute the flow angle in the region behind the expansion fan (region 3)
(_,nu[n],mu[n]) = aerotbx.flowprandtlmeyer(gamma=gamma,M=M[n])
phi[n] = nu[n] - nu[0] + phi[0]

# Compute the flow angles and Mach in the fan
delta_phi_en = (phi[n]-phi[0])/(float(n)) # Since we take three characteristics

# Expansion wave at exit ABC
phi[1:n] = np.arange(1,n)*delta_phi_en
nu[1:n] = phi[1:n] + nu[0]

# Region BCE
count = n+1
out_count = np.zeros(n) # Extract region number going out
for i in range(n,0,-1):
    for j in range(i):
        # if j == 0 we have a wall problem with phi = 0
        if j == 0:
            # phi[count] = 0.0
            nu[count] = nu[count-i] + phi[count-i]
        # else we solve a piston-like problem with the adjacent regions
        else:
            phi[count] = 0.5*(nu[count-i] + phi[count-i] + phi[count-1] - nu[count-1])
            nu[count] = 0.5 * (nu[count - i] + phi[count - i] - phi[count - 1] + nu[count - 1])
        count += 1
    out_count[i - 1] = count-1
out_count = np.flipud(out_count)

phi_outer = np.zeros(n+1) # For later determining the intersection points
phi_outer[0] = phi[n]
region_2_count = np.zeros(n) # Extract region number going out
# Region DFG
for i in range(n,0,-1):
    for j in range(i):
        # For first interaction fan, correct the counting
        if i == n:
            # if j == 0 we have a pressure ratio problem
            if j == 0:
                # phi[count] = 0.0
                M = aerotbx.flowisentropic(gamma=gamma,p=P_aP_t)[0][0]
                nu[count] = aerotbx.flowprandtlmeyer(gamma=gamma,M=M)[1][0]
                phi[count] = nu[count] - nu[out_count[j]] + phi[out_count[j]]
                phi_outer[1] = phi[count]
            # else we solve a piston-like problem with the adjacent regions
            else:
                phi[count] = 0.5*(-nu[out_count[j]] + phi[out_count[j]] + phi[count-1] + nu[count-1])
                nu[count] = 0.5 * (nu[out_count[j]] - phi[out_count[j]] + phi[count - 1] + nu[count - 1])
        # For the rest of the interaction fan
        else:
            if j == 0:
                # phi[count] = 0.0
                M = aerotbx.flowisentropic(gamma=gamma,p=P_aP_t)[0][0]
                nu[count] = aerotbx.flowprandtlmeyer(gamma=gamma,M=M)[1][0]
                phi[count] = nu[count] - nu[count-i] + phi[count-i]
                phi_outer[n-i+1] = phi[count]
            else:
                phi[count] = 0.5 * (nu[count - 1] + phi[count - 1] + phi[count - i] - nu[count - i])
                nu[count] = 0.5 * (nu[count - i] - phi[count - i] + phi[count - 1] + nu[count - 1])
        if j == i-1: # Extrect region
            region_2_count[i-1] = count
        count += 1

region_2_count = np.flipud(region_2_count)

# Region HIK
for i in range(n,0,-1):
    for j in range(i):
        if i == n:
            if j == 0:
                # phi[count] = 0
                nu[count] = nu[region_2_count[j]] + phi[region_2_count[j]]
            else:
                phi[count] = 0.5*(nu[region_2_count[j]] + phi[region_2_count[j]] + phi[count-1] - nu[count-1])
                nu[count] = 0.5 * (nu[region_2_count[j]] + phi[region_2_count[j]] - phi[count - 1] + nu[count - 1])
        else:
            # if j == 0 we have a wall problem with phi = 0
            if j == 0:
                # phi[count] = 0.0
                nu[count] = nu[count-i] + phi[count-i]
            # else we solve a piston-like problem with the adjacent regions
            else:
                phi[count] = 0.5*(nu[count-i] + phi[count-i] + phi[count-1] - nu[count-1])
                nu[count] = 0.5 * (nu[count - i] + phi[count - i] - phi[count - 1] + nu[count - 1])
        count += 1

# Convert angles to flow characteristics
(Mach,_,mu) = aerotbx.flowprandtlmeyer(gamma=gamma,nu=nu) # Determine M using Prandtl-Meyer Fucntion
(_,T,PP_t,rho,_) = aerotbx.flowisentropic(gamma=gamma,M=Mach)  # Isentropic relations

# Define pressure ratio of each region with respect to the exit
PP_e = PP_t / P_eP_t

phi_outer = np.tan(phi_outer*np.pi/180.) # Convert phi in deg to slopes

# Plotting of characteristics
# Characteristic pattern
char_min = phi - mu
char_plus = phi + mu

angles_negative = np.zeros((n+1)*n + (n)*(n-1)*0.5) # Number of negative slopes

# Negative slopes in ABC
count = 0
for i in range(n):
    angles_negative[i] = (char_min[i] + char_min[i+1]) * 0.5
    count += 1
# Negative slopes in BCE
for i in range(n-1,0,-1):
    for j in range(i):
        angles_negative[count] = (char_min[count+n-i] + char_min[count+n-i+1]) * 0.5
        count += 1
# Negative slopes in DFG
count_regions = count + n + 1 # Skipped n+1 regions before next region with a negative slope line
for i in range(n,0,-1):
    for j in range(i):
        if i == n:
            angles_negative[count] = (char_min[out_count[j]] + char_min[count_regions]) * 0.5
        else:
            angles_negative[count] = (char_min[count_regions] + char_min[count_regions-i]) * 0.5
        count += 1
        count_regions += 1
# Negative slopes in HIK
for i in range(n - 1, 0, -1):
    for j in range(i):
        angles_negative[count] = (char_min[count + 2*n - i] + char_min[count + 2*n - i + 1]) * 0.5
        count += 1

## Positive slopes
angles_positive = np.zeros((n+1)*n-n+(n+1)*n*0.5)

count_regions = n+1
count = 0
count_plus_char = np.zeros(n) # For numbering gamma pluses when defining lines
# Positive slopes in BCE
for i in range(n,0,-1):
    for j in range(i):
        angles_positive[count] = 0.5 * (char_plus[count_regions]+char_plus[count_regions-i])
        count += 1
        count_regions += 1
    count_plus_char[n-i] = count-1

# Positive slopes in DFG
for i in range(n-1,0,-1):
    for j in range(i):
        angles_positive[count] = 0.5 * (char_plus[count_regions + n - i - 1] + char_plus[count_regions + n - i])
        count += 1
        count_regions += 1

# Positive slopes in HIK
count_regions += n
count_slopos_3 = np.zeros(n) # Extract slopes pos going out of HIK for latter use
for i in range(n, 0, -1):
    for j in range(i):
        if i == n:
            angles_positive[count] = 0.5 * (char_plus[count_regions] + char_plus[region_2_count[j]])
            count_slopos_3[j] = count
        else:
            angles_positive[count] = 0.5 * (char_plus[count_regions] + char_plus[count_regions - i])
        count += 1
        count_regions += 1

# Turn angles into slopes
slopes_negative = np.tan(angles_negative*np.pi/180.)
slopes_positive = np.tan(angles_positive*np.pi/180.)

## Compute intersection points
# Number of points
x = np.zeros(n*(n+1)+2+n + n*(n-1)*0.5)
y = np.zeros(n*(n+1)+2+n + n*(n-1)*0.5)

# Setting initial values and initial numbering values
y[0] = h
count = 1
count_pos  = count - 2
count_neg = count - 1
out_count_p = np.zeros(n) # Extract last all last points of row
color = 'b'
shear_col = 'g--'
#  Region BCE
for i in range(n,0,-1):
    for j in range(i):
        if i == n: # The first row of points
            if j == 0: # Points at centre line
                x[count] = (y[count] - y[0]) / slopes_negative[count_neg] + x[0]
                if plot_chars == True: # Plot
                    plt.plot(x[:count+1],y[:count+1])
            else: # Points not on centre line
                x[count] = (y[0] - y[count-1] + slopes_positive[count_pos]*x[count-1] - slopes_negative[count_neg]*x[0])\
                           /(slopes_positive[count_pos]-slopes_negative[count_neg])
                y[count] = slopes_positive[count_pos]*(x[count]-x[count-1]) + y[count-1]
                if plot_chars == True:
                    plt.plot((x[0],x[count]),(y[0],y[count]),color)
                    plt.plot((x[count - 1], x[count]), (y[count - 1], y[count]), color)
        else: # Not the first row
            if j == 0:
                x[count] = (y[count] - y[count-i]) / slopes_negative[count_neg] + x[count-i]
                if plot_chars == True:
                    plt.plot((x[count-i], x[count]), (y[count-i], y[count]), color)
            else:
                x[count] = (y[count-i] - y[count - 1] + slopes_positive[count_pos] * x[count - 1] - slopes_negative[
                    count_neg] * x[count-i]) / (slopes_positive[count_pos] - slopes_negative[count_neg])
                y[count] = slopes_positive[count_pos] * (x[count] - x[count - 1]) + y[count - 1]
                if plot_chars == True:
                    plt.plot((x[count - i], x[count]), (y[count - i], y[count]), color)
                    plt.plot((x[count-1], x[count]), (y[count-1], y[count]), color)
        count += 1
        count_pos += 1
        count_neg += 1
    out_count_p[n-i] = count - 1

# Region DFG
# Extract point for latter use
count_finals = np.zeros(n)
count_out_2 = np.zeros(n+1)
count_out_3 = np.zeros(n)
count_neg -= 1
for i in range(n+1,1,-1): # Same loop build up as the previous one
    for j in range(i):
        if i == n+1:
            if j == 0:
                x[count] = (y[0] - y[out_count_p[j]] + slopes_positive[count_plus_char[j]]*x[out_count_p[j]]
                            - phi_outer[j]*x[0]) / (slopes_positive[count_plus_char[j]] - phi_outer[j])
                y[count] = slopes_positive[count_plus_char[j]]*(x[count] - x[out_count_p[j]]) + y[out_count_p[j]]
                if plot_chars == True:
                    plt.plot((x[out_count_p[j]], x[count]), (y[out_count_p[j]], y[count]), color)
                    plt.plot((x[0], x[count]), (y[0], y[count]), shear_col)
            elif j == n:
                x[count] = -y[count-1]/slopes_negative[count_neg] + x[count-1]
                count_out_3[i-2] = count
                if plot_chars == True:
                    plt.plot((x[count-1], x[count]), (y[count - 1], y[count]), color)
            else:
                x[count] = (y[out_count_p[j]] - y[count-1] + slopes_negative[count_neg]*x[count-1]
                            - slopes_positive[count_plus_char[j]]*x[out_count_p[j]]) / (slopes_negative[count_neg] -
                            slopes_positive[count_plus_char[j]])
                y[count] = slopes_negative[count_neg]*(x[count] - x[count-1]) + y[count-1]
                if plot_chars == True:
                    plt.plot((x[count - 1], x[count]), (y[count - 1], y[count]), color)
                    plt.plot((x[out_count_p[j]], x[count]), (y[out_count_p[j]], y[count]), color)
            count_out_2[j] = count
        else:
            if j == 0:
                x[count] = (y[count-i-1] - y[count-i] + slopes_positive[count - (2*n-i+2)]*x[count-i]
                            - phi_outer[n-i+1]*x[count-i-1])/(slopes_positive[count - (2*n-i+2)] - phi_outer[n-i+1])
                y[count] = slopes_positive[count - (2*n-i+2)]*(x[count] - x[count-i]) + y[count-i]
                if plot_chars == True:
                    plt.plot((x[count - i - 1], x[count]), (y[count - i - 1], y[count]), shear_col)
                    plt.plot((x[count - i], x[count]), (y[count - i], y[count]), color)
            elif j == i-1:
                count_out_3[i-2] = count
            else:
                x[count] = (y[count-1] - y[count-i] + slopes_positive[count - (2*n-i+2)]*x[count-i]
                            - slopes_negative[count - (3+n-i)]*x[count-1])/(slopes_positive[count - (2*n-i+2)]
                                                                            - slopes_negative[count-(3+n-i)])
                y[count] = y[count-i] + slopes_positive[count - (2*n-i+2)]*(x[count] - x[count-i])
                if plot_chars == True:
                    plt.plot((x[count - 1], x[count]), (y[count - 1], y[count]), color)
                    plt.plot((x[count - i], x[count]), (y[count - i], y[count]), color)
        if j == i - 2:
            count_finals[n-i+1] = count
        count += 1
        count_neg += 1
count_out_3 = np.flipud(count_out_3)

# Region HIK
limit = False
for i in range(n,0,-1):
    for j in range(i):
        if i == n:
            if j == 0: # This point is already determined in the previous loop
                continue
            else:
                x[count_out_3[j]] = (y[count_out_3[j]-1] - y[count_out_3[j-1]] + slopes_positive[count_slopos_3[j-1]]
                        *x[count_out_3[j-1]] - slopes_negative[count_out_3[j]-2-j]*x[count_out_3[j]-1])\
                        /(slopes_positive[count_slopos_3[j-1]]-slopes_negative[count_out_3[j]-2-j])
                y[count_out_3[j]] = slopes_positive[count_slopos_3[j-1]]*(x[count_out_3[j]]-x[count_out_3[j-1]]) \
                        + y[count_out_3[j-1]]
                if plot_chars == True:
                    plt.plot((x[count_out_3[j]-1],x[count_out_3[j]]),(y[count_out_3[j]-1],y[count_out_3[j]]),color)
                    plt.plot((x[count_out_3[j-1]], x[count_out_3[j]]), (y[count_out_3[j-1]], y[count_out_3[j]]), color)
        elif i == n-1:
            if j == 0:
                x[count] = (y[count] - y[count_out_3[j+1]]) / slopes_negative[count-n-1] + x[count_out_3[j+1]]
                if plot_chars == True:
                    plt.plot((x[count_out_3[j+1]], x[count]), (y[count_out_3[j+1]], y[count]), color)
            else:
                x[count] = (y[count_out_3[j+1]] - y[count - 1] + slopes_positive[count-n-2] * x[count - 1] - slopes_negative[
                    count-n-1] * x[count_out_3[j+1]]) / (slopes_positive[count-n-2] - slopes_negative[count-n-1])
                y[count] = slopes_positive[count-n-2] * (x[count] - x[count - 1]) + y[count - 1]
                # Have the characteristics intersected?
                if x[count] < x[count_out_3[j+1]] and not limit and n>10: # Yes? Set x max
                    x_max = x[count-1]
                    limit = True # Loop must not be stoped after point, however point must be safed
                if plot_chars == True:
                    plt.plot((x[count-1], x[count]), (y[count-1], y[count]), color)
                    plt.plot((x[count_out_3[j + 1]], x[count]), (y[count_out_3[j + 1]], y[count]), color)
        else:
            if j == 0:
                x[count] = (y[count] - y[count - i]) / slopes_negative[count - n - 1] + x[count - i]
                if plot_chars == True:
                    plt.plot((x[count-i], x[count]), (y[count-i], y[count]), color)
            else:
                x[count] = (y[count - i] - y[count - 1] + slopes_positive[count - n - 2] * x[count - 1] -
                            slopes_negative[
                                count - n - 1] * x[count-i]) / (
                           slopes_positive[count - n - 2] - slopes_negative[count - n - 1])
                y[count] = slopes_positive[count - n - 2] * (x[count] - x[count - 1]) + y[count - 1]
                if plot_chars == True:
                    plt.plot((x[count], x[count-1]), (y[count], y[count-1]), color)
                    plt.plot((x[count-i], x[count]), (y[count-i], y[count]), color)
        if i < n:
            count += 1
        count_pos += 1
        count_neg += 1

# Define last shear boundary point
x[count] = x[count-1]
y[count] = y[-2+n*(n+1)*0.5+(n+1)*(n+2)*0.5] + phi_outer[-1]*(x[count] - x[-2+n*(n+1)*0.5+(n+1)*(n+2)*0.5])
if plot_chars == True:
    plt.plot((x[-2+n*(n+1)*0.5+(n+1)*(n+2)*0.5], x[count]), (y[-2+n*(n+1)*0.5+(n+1)*(n+2)*0.5], y[count]), shear_col)
    if plot_regions == False and plot_pres_s == False and plot_pres_s == False:
        if n > 10:
            plt.xlim((0, x_max))
        plt.show()

# Plot Mach distribution in the regions
if plot_regions == True:
    origin = [0.,0.]
    weight = Mach # One can choise between pressure PP_e,rho,T etc.

    # Setting color map
    cmap = plt.get_cmap("jet")
    # Take input array of weights to color by (in our case Mach numbers) and extract a color based on the Mach number norm
    N = (weight-weight.min())/(weight.max()-weight.min())
    weight_color = cmap(N)

    # First triangle - region 0
    plt.fill((origin[0],x[0],x[1],origin[0]),(origin[1],y[0],y[1],origin[1]),color=weight_color[0])

    # Loop over ABC
    for i in range(1,n+1):
        if i == n:
            plt.fill((x[0], x[i], x[n*(n+1)*0.5+1], x[0]), (y[0], y[i], y[n*(n+1)*0.5+1], y[0]), color=weight_color[i])
        else:
            plt.fill((x[0], x[i], x[i+1], x[0]), (y[0], y[i], y[i+1], y[0]), color=weight_color[i])

    # Loop over BCE
    count_regions = n+1
    for i in range(n,0,-1):
        for j in range(i):
            if j == 0:
                if i == 1:
                    plt.fill(
                        (x[count_regions - n], x[count_regions], x[count_regions + 1], x[count_regions - n]),
                        (y[count_regions - n], y[count_regions], y[count_regions + 1], y[count_regions - n]),
                        color=weight_color[count_regions])
                else:
                    plt.fill((x[count_regions-n], x[count_regions-n+1], x[count_regions+i-n], x[count_regions-n]),
                         (y[count_regions - n], y[count_regions - n + 1], y[count_regions+i-n], y[count_regions - n]),
                          color=weight_color[count_regions])
            elif j == i-1:
                plt.fill((x[count_regions - n], x[count_out_2[n-j-1]], x[count_out_2[n-j]],
                          x[count_regions - (n-i+1)], x[count_regions - n]),
                         (y[count_regions - n], y[count_out_2[n-j-1]], y[count_out_2[n-j]],
                          y[count_regions - (n-i+1)], y[count_regions - n]),
                         color=weight_color[count_regions])
            else:
                plt.fill((x[count_regions - n], x[count_regions - n + 1], x[count_regions + i - n],
                          x[count_regions + i - n - 1], x[count_regions - n]),
                         (y[count_regions - n], y[count_regions - n + 1], y[count_regions + i - n],
                          y[count_regions + i - n - 1], y[count_regions - n]),
                         color=weight_color[count_regions])
            count_regions += 1

    # Loop over DFG
    for i in range(n,0,-1):
        for j in range(i):
            if j == 0:
                if i == 1:
                    plt.fill(
                        (x[count_regions - i], x[count_regions - i + 1], x[-1], x[count_regions - i]),
                        (y[count_regions - i], y[count_regions - i + 1], y[-1], y[count_regions - i]),
                        color=weight_color[count_regions])
                else:
                    plt.fill((x[count_regions - i], x[count_regions - i + 1], x[count_regions + 1], x[count_regions - i]),
                             (y[count_regions - i], y[count_regions - i + 1], y[count_regions + 1], y[count_regions - i]),
                             color=weight_color[count_regions])

            else:
                plt.fill((x[count_regions], x[count_regions+1], x[count_regions-i+1],
                          x[count_regions-i], x[count_regions]),
                         (y[count_regions], y[count_regions+1], y[count_regions-i+1],
                          y[count_regions-i], y[count_regions]),
                         color=weight_color[count_regions])

            count_regions+=1

    # Loop over HIK
    for i in range(n-1, 0, -1):
        for j in range(i):
            if j == 0:
                if i == 1:
                    plt.fill(
                        (x[count_regions - n], x[count_regions - n + 1], x[count_regions - n + 2], x[count_regions - n]),
                        (y[count_regions - n], y[count_regions - n + 1], y[count_regions - n + 2], y[count_regions - n]),
                        color=weight_color[count_regions])
                elif i == n-1:
                    plt.fill((x[count_out_3[j]], x[count_out_3[j+1]], x[count_regions + i + 1 - n],
                              x[count_out_3[j]]),
                             (y[count_out_3[j]], y[count_out_3[j+1]], y[count_regions + i + 1 - n],
                              y[count_out_3[j]]), color=weight_color[count_regions])
                else:
                    plt.fill((x[count_regions - n], x[count_regions - n + 1], x[count_regions + i + 1 - n],
                              x[count_regions - n]),
                             (y[count_regions - n], y[count_regions - n + 1], y[count_regions + i + 1 - n],
                              y[count_regions - n]), color=weight_color[count_regions])
            else:
                if i == n-1:
                    plt.fill((x[count_out_3[j]], x[count_out_3[j+1]], x[count_regions + 1 + i - n],
                              x[count_regions + i - n], x[count_out_3[j]]),
                             (y[count_out_3[j]], y[count_out_3[j+1]], y[count_regions + 1 + i - n],
                              y[count_regions + i - n], y[count_out_3[j]]),
                             color=weight_color[count_regions])
                else:
                    plt.fill((x[count_regions - n], x[count_regions - n + 1], x[count_regions + i + 1 - n],
                             x[count_regions + i - n], x[count_regions - n]),
                             (y[count_regions - n], y[count_regions - n + 1], y[count_regions + i + 1 - n],
                             y[count_regions + i - n], y[count_regions - n]),
                             color=weight_color[count_regions])
            count_regions += 1
        count_regions += 1
    # Build up graph
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(weight)
    cb = plt.colorbar(m)
    cb.set_label("Mach number [-]")
    if n >10:
        plt.xlim((0,x_max))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

# Streamlines at 0.0*H and 0.5*H
slopes_phi = np.tan(phi*np.pi/180.)

## 0.0*h (SO 0.50*H of the entire exit)
# Check which points we need
x_s_c_i = np.where(y<10E-14)
x_s_c = np.hstack((0.,x[x_s_c_i]))
y_s_c = np.zeros(2.*n+1)

# Pressure distribution on centreline
# Extract pressure in selected points
phi_where_2 = np.where(abs(phi[(n+1)*n+n+1:])<10E-14)
len_phi_later = len(phi_where_2)
phi_later = np.ones(len_phi_later)*(1+(n+1)*n+n)
phi_0_i = np.hstack(( np.where(phi[:(n+1)*n*0.5+n+1]<10E-14), phi_later + phi_where_2 ))[0]

P_center = np.zeros(len(phi_0_i))
for i in range(len(P_center)):
    P_center[i] = PP_e[phi_0_i[i]]

# Plot pressure distribution
if plot_pres_c == True:
    for i in range(len(x_s_c)-1):
        plt.plot((x_s_c[i],x_s_c[i+1]),(P_center[i],P_center[i]),'b-')
        if n > 10:
            plt.xlim((0,x_max))
        if i<len(x_s_c)-2:
            # print i
            plt.plot((x_s_c[i+1],x_s_c[i+1]),(P_center[i],P_center[i+1]),'b-')
    plt.show()

## 0.5*h (SO 0.75*H of the entire exit) - Assuming it does not go throught DFG
x_s = np.zeros(3*n+2)
y_s = np.zeros(3*n+2)

# Initialise
y_s[0] = 0.5*h
color_s = 'g'
count = 1
# Particle through ABC
for i in range(n):
    x_s[count] = (y[0] - y_s[count - 1] + slopes_phi[count-1] * x_s[count - 1] - slopes_negative[count-1] * x[0]) \
               / (slopes_phi[count-1] - slopes_negative[count-1])
    y_s[count] = slopes_phi[count-1] * (x_s[count] - x_s[count - 1]) + y_s[count - 1]
    if plot_pres_s == True:
        plt.plot((x_s[count-1],x_s[count]),(Mach[count-1],Mach[count-1]), color_s)
        plt.plot((x_s[count], x_s[count]), (Mach[count], Mach[count-1]), color_s)
    count +=1

# Particle through CDEF
for i in range(n):
    if i == 0:
        x_s[count] = (y_s[count-1] - y[out_count_p[i]] + slopes_positive[count_plus_char[i]] * x[out_count_p[i]]
                - slopes_phi[count-1] * x_s[count-1]) / (slopes_positive[count_plus_char[i]] - slopes_phi[count-1])
        y_s[count] = slopes_positive[count_plus_char[i]] * (x_s[count] - x[out_count_p[i]]) + y[out_count_p[i]]
    else:
        x_s[count] = (y_s[count - 1] - y[out_count_p[i]] + slopes_positive[count_plus_char[i]] * x[out_count_p[i]]
                      - slopes_phi[out_count[i]] * x_s[count - 1]) / (
                     slopes_positive[count_plus_char[i]] - slopes_phi[out_count[i]])
        y_s[count] = slopes_positive[count_plus_char[i]] * (x_s[count] - x[out_count_p[i]]) + y[out_count_p[i]]
    if plot_pres_s == True:
        plt.plot((x_s[count-1],x_s[count]),(Mach[out_count[i]-1],Mach[out_count[i]-1]), color_s)
        plt.plot((x_s[count], x_s[count]), (Mach[out_count[i]], Mach[out_count[i]-1]), color_s)
    count += 1

# Particle through FGHI
i = n+1 # Maximum n+1 points
while x_s[count] < x[-1] and i>1:
    x_s[count] = (y[count_finals[n-i+1]] - y_s[count - 1] + slopes_phi[count_finals[n-i+1]] * x_s[count - 1]
                  - slopes_negative[count_finals[n-i+1]+(i-n-2)] * x[count_finals[n-i+1]]) \
                 / (slopes_phi[count_finals[n-i+1]] - slopes_negative[count_finals[n-i+1]+(i-n-2)])
    y_s[count] = slopes_phi[count_finals[n-i+1]] * (x_s[count] - x_s[count-1]) + y_s[count-1]
    if plot_pres_s == True:
        plt.plot((x_s[count-1],x_s[count]),(Mach[count_finals[n-i+1]],Mach[count_finals[n-i+1]]), color_s)
        plt.plot((x_s[count], x_s[count]), (Mach[count_finals[n-i+1]], Mach[count_finals[n-i+1]-1]), color_s)
    count += 1
    i-=1

# Last point
x_s[-1] = x[-1]
y_s[-1] = y_s[-2] + phi_outer[-1]*(x_s[-1] - x_s[-2])
# Plot pressure distribution
if plot_pres_s == True:
    plt.plot((x_s[-1], x_s[-2]), (Mach[count_finals[-1]+1], Mach[count_finals[-1]+1]), color_s)
    if n > 10:
        plt.xlim((0,x_max))
    plt.show()

# Plot stream line
if plot_s == True:
    plt.plot(x_s,y_s)
    if n>10:
        plt.xlim((0,x_max))
    plt.show()