import aerotbx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

p0_pa=2.0
gamma=1.4
M0=2.0
phi0=0.0 # Degrees
n=19 # Number of characteristics
h=1

#Plotting
plot_pres_c=False
plot_pres_H=True
plot_mach=False


phi=np.zeros(n+1+3*round(n*(n+1)/2.))
nu=np.zeros(n+1+3*round(n*(n+1)/2.))
M=np.zeros(n+1+3*round(n*(n+1)/2.))
mu=np.zeros(n+1+3*round(n*(n+1)/2.))



(M[0],nu[0],mu[0])=aerotbx.flowprandtlmeyer(gamma=gamma,M=M0)   # Output: [m,nu,mu]

# # First calculate the Mach number in region n using the pressure ratio since the expension is isentropic
p_0_p_t=aerotbx.flowisentropic(gamma=gamma, M=M0)[2] # Output: [mach,T,P,rho,area]
p_n_p_t=p_0_p_t*(1/p0_pa) # Calculate ratio total pressure vs pressure in 3
p_a_p_t=p_n_p_t

#
# # Using this total pressure ratio the mach number in section 3 can be determined
M[n]=aerotbx.flowisentropic(gamma=gamma, p=p_n_p_t)[0][0]  # Returns an array of mach numbers
(_,nu[n],mu[n])=aerotbx.flowprandtlmeyer(gamma=gamma,M=M[n])

phi[n]=nu[n]-nu[0]+phi[0]
dphi_A=(phi[n]-phi[0])/float(n)

phi[1:n]=np.arange(1,n)*dphi_A

nu[1:n]=phi[1:n]+nu[0]

phi_outer=np.zeros(n+1)
phi_outer[0]=phi[n]
index_result_A=np.zeros(n)
sub=n+1

#First region
for i in range(n,0,-1):
     for j in range(i):
         if j==0:
             nu[sub] = nu[sub - i] + phi[sub - i]
         else:
             nu[sub] =  0.5 * (nu[sub-i] + phi[sub-i] - phi[sub-1] + nu[sub-1])
             phi[sub] = 0.5 * (nu[sub-i] + phi[sub-i] + phi[sub-1]- nu[sub-1])

         sub=sub+1
     index_result_A[i-1] = sub-1

index_result_A=np.flipud(index_result_A)

region_II_count=np.zeros(n)
#Second region
for i in range(n,0,-1):
     for j in range(i):
         if i==n:
             if j==0:
                 M=aerotbx.flowisentropic(gamma=gamma,p=p_a_p_t)[0][0]
                 nu[sub] = aerotbx.flowprandtlmeyer(gamma=gamma,M=M)[1][0]
                 phi[sub] = nu[sub] - nu[index_result_A[j]] + phi[index_result_A[j]]
                 phi_outer[1]=phi[sub]
             else:
                 nu[sub] =  0.5 * (nu[index_result_A[j]] - phi[index_result_A[j]] + phi[sub-1] + nu[sub-1])
                 phi[sub] = 0.5 * (-nu[index_result_A[j]] + phi[index_result_A[j]] + phi[sub-1] + nu[sub-1])

         else:
             if j==0:
                 M = aerotbx.flowisentropic(gamma=gamma, p=p_a_p_t)[0][0]
                 nu[sub] = aerotbx.flowprandtlmeyer(gamma=gamma, M=M)[1][0]
                 phi[sub] = nu[sub] - nu[sub-i] + phi[sub-i]
                 phi_outer[n-i+1]=phi[sub]
             else:
                 nu[sub] = 0.5 * (nu[sub - i] - phi[sub - i] + phi[sub - 1] + nu[sub - 1])
                 phi[sub] = 0.5 * (-nu[sub - i] + phi[sub - i] + phi[sub - 1] + nu[sub - 1])
         if j==i-1:
             region_II_count[i-1]=sub

         sub=sub+1
phi_outer=np.tan(phi_outer*np.pi/180.)

region_II_count=np.flipud(region_II_count)

# Third region
for i in range(n,0,-1):
     for j in range(i):
        if i==n:
            if j==0:
                nu[sub] = nu[region_II_count[j]] + phi[region_II_count[j]]
            else:
                nu[sub] = 0.5 * (nu[region_II_count[j]] + phi[region_II_count[j]] - phi[sub - 1] + nu[sub - 1])
                phi[sub] = 0.5 * (nu[region_II_count[j]] + phi[region_II_count[j]] + phi[sub - 1] - nu[sub - 1])
        else:
             if j==0:
                 nu[sub] = nu[sub - i] + phi[sub - i]
             else:
                 nu[sub] =  0.5 * (nu[sub-i] + phi[sub-i] - phi[sub-1] + nu[sub-1])
                 phi[sub] = 0.5 * (nu[sub-i] + phi[sub-i] + phi[sub-1]- nu[sub-1])

        sub=sub+1

(mach,_,mu)=aerotbx.flowprandtlmeyer(gamma=gamma,nu=nu)
(_,T,P_Pt,rho,_)=aerotbx.flowisentropic(gamma=gamma,M=mach)
P_Pe=P_Pt/p_0_p_t                #Pressure ratio between pressure in each region and exit pressure

#Plotting characteristics

#Characteristic pattern

char_min=phi-mu
char_plus=phi+mu


angles_neg=np.zeros(2*round(n*(n+1)/2.)+1*round((n-1)*(n)/2.))
angles_pos=np.zeros(3*round(n*(n+1)/2.)-n)

# # Negative slopes of first expansion
sub=0
for i in range(n):
    angles_neg[i]=(char_min[i]+char_min[i+1])*0.5
    sub+=1
# Negative slopes interaction I
for i in range(n-1,0,-1):
    for j in range(i):
        angles_neg[sub]=(char_min[sub+n-i] + char_min[sub+1+n-i])*0.5
        sub+=1
# Negative slopes interaction II
sub_regions=sub+n+1  # Skipped three regions
for i in range(n, 0, -1):
    for j in range(i):
        if i==n:
            angles_neg[sub]=(char_min[sub_regions]+char_min[index_result_A[j]])*0.5
        else:
            angles_neg[sub]=(char_min[sub_regions]+char_min[sub_regions-i])*0.5
        sub+=1
        sub_regions+=1

# For region 3
for i in range(n-1,0,-1):
    for j in range(i):
        angles_neg[sub]=(char_min[sub+n-i+n] + char_min[sub+1+n-i+n])*0.5
        sub+=1

# Now for the gamma+es
index_gamma_plus_I_II=np.zeros(n)
sub_regions=n+1
sub=0
# For the first region
for i in range(n,0,-1):
    for j in range(i):
        angles_pos[sub]=(char_plus[sub_regions]+char_plus[sub_regions-i])*0.5
        sub_regions+=1
        sub+=1
    index_gamma_plus_I_II[n-i]=sub-1

# Second region
for i in range(n-1,0,-1):
    for j in range(i):
        angles_pos[sub]=(char_plus[sub_regions+n-i-1]+char_plus[sub_regions+n-i])*0.5
        sub_regions+=1
        sub+=1

# Interaction 3
sub_regions=sub_regions+n
index_slopes_pos_II=np.zeros(n)
for i in range(n,0,-1):
    for j in range(i):
        if i==n:
            angles_pos[sub] = (char_plus[sub_regions] + char_plus[region_II_count[j]]) * 0.5
            index_slopes_pos_II[j]=sub
        else:
            angles_pos[sub]=(char_plus[sub_regions]+char_plus[sub_regions-i])*0.5
        sub_regions+=1
        sub+=1

# Turn angles into slopes
slopes_neg=np.tan(angles_neg*np.pi/180.)
slopes_pos=np.tan(angles_pos*np.pi/180.)

x=np.zeros(2*round(n*(n+1)/2.)+2+n+round((n-1)*(n)/2.))
y=np.zeros(2*round(n*(n+1)/2.)+2+n+round((n-1)*(n)/2.))

x[0]=0
y[0]=h

#First Fan
color='b-'
index_outer_numbers=np.zeros(n) # To check which gamma+ goes directly to II
sub=1
sub_pos=sub-2
sub_neg=sub-1
for i in range(n,0,-1):
    for j in range(i):
        if i==n:
            if j==0:
                x[sub] = (y[sub] - y[0] )/ slopes_neg[sub_neg] + x[0]
                plt.plot(x[:sub+1],y[:sub+1])
            else:
                x[sub] = ( y[0] - y[sub-1] + slopes_pos[sub_pos]*x[sub-1] -
                slopes_neg[sub_neg]*x[0])/(slopes_pos[sub_pos] - slopes_neg[sub_neg])
                y[sub] = slopes_pos[sub_pos]*(x[sub]-x[sub-1])+y[sub-1]
                plt.plot((x[0],x[sub]),(y[0],y[sub]),color)
                plt.plot((x[sub-1], x[sub]), (y[sub-1], y[sub]), color)
        else:
            if j==0:
                x[sub] = (y[sub] - y[sub-i]) / slopes_neg[sub_neg] + x[sub-i]
                plt.plot((x[sub-i], x[sub]), (y[sub-i], y[sub]), color)
            else:
                x[sub] = (y[sub-i] - y[sub - 1] + slopes_pos[sub_pos] * x[sub - 1] -
                          slopes_neg[sub_neg] * x[sub-i]) / (slopes_pos[sub_pos] - slopes_neg[sub_neg])
                y[sub] = slopes_pos[sub_pos] * (x[sub] - x[sub - 1]) + y[sub - 1]
                plt.plot((x[sub-i], x[sub]), (y[sub-i], y[sub]), color)
                plt.plot((x[sub - 1], x[sub]), (y[sub - 1], y[sub]), color)
        sub += 1
        sub_pos += 1
        sub_neg += 1
    index_outer_numbers[n-i]=sub-1

# Second fan
shearcolor='g--'
index_outer_II=np.zeros(n+1)
# SUB = 7 for n=3 :):)
sub_neg -= 1 # Due to point deffinition
index_final_B=np.zeros(n)  # Define last point of fan II
index_number_II_III=np.zeros(n)
for i in range(n+1,1,-1):
    for j in range(i):
        if i==n+1:
            if j==0:
                x[sub]=(y[0]-y[index_outer_numbers[j]]+slopes_pos[index_gamma_plus_I_II[j]]*x[index_outer_numbers[j]]-phi_outer[j]*x[j])/(slopes_pos[index_gamma_plus_I_II[j]]-phi_outer[j])
                y[sub]=slopes_pos[index_gamma_plus_I_II[j]]*(x[sub]-x[index_outer_numbers[j]]) + y[index_outer_numbers[j]]
                plt.plot((x[0], x[sub]), (y[0], y[sub]),shearcolor )
                plt.plot((x[index_outer_numbers[j]], x[sub]), (y[index_outer_numbers[j]], y[sub]), color)
            elif j==n:
                x[sub]=-y[sub-1]/slopes_neg[sub_neg]+x[sub-1]
                plt.plot((x[sub-1], x[sub]), (y[sub-1], y[sub]), color)
                index_number_II_III[i-2]=sub
            else:
                x[sub]=(y[index_outer_numbers[j]]-y[sub-1]+slopes_neg[sub_neg]*x[sub-1]-slopes_pos[index_gamma_plus_I_II[j]]
                        *x[index_outer_numbers[j]])/(slopes_neg[sub_neg]-slopes_pos[index_gamma_plus_I_II[j]])
                y[sub]=slopes_neg[sub_neg]*(x[sub]-x[sub-1])+y[sub-1]
                plt.plot((x[sub-1], x[sub]), (y[sub-1], y[sub]), color)
                plt.plot((x[index_outer_numbers[j]], x[sub]), (y[index_outer_numbers[j]], y[sub]), color)
            index_outer_II[j]=sub
        else:
            if j==0:
                x[sub]=(y[sub-i-1]-y[sub-i]+slopes_pos[sub-(2*n-i+2)]*x[sub-i]-phi_outer[n-i+1]*x[sub-i-1])/(slopes_pos[sub-(2*n-i+2)]-phi_outer[n-i+1])
                y[sub]= slopes_pos[sub-(2*n-i+2)]*(x[sub]-x[sub-i])+y[sub-i]
                plt.plot((x[sub-i-1], x[sub]), (y[sub-i-1], y[sub]), shearcolor)
                plt.plot((x[sub-i], x[sub]), (y[sub-i], y[sub]), color)

            elif j==i-1:
                index_number_II_III[i-2]=sub
            else:
                x[sub]=(y[sub-1]-y[sub-i] + slopes_pos[sub-(2*n-i+2)]*x[sub-i]-slopes_neg[sub-(3+n-i)]*x[sub-1] )/(slopes_pos[sub-(2*n-i+2)]-slopes_neg[sub-(3+n-i)])
                y[sub]= y[sub-i]+slopes_pos[sub-(2*n-i+2)]*(x[sub]-x[sub-i])
                plt.plot((x[sub-1], x[sub]), (y[sub-1], y[sub]), color)
                plt.plot((x[sub-i], x[sub]), (y[sub-i], y[sub]), color)
        if j == i-2:
            index_final_B[n-i+1]=sub
        sub+=1
        sub_neg +=1

index_number_II_III=np.flipud(index_number_II_III)


limit=False

# Determine x_max
for i in range(n,0,-1):
    for j in range(i):
        if i==n:
            if j==0:
                continue
            else:
                x[index_number_II_III[j]] = (y[index_number_II_III[j]-1] - y[index_number_II_III[j-1]] + slopes_pos[index_slopes_pos_II[j-1]]*x[index_number_II_III[j-1]] -
                slopes_neg[index_number_II_III[j]-2-j]*x[index_number_II_III[j]-1])/(slopes_pos[index_slopes_pos_II[j-1]] - slopes_neg[index_number_II_III[j]-2-j])
                y[index_number_II_III[j]] = slopes_pos[index_slopes_pos_II[j-1]]*(x[index_number_II_III[j]]-x[index_number_II_III[j-1]])+y[index_number_II_III[j-1]]
                plt.plot((x[index_number_II_III[j]-1],x[index_number_II_III[j]]),(y[index_number_II_III[j]-1],y[index_number_II_III[j]]),color)
                plt.plot((x[index_number_II_III[j-1]], x[index_number_II_III[j]]), (y[index_number_II_III[j-1]], y[index_number_II_III[j]]), color)
        elif i==n-1:
            if j==0:
                x[sub] = (y[sub] - y[index_number_II_III[j+1]]) / slopes_neg[sub-n-1] + x[index_number_II_III[j+1]]
                plt.plot((x[index_number_II_III[j+1]], x[sub]), (y[index_number_II_III[j+1]], y[sub]), color)
            else:
                x[sub] = (y[index_number_II_III[j+1]] - y[sub - 1] + slopes_pos[sub-n-2] * x[sub - 1] -
                          slopes_neg[sub-n-1] * x[index_number_II_III[j+1]]) / (slopes_pos[sub-n-2] - slopes_neg[sub-n-1])
                y[sub] = slopes_pos[sub-n-2] * (x[sub] - x[sub - 1]) + y[sub - 1]
                plt.plot((x[sub], x[sub-1]), (y[sub], y[sub-1]), color)
                plt.plot((x[index_number_II_III[j + 1]], x[sub]), (y[index_number_II_III[j + 1]], y[sub]), color)
                if x[sub]<x[index_number_II_III[j+1]] and not limit:
                    x_max=x[sub-1]
                    limit = True

        else:
            if j == 0:
                x[sub] = (y[sub] - y[sub - i]) / slopes_neg[sub - n - 1] + x[sub - i]
                plt.plot((x[sub-i], x[sub]), (y[sub-i], y[sub]), color)
            else:
                x[sub] = (y[sub-i] - y[sub - 1] + slopes_pos[sub - n - 2] * x[sub - 1] -
                          slopes_neg[sub - n - 1] * x[sub-i]) / (
                          slopes_pos[sub - n - 2] - slopes_neg[sub - n - 1])
                y[sub] = slopes_pos[sub - n - 2] * (x[sub] - x[sub - 1]) + y[sub - 1]
                plt.plot((x[sub-i], x[sub]), (y[sub-i], y[sub]), color)
                plt.plot((x[sub - 1], x[sub]), (y[sub - 1], y[sub]), color)
        if i<n:
            sub += 1
            sub_pos += 1
            sub_neg += 1


# Define the last point where the shear layer
x[sub]=x[sub-1]  # Define last x coordinate
y[sub]=y[-2+n*(n+1)*0.5+(n+1)*(n+2)*0.5]+phi_outer[-1]*(x[sub]-x[-2+n*(n+1)*0.5+(n+1)*(n+2)*0.5])
plt.plot((x[-2+n*(n+1)*0.5+(n+1)*(n+2)*0.5], x[sub]), (y[-2+n*(n+1)*0.5+(n+1)*(n+2)*0.5], y[sub]), shearcolor)

plt.scatter(x,y)
if n>10:
    plt.xlim(0,x_max)



if plot_mach==True:
    origin = [0.,0.]
    weight=mach
    # Color plot of Mach number
    # First setting colorbar from Mmax to Mmin
    # Use fill to fill between points

    #plt.fill(x,y,colour) x and y of the first have to be specified twice
    cmap=plt.get_cmap("jet")
    N=(weight-weight.min())/(weight.max()-weight.min())
    mach_color = cmap(N)

    # First traingle
    plt.fill((x[0],origin[0],x[1],x[0]),(y[0],origin[1],y[1],y[0]),color=mach_color[0])

    #First expension fan
    for i in range(1,n+1):
        if i==n:
            plt.fill((x[0], x[i] , x[n*(n+1)*0.5+1] , x[0]), (y[0], y[i] , y[n*(n+1)*0.5+1] , y[0]), color=mach_color[i])

        else:
            plt.fill((x[0], x[i] , x[i+1] , x[0]), (y[0], y[i] , y[i+1] , y[0]), color=mach_color[i])

    # Region I
    region_count=n+1
    for i in range(n,0,-1):
        for j in range(i):
            if j==0:
                if i==1:
                    plt.fill((x[region_count - n], x[region_count], x[region_count+1], x[region_count - n]),
                             (y[region_count - n], y[region_count], y[region_count + 1], y[region_count - n]),
                             color=mach_color[region_count])

                else:
                    plt.fill((x[region_count-n], x[region_count-n+1], x[region_count+i-n], x[region_count-n]),
                         (y[region_count-n], y[region_count-n+1], y[region_count+i-n], y[region_count-n]), color=mach_color[region_count])
            elif j==i-1: # For this loop one needs the outer points of II
                plt.fill((
                         x[region_count - n], x[index_outer_II[n-j-1]], x[index_outer_II[n-j]], x[region_count-(n-i+1)],
                         x[region_count - n]),
                         (
                         y[region_count - n], y[index_outer_II[n-j-1]], y[index_outer_II[n-j]], y[region_count-(n-i+1)],
                         y[region_count - n]),
                         color=mach_color[region_count])
            else:
                plt.fill((x[region_count - n], x[region_count - n + 1], x[region_count + i - n],x[region_count + i - n-1], x[region_count - n]),
                         (y[region_count - n], y[region_count - n + 1], y[region_count + i - n],y[region_count + i - n-1], y[region_count - n]),
                         color=mach_color[region_count])

            region_count+=1

    # The same for region II
    for i in range(n,0,-1):
        for j in range(i):
            if j==0:
                if i==1:
                    plt.fill((x[region_count-1], x[region_count], x[-1], x[region_count - 1]),
                             (y[region_count - 1], y[region_count], y[-1], y[region_count - 1]),
                             color=mach_color[region_count])
                else:
                    plt.fill((x[region_count - i], x[region_count-i+1], x[region_count + 1], x[region_count - i]),
                         (y[region_count - i], y[region_count-i+1], y[region_count + 1], y[region_count - i]),
                         color=mach_color[region_count])
            else:
                plt.fill((
                    x[region_count], x[region_count+1], x[region_count-i+1],
                    x[region_count-i],
                    x[region_count]),
                    (
                        y[region_count], y[region_count+1], y[region_count-i+1],
                        y[region_count-i],
                        y[region_count]),
                    color=mach_color[region_count])

            region_count+=1
    # Third reflection

    for i in range(n-1, 0, -1):
        for j in range(i):
            if j == 0:
                if i == 1:
                    plt.fill((x[region_count - n], x[region_count-n+1], x[region_count -n + 2], x[region_count - n]),
                             (y[region_count - n], y[region_count-n+1], y[region_count -n + 2], y[region_count - n]),
                             color=mach_color[region_count])

                elif i==n-1:
                    plt.fill((x[index_number_II_III[j]], x[index_number_II_III[j+ 1]], x[region_count + i - n + 1],
                              x[index_number_II_III[j]]),
                             (y[index_number_II_III[j]], y[index_number_II_III[j + 1]], y[region_count + i - n + 1],
                              y[index_number_II_III[j]]), color=mach_color[region_count])
                else:
                    plt.fill((x[region_count - n], x[region_count - n + 1], x[region_count + i - n +1],
                              x[region_count - n]),
                             (y[region_count - n], y[region_count - n + 1], y[region_count + i - n+1],
                              y[region_count - n]), color=mach_color[region_count])

            else:
                if i==n-1:
                    plt.fill((x[index_number_II_III[j]], x[index_number_II_III[j+1]], x[region_count+i-n+1],
                          x[region_count + i - n], x[index_number_II_III[j]]),
                         (y[index_number_II_III[j]], y[index_number_II_III[j+1]], y[region_count+i-n+1],
                          y[region_count + i - n ], y[index_number_II_III[j]]),
                         color=mach_color[region_count])
                else:
                    plt.fill((x[region_count - n], x[region_count - n + 1], x[region_count + i - n +1],x[region_count + i - n],
                              x[region_count - n]),
                             (y[region_count - n], y[region_count - n + 1], y[region_count + i - n+1],y[region_count + i - n],
                              y[region_count - n]), color=mach_color[region_count])

            region_count += 1
        region_count +=1

    # Set color bar
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(weight)
    cb=plt.colorbar(m)
    cb.set_label('Mach [-]')
    if n>10:
        plt.xlim(0,x_max)
    plt.show()


#Streamline through centre line (y=h/2)

y_s_c=np.zeros(2.*n+1.)
x_s_c_i=np.where(y<10E-14)
x_s_c=x[x_s_c_i]
x_s_c=np.hstack(([0],x_s_c))

regions_c1=np.where(abs(phi[:n*(n+1)/2+1+n])<1E-14)
regions_c2=np.where(abs(phi[(n+1)*n+1+n:])<1E-14)
regions_c2=regions_c2+np.ones(len(regions_c2))*((n+1)*n+1+n)

regions_c=np.hstack((regions_c1,regions_c2))
regions_c=regions_c[0]
p_c=np.zeros(len(regions_c))
for i in range(len(regions_c)):
    p_c[i]=P_Pe[regions_c[i]]
print p_c


if plot_pres_c == True:
    for i in range(len(x_s_c)-1):
        plt.plot((x_s_c[i],x_s_c[i+1]),(p_c[i],p_c[i]),'b-')
        if i<len(x_s_c)-2:
            plt.plot((x_s_c[i+1],x_s_c[i+1]),(p_c[i],p_c[i+1]),'b-')
    if n>10:
        plt.xlim(0,x_max)
    plt.show()
if plot_pres_H==True:
    # Streamlines for y=3H/4
    x_s=np.zeros(3*n+2)
    y_s=np.zeros(3*n+2)

    x_s[0]=0.
    y_s[0]=0.5*h

    slopes_phi=np.tan(phi*np.pi/180.)
    # First expansion
    count=1

    for i in range(n):
        x_s[count] = (y[0] - y_s[count - 1] + slopes_phi[count-1] * x_s[count - 1] -
                  slopes_neg[count-1] * x[0]) / (slopes_phi[count-1] - slopes_neg[count-1])
        y_s[count] = slopes_phi[count-1] * (x_s[count] - x_s[count - 1]) + y_s[count - 1]
        plt.plot((x_s[count-1],x_s[count]),(P_Pe[count-1],P_Pe[count-1]),color)
        plt.plot((x_s[count], x_s[count]), (P_Pe[count - 1], P_Pe[count]), color)
        count+=1

    # Through reflection
    for i in range(n):
        if i==0:
            x_s[count] = (y_s[count-1] - y[index_outer_numbers[i]] + slopes_pos[index_gamma_plus_I_II[i]] * x[index_outer_numbers[i]] -
                      slopes_phi[count-1] * x_s[count-1]) / (slopes_pos[index_gamma_plus_I_II[i]] - slopes_phi[count-1])
            y_s[count] = slopes_pos[index_gamma_plus_I_II[i]] * (x_s[count] - x[index_outer_numbers[i]]) + y[index_outer_numbers[i]]
        else:
            x_s[count] = (y_s[count - 1] - y[index_outer_numbers[i]] + slopes_pos[index_gamma_plus_I_II[i]] * x[
                index_outer_numbers[i]] -
                          slopes_phi[index_result_A[i]] * x_s[count - 1]) / (
                         slopes_pos[index_gamma_plus_I_II[i]] - slopes_phi[index_result_A[i]])
            y_s[count] = slopes_pos[index_gamma_plus_I_II[i]] * (x_s[count] - x[index_outer_numbers[i]]) + y[
                index_outer_numbers[i]]

        plt.plot((x_s[count - 1], x_s[count]), (P_Pe[index_result_A[i]-1], P_Pe[index_result_A[i]-1]), color)
        plt.plot((x_s[count], x_s[count]), (P_Pe[index_result_A[i]], P_Pe[index_result_A[i]-1]), color)

        count+=1

    #Through second reflection (simple)
    i=n+1 #Maximum n+1 points in this region
    while x_s[count] < x[-1] and i>1:
        x_s[count] = (y[index_final_B[n-i+1]] - y_s[count - 1] + slopes_phi[index_final_B[n-i+1]] * x_s[count - 1] -
                      slopes_neg[index_final_B[n-i+1]+(i-n-2)] * x[index_final_B[n-i+1]]) / (slopes_phi[index_final_B[n-i+1]] - slopes_neg[index_final_B[n-i+1]+(i-n-2)])
        y_s[count] = slopes_phi[index_final_B[n-i+1]] * (x_s[count] - x_s[count-1]) + y_s[count-1]
        plt.plot((x_s[count - 1], x_s[count]), (P_Pe[index_final_B[n-i+1]], P_Pe[index_final_B[n-i+1]]), color)
        plt.plot((x_s[count], x_s[count]), (P_Pe[index_final_B[n-i+1]-1], P_Pe[index_final_B[n-i+1]]), color)
        count+=1
        i-=1

    print index_final_B
    # Now append the last point
    x_s[-1]=x[-1]
    y_s[-1]=(x[-1]-x_s[-2])*phi_outer[-1]+y_s[-2]
    plt.plot((x_s[-2], x_s[-1]), (P_Pe[index_final_B[-1]+1], P_Pe[index_final_B[-1]+1]), color)
    if n>10:
        plt.xlim(0,x_max)
    plt.show()






