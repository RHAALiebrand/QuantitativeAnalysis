clear all
close all
clc
 
% The system that you need to solve will be singular. Matlab gives you a
% warning at each time step. To switch of this warning, remove the comment
% in the next line
 
warning off
 
% This file contains the skeleton of the program which will solve the lid
% driven cavity problem on a unit square. The parts that have to be
% supplemented are described in the assignment.
%
% The pieces that need to be filled in are indicated by dots: ............
%
 
%
% When running the code, determine a suitable time step. A too small time
% step will make the calculation very long, while for a too large time step
% the solution will blow up due to numerical instability.
%

iter_max=100000000;
tol = 1e-9;    % tol determines when steady state is reached and the program terminates
%dt = 0.0048    % time step

Re = 1000;             % Reynolds number
N =48;                 % Number of volumes in the x- and y-direction
Delta = 1/N;           % uniform mesh width in the mapping to a cosine grid
t = 0.0001;            % time
             
% wall velocities
U_wall_top = -1;
U_wall_bot =0;
U_wall_left = 0;
U_wall_right = 0;
V_wall_top = 0;
V_wall_bot = 0;
V_wall_left = 0;
V_wall_right = 0;
 
%
%   Generation of a non-uniform mesh
%
 
%
%   tx are the coordinates of the nodal points on the outer-oriented grid.
%   tx is an abbreviation for \tilde{x}
%
tx = zeros(1,N+1);
for i=1:N+1
    xi = (i-1)*Delta;
    tx(i) = 0.5*(1. - cos(pi*xi));
end 

% Local mesh size on outer oriented grid
th = tx(2:N+1) - tx(1:N);
%  x are the coordinates of the nodal points on the inner-oriented grid  (including
%  endpoints 0 and 1)
%  h contains the edge lengths on the inner-oriented grid
x = 0.5*(tx(1:N) + tx(2:N+1));
x = [0 x 1];
h = x(2:N+2) - x(1:N+1);
%   Initial condition u=v=0
%
%   Both u and v will be stored in one big vector called 'u'
%
%   The vector u only contains the true unknowns, not the velocities that
%   are prescribed by the boundary conditions
%
%   The vector u contains the *inner-oriented* fluxes as unknowns
%
u = zeros(2*N*(N-1),1);
% Set up the Incindence matrix 'tE21' which connects the fluxes to the
% volumes. Use the orientation described in the assignment.
tE21=sparse(N^2,2*N*(N+1));  % Number of edges is determined to be 2N(N+1)
tE21pattern=[-1 zeros(1,N-1) -1 1 zeros(1,N-1) 1];
edge=0;
cell=0;
for i=[1:N]
    for j=[1:N]
        edge=edge+1;
        cell=cell+1;
        tE21(cell,edge:edge+length(tE21pattern)-1)=tE21pattern;   
    end
    edge=edge+(N+1);
end
full_tE21=tE21;
%  Inserting boundary conditions for normal velocity components
tbndu=zeros(4*N,1);  % Set vector containing all boundary conditions
tbndmat=sparse(N^2,4*N);  % Set matrix containing boundary edges

tbndu(1:N,1)=V_wall_bot.*th;
tbndmat(:,1:N)=tE21(:,1:N); % Replace corresponding columns in tbndmat 
tbndu(end-N+1:end,1)=V_wall_top.*th;
tbndmat(:,end-N+1:end)=tE21(:,end-N+1:end); % Replace again
% Now we have inserted the first and last boundaries, no the middel ones
Eleft=N+1; % Count in tbndu
Eright=N+2;
subleft=N+1;  % Count in tE21
subright=2*N+1;
cols2remove=[];
for i=[1:N] % Loop over rows with side element
    tbndu(Eleft,1) = U_wall_left*th(i);
    tbndmat(:,Eleft)=tE21(:,subleft);
    
    tbndu(Eright,1)= U_wall_right*th(i);
    tbndmat(:,Eright)=tE21(:,subright);
    
    cols2remove=[cols2remove subleft subright];
    Eleft= Eleft + 2;
    Eright= Eright + 2;
    subleft=subleft + 2* N +1 ;
    subright=subright + 2* N +1 ;
end
tu_norm=tbndmat*tbndu; % Calculate vector which represents mass flow out of cells
%   Remove columns associated with prescribed normal velocities from
%   Incidence matrix 'tE21'
tE21(:,cols2remove) = [] ;% Remove left right
tE21(:,1:N)=[]; % Remove bottom
tE21(:,end-N+1:end) = []; % Remove top
% Setting up simple Hodge matrix which converts the fluxes on the
% outer-oriented grid to circulation on the inner-oriented grid. Assume
% that the fluxes and circulation are constant over each 1-cell. This will
% give a diagonal Hodge matrix. Call this Hodge matrix 'H1t1'
H1t1=[];
for i=[1:N]
    H1t1 = [H1t1 h(i)./th];
    H1t1 = [H1t1 h./th(i)];
end
H1t1 = [H1t1 H1t1(1:N)];
full_H1t1=sparse(diag(H1t1));
full_Ht11=sparse(diag(1./H1t1));
Hu_norm = zeros(2*N*(N+1),1); % Initially contains the total number of edges
% Hu_norm is the vector which will contain the Hodge of the prescribed
% normal fluxes. Calculate the vector 'Hu_norm'
Hu_norm(1:N)=H1t1(1:N); % Insert  bottom
Hu_norm(end-N+1:end)=H1t1(end-N+1:end); % Insert top
% Now start appending left and right side
subleft=N+1;
subright=2*N+1;
for i=[1:N]
    Hu_norm(subleft)=H1t1(subleft);
    Hu_norm(subright)=H1t1(subright);
    subleft=subleft + 2* N +1;
    subright=subright + 2* N +1;
end
%  Remove corresponding row and columns from the Hodge matrix and also
%  remove the corresp0onding 'rows' from Hu_norm
intindex=find(Hu_norm==0) ;
bndindex=find(Hu_norm ~= 0)  ;
% Removing boundaries edges and vise versa
Hu_norm(intindex)=[];  % Set interiors to zero
H1t1(bndindex)=[]; % Set boundaries to zero

% Now calculate u_norm
Hu_norm=sparse(diag(Hu_norm));
u_norm=Hu_norm*tbndu;
% Set up the incidence E^{2,1} between 1-cochain circulation and 2-cochain vorticity on
% the inner-oriented (extended) grid
% This incidence matrix will be called 'E21' in the program
E21=sparse((N+1)^2,2*(N+2)*(N+1));  % Number of edges is determined to be 2N(N+1)
E21pattern=[1 zeros(1,N) -1 1 zeros(1,N) -1];
edge=0;
cell=0;
for i=[1:N+1]
    for j=[1:N+1]
        edge=edge+1;
        cell=cell+1;
        E21(cell,edge:edge+length(E21pattern)-1)=E21pattern;   
    end
    edge=edge+(N+2);
end
full_E21=E21;
% Inserting prescribed tangential boundary conditions
% Initialize vector, matrix and cols2remove
bndu=zeros(4*(N+1)+4*N,1);  % Set vector containing all boundary conditions
bndmat=sparse((N+1)^2,4*(N+1)+4*N);  % Set matrix containing boundary edges
cols2remove=[];

% Insert bot boundary
bndu(1:(N+1),1)=h.*U_wall_bot;
bndmat(:,1:(N+1))=E21(:,1:(N+1));

% Insert top boundary
bndu(end-N:end,1)=h.*U_wall_top;
bndmat(:,end-N:end)=E21(:,end-N:end);

Eleft=N+2; 
Eright=2*(N+1)+1;
subleft=N+2;
subright=2*(N+1)+1;

% First two tangential velocities (5 and 9 in n=3 case)
bndu(Eleft,1) = h(1)*V_wall_left;
bndmat(:,Eleft)=E21(:,subleft);   
bndu(Eright,1)=h(1)*V_wall_right;
bndmat(:,Eright)=E21(:,subright);
cols2remove=[cols2remove subleft subright];

% Last two tangential velocities (32 and 36)
bndu(end-Eleft+1,1)=h(end)*V_wall_right;
bndu(end-Eright+1,1)=h(end)*V_wall_left;
bndmat(:,end-Eleft+1)=E21(:,end-subleft+1);
bndmat(:,end-Eright+1)=E21(:,end-subright+1);
cols2remove=[cols2remove length(E21(1,:))-subleft+1 length(E21(1,:))-subright+1];

% Pick out inner tangentails
Eleft=Eright+3;
Eright=Eright+4;
subleft=subleft+2*(N+1)+1;
subright=subright+2*(N+1)+1;

for i=[1:N-1]
    bndu(Eleft,1) = h(i+1)*V_wall_left;
    bndmat(:,Eleft)=E21(:,subleft);
    
    bndu(Eright,1)= h(i+1)*V_wall_right;
    bndmat(:,Eright)=E21(:,subright);
    
    cols2remove=[cols2remove subleft subright];
    Eleft= Eleft + 4;
    Eright= Eright + 4;
    subleft=subleft + 2* (N+1) +1 ;
    subright=subright + 2* (N+1) +1 ;
end

% Now input the normal velcoities
% Insert bot normal boundaries
bndu(N+3:N+3+N-1,1)=u_norm(1:N);
bndmat(:,N+3:N+3+N-1) = E21(:,N+3:N+3+N-1);
% Insert top boundary
bndu(end-(N+3+N-1)+1:end-(N+3)+1,1)=u_norm(end-N+1:end);
bndmat(:,end-(N+3+N-1)+1:end-(N+3)+1)=E21(:,length(E21(1,:))-2*N-1:length(E21(1,:))-N-2);

Eleft=N+1+N+2+1;
Eright=Eleft+1;
subleft=Eleft;
subright=subleft+N;
normindex=N+1;
for i=[1:N]
    bndu(Eleft,1)=u_norm(normindex);
    bndu(Eright,1)= u_norm(normindex+1);
    
    bndmat(:,Eleft)=E21(:,subleft);
    bndmat(:,Eright)=E21(:,subright);
    cols2remove=[cols2remove subleft subright];
    
    normindex=normindex+2;
    subleft=subleft+2*(N+1)+1;
    subright=subright+2*(N+1)+1;
    Eleft=Eleft+4;
    Eright=Eright+4;
end

% Remove columns from the incidence matrix E21 corresponding to both the
% prescribed tangential velocities and normal velocities
E21(:,cols2remove) = [] ;% Remove left right
E21(:,1:N+1)=[]; % Remove bottom
E21(:,end-N:end) = []; % Remove top
% FIXME: Remove columns responsible for first and last 3 normals
E21(:,1:N)=[] ; % Remove normals bot
E21(:,end-N+1:end)=[];
% Store the prescribed normal and tangential velocities in the vector
% 'xi_pres'
xi_pres=bndmat*bndu;
dlmwrite(strcat('xipres_N=',int2str(N),'.dat'),xi_pres)
% Set up the Hodge matrix which maps inner-oriented 2-cochains to
% outer-oriented 0-cochains. Assume that the vorticity is constant in the
% inner-oriented 2-cells. This will give a diagonal Hodge matrix. Call this
% Hodge matrix 'Ht02'
Ht02=[];
for i=[1:N+1]
    Ht02=[Ht02 1./(h.*h(i))];
end
Ht02=sparse(diag(Ht02));
% Set up the Hodge matrix which maps inner-oriented 1-cochains to
% outer-oriented 1-cochains. Call this Hodge matrix 'Ht11'. Assume again
% that everything is constant over the 1-cells, which will then yield a
% diagonal Hodge matrix.
Ht11=1./H1t1;
H1t1=sparse(diag(H1t1));
Ht11=sparse(diag(Ht11));
% The prescribed velocities will play a role in the momentum equation
xi_pres = H1t1*E21'*Ht02*xi_pres;

% Now all matrices are set up and the time stepping can start. 'iter' will
% record the number of time steps. This allows you to give output after a
% preselected number of time steps.
%
% 'diff' will be the maximal du/dt or dv/dt. If 'diff' is sufficiently
% small, steady state has been reached. Determine a suitable value for
% 'tol'
%
 
diff = 1;
iter = 1;
A = -tE21*Ht11 *tE21';
dt=min(min(h),0.5*Re*min(h)^2);
% while (diff > tol && iter <= iter_max)
% %       Calculate the convective terms using the vorticity and the local
% %       velocity field. Store the convective terms in a vector
% %       'convective'.
%     
% %      Note that you only have to set up the convective terms for true
% %      velocity unknowns and not for those inner-oriented circulations for
% %      which the value is already known.
%     
%     xi=E21*u;
%     txi=Ht02*xi;     %Vorticity in point of primal mesh
% %     Start for y direction, so integrate over dx
% %     -v*xi*dx
% %     Averaging to orange points
%     convection=zeros(length(u),1);
%     vavg=zeros((N+1),(N-1));
%     count=N;
%     countxi=N+1+2;
%     for i=[1:N+1] % Loop over rows
%         if i==1
%             for j=[1:N-1] % Loop over columns
%                 vavg(i,j)=-(u_norm(j)+u_norm(j+1))/(2*h(i))*txi(j+1);
%             end
%         elseif i==N+1    
%             for j=[1:N-1] % Loop over columns
%                 vavg(i,j)=-(u_norm(N+2*N+j)+u_norm(N+2*N+j+1))/(2*h(i))*txi(N*(N+1)+1+j);
%             end
%         else
%              for j=[1:N-1] % Loop over columns    
%                 vavg(i,j)=-(u(count+j-1)+u(count+j))/(2*h(i))*txi(countxi+j-1);
% %                 i,j
% %                 left=count+j-1
% %                 right=count+j
% %                 mid=countxi+j-1
%              end 
%              count=count+N+N-1;
%              countxi=countxi+N-1+2;
%         end     
%     end
%     convcount=1;
%     for i=[1:N]
%         for j=[1:N-1]
%             convection(convcount+j-1)=(vavg(i+1,j)+vavg(i,j))/2*h(j+1);
% %             lower=i,j
% %             upper= i+1,j
% %             conv=convcount+j-1
%         end 
%         convcount=convcount+N-2+N+1;
%     end
%     
%     
%     uavg=zeros((N+1),(N-1));
%     for i=[1:N+1] % Loop over rows
%         if i==1
%             normcount=N+1;
%             countxi=N+2;
%             for j=[1:N-1] % Loop over columns
%                 uavg(i,j)=(u_norm(normcount)+u_norm(normcount+2))/(2*h(i))*txi(countxi);
% %                 i,j
% %                 top=normcount+2
% %                 bot=normcount
% %                 mid=countxi
%                 normcount=normcount+2;
%                 countxi=countxi+N+1;
%             end
%          elseif i==N+1
%              normcount=N+2;
%              countxi=2*(N+1);
%              for j=[1:N-1] % Loop over columns
%                  uavg(i,j)=(u_norm(normcount)+u_norm(normcount+2))/(2*h(i))*txi(countxi);
% %                  i,j
% %                  bot=normcount
% %                  top=normcount+2
% %                  mid=countxi
%                  normcount=normcount+2;
%                  countxi=countxi + N + 1;
%              end
%         else    
%               count = i-1;
%               deltacount= 2*N-1;
%               countxi=i+N+1;
%               deltaxi=N+1;
%               for j=[1:N-1] % Loop over columns    BUT WHICH U TO USE??? U of U_pres??
%                  uavg(i,j)=(u(count)+u(count+deltacount))/(2*h(i))*txi(countxi);
% %                  i,j
% %                  bot=count
% %                  top=count+deltacount
% %                  mid=countxi
%                  count=count+deltacount;
%                  countxi=countxi+deltaxi;
%               end 
%         end     
%     end
%     
%     convcount=N;
%     for i=[1:N-1]
%         for j=[1:N]
%             convection(convcount+j-1)=(uavg(j,i)+uavg(j+1,i)) / 2*h(i+1);
% %             conv=convcount+j-1
% %             left=j,i
% %             right=j+1,i
% %             lenggth=i+1
%         end
%         convcount=convcount+2*N-1;
%     end
% %     Set up the right hand side for the Poisson equation for the pressure
%     
%     rhs_Poisson  =   tE21*Ht11*(u/dt  - convection - H1t1*E21'*Ht02*E21*u/Re - xi_pres/Re) + tu_norm/dt;
%     
% %     Set up the matrix for the Poisson equation
%     
% %     Folve for the new pressure
%     
%     p = A\rhs_Poisson;
%     
% %     Store the velocity from the previous time step in the vector u_old
%     
%     uold = u;
%     
% %     Update the velocity field
%     
%     u = u - dt* (convection - tE21'*p + H1t1 *E21'*Ht02*E21*u/Re + xi_pres/Re);
%     
%     
% %      Every other 1000 iterations check whether you approach steady state
% %      and check whether you satisfy conservation of mass. The largest
% %      rate at which mass is destroyed or created is denoted by 'maxdiv'.
% %      This number should be very small, in the order of machine precision.
%     
%     if mod(iter,1000) == 0
%         maxdiv = max(tE21*Ht11*u + tu_norm);
%         diff = max(abs(u-uold))/dt
%         
%     end
%         
%     iter = iter + 1;
% end    
% %  Produce desired output to compare your results with those given in the
% %  reference by Botella and Peyret
% dlmwrite(strcat('u_N=',int2str(N),'.dat'),u)
% dlmwrite(strcat('xi_N=',int2str(N),'Re=',int2str(Re),'.dat'),xi)
% dlmwrite(strcat('txi_N=',int2str(N),'.dat'),txi)
% dlmwrite(strcat('p_N=',int2str(N),'.dat'),p)

u=dlmread(strcat('u_N=',int2str(N),'.dat'));
xi=dlmread(strcat('xi_N=',int2str(N),'.dat'));
txi=dlmread(strcat('txi_N=',int2str(N),'.dat'));
p=dlmread(strcat('p_N=',int2str(N),'.dat'));

v_inner=zeros(2*N*(N-1),1);
vx_inner=[];
vy_inner=[];
x_inner=zeros(2*N*(N-1),1);
xx_inner=[];
xy_inner=[];
y_inner=zeros(2*N*(N-1),1);
yy_inner=[];
yx_inner=[];
count=1;
for i=[1:N]
    v_inner(count:count+N-2)=u(count:count+N-2)./h(2:N)';
    vx_inner=[vx_inner v_inner(count:count+N-2)'];
    x_inner(count:count+N-2)=(x(3:N+1)+x(2:N)).*0.5;
    xx_inner=[xx_inner x_inner(count:count+N-2)'];
    y_inner(count:count+N-2)=x(i+1);
    xy_inner=[xy_inner y_inner(count:count+N-2)'];
    count=count+2*N-1;
end

count=N;
for i=[1:N-1]
    v_inner(count:count+N-1)=u(count:count+N-1)./h(i+1)';
    vy_inner=[vy_inner v_inner(count:count+N-1)'];
    x_inner(count:count+N-1)=x(2:N+1);
    yx_inner=[yx_inner x_inner(count:count+N-1)'];
    y_inner(count:count+N-1)=(x(i+2)+x(i+1))*0.5;
    yy_inner=[yy_inner y_inner(count:count+N-1)'];
    count=count+2*(N-1)+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour values
cont_val_stream=[-1.5e-3 -1e-3 -5e-4 -2.5e-4 -1.0e-4 -5.0e-5 -1e-5 -1e-6 0.0 1e-10 1e-5 1e-4 1e-2 3e-2 5e-2 7e-2 9e-2 0.1 0.11 0.115 0.1175];
cont_val_vort=[-3.0 -2.0 -1.0 -0.5 0.0 1.0 2.0 3.0 4.0 5.0];
cont_val_pres=[-0.002 0.0 0.02 0.05 0.07 0.09 0.11 0.12 0.17 0.3 ];

% Stream function
[tX,tY]=meshgrid(tx,tx);
psi = (E21'\Ht11*u)';
psi_martin = psi;
psi=flipud(reshape(psi,[N+1,N+1])');
figure
contour(tX,flipud(tY),psi,cont_val_stream);
colorbar
xlabel('$x$','fontsize',14,'Interpreter','LaTex')
ylabel('$y$','fontsize',14,'Interpreter','LaTex')
print(strcat('plots/psi_N',int2str(N)),'-depsc')

% Vorticity
txi=flipud(reshape(txi,[N+1,N+1])');
figure;
contour(tX,flipud(tY),txi,cont_val_vort);
colorbar
xlabel('$x$','fontsize',14,'Interpreter','LaTex')
ylabel('$y$','fontsize',14,'Interpreter','LaTex')
% xlim([0.95 1])
% ylim([0.95 1])
print(strcat('plots/txi_N',int2str(N)),'-depsc')

% Pressure
Fx=scatteredInterpolant(xx_inner',xy_inner',vx_inner');
Fy=scatteredInterpolant(yx_inner',yy_inner',vy_inner');
[X,Y]=meshgrid(x(2:N+1),x(2:N+1));
Vx=flipud(Fx(X,Y));
Vy=flipud(Fy(X,Y));
p_dual=flipud(reshape(p,[N,N])');
p_res=p_dual-0.5*(Vx.^2+Vy.^2);
p_ref=(p_res(N/2,N/2)+p_res(N/2,N/2+1)+p_res(N/2+1,N/2)+p_res(N/2+1,N/2+1))*0.25;
p_res=p_res-p_ref;

figure
contour(X,flipud(Y),p_res,cont_val_pres);
colorbar
xlabel('$x$','fontsize',14,'Interpreter','LaTex')
ylabel('$y$','fontsize',14,'Interpreter','LaTex')
xlim([0 1])
ylim([0 1])
print(strcat('plots/p_N',int2str(N)),'-depsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vertical centre line horizontal velocities
% Reference data
y_ref = [0.0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5000 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1.000];
vh_ref = [0.000 0.1812881 0.2023300 0.2228955 0.3004561 0.3885691 0.2803696 0.1081999 0.0620561 -0.0570178 -0.1886747 -0.3372212 -0.4723329 -0.5169277 -0.5808359 -0.6644227 -1.0000];
p_v_ref = [0.110591 0.109689 0.109200 0.108566 0.104187 0.081925 0.040377 0.004434 0.00000 -0.000827 0.012122 0.034910 0.050329 0.050949 0.051514 0.052009 0.052987];
xi_v_ref=[-4.16648 -2.44960 -2.31786 -2.20175 -1.63436 1.05467 2.26772 2.06215 2.06722 2.06539 2.09121 1.76200 4.85754 6.95968 9.49496 12.0570 14.7534];

vx_v=zeros(N,1);
count=N/2;
for i=[1:N]
    vx_v(i)=v_inner(count);
    count=count+N/2-1+N+N/2;
end

% Pressure
Fp=scatteredInterpolant(reshape(X,(N)^2,1),reshape(Y,(N)^2,1),reshape(p_res,(N)^2,1));
p_v=Fp(0.5.*ones(length(x),1),x');
figure
plot(flipud(p_v),x',p_v_ref,y_ref)
xlabel('$p$','fontsize',14,'Interpreter','LaTex')
ylabel('$y$','fontsize',14,'Interpreter','LaTex')
legend('Program','Reference');
print(strcat('plots/p_v_N',int2str(N)),'-depsc')


% Horizontal velocity
figure
plot(vx_v,x(2:N+1),vh_ref,y_ref)
legend('Program','Reference')
xlabel('$v_x$','fontsize',14,'Interpreter','LaTex')
ylabel('$y$','fontsize',14,'Interpreter','LaTex')
print(strcat('plots/vx_v_N',int2str(N)),'-depsc')

% Vertical velocity
vy_v=Fy(0.5.*ones(length(x),1),x');
figure
plot(vy_v,x')
xlabel('$v_y$','fontsize',14,'Interpreter','LaTex')
ylabel('$y$','fontsize',14,'Interpreter','LaTex')
print(strcat('plots/vy_v_N',int2str(N)),'-depsc')

% Vorticity
Fxi=scatteredInterpolant(reshape(tX',(N+1)^2,1),reshape(tY',(N+1)^2,1),reshape(flipud(txi)',(N+1)^2,1));
xi_v=Fxi(0.5.*ones(length(tx(2:end-1)),1),tx(2:end-1)');
figure
plot(xi_v,tx(2:end-1)',xi_v_ref,y_ref)
legend('Program','Reference','Location','southeast')
xlabel('$\xi$','fontsize',14,'Interpreter','LaTex')
ylabel('$y$','fontsize',14,'Interpreter','LaTex')
print(strcat('plots/xi_v_N',int2str(N)),'-depsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Horizontal centre line
x_ref = [0.000 0.0312 0.0391 0.0469 0.0547 0.0937 0.1406 0.1953 0.5000 0.76556 0.7734 0.8437 0.9062 0.9219 0.9297 0.9375 1.0000];
vv_ref = [0.00 -0.2279225 -0.2936869 -0.3553213 -0.4103754 -0.5264392 -0.4264545 -0.3202137 0.0257995 0.3253592 0.3339924 0.3769189 0.3330442 0.3099097 0.2962703 0.2807056 0.0000];
p_h_ref = [0.077455 0.078837 0.078685 0.078148 0.077154 0.065816 0.049029 0.034522 0.0000 0.044848 0.047260 0.069511 0.084386 0.086716 0.087653 0.088445 0.090477];
xi_h_ref=[-5.46217 -8.44350 -8.24616 -7.58524 -6.50867 0.92291 3.43016 2.21171 2.06722 2.06122 2.00174 0.74207 -0.82398 -1.23991 -1.50306 -1.83308 -7.66369];

%Vertical velocities
count = (N-1)*N/2 + N*(N/2-1)+1;
vy_h= v_inner(count:count+N-1);
figure;
plot(x(2:N+1),vy_h,x_ref,vv_ref);
legend('Program','Reference','Location','southeast')
xlabel('$x$','fontsize',14,'Interpreter','LaTex')
ylabel('$v_y$','fontsize',14,'Interpreter','LaTex')
print(strcat('plots/vy_h_N',int2str(N)),'-depsc')

% Horizontal velocity
figure
vx_h=Fx(x',0.5.*ones(length(x),1));
plot(x,vx_h)
xlabel('$x$','fontsize',14,'Interpreter','LaTex')
ylabel('$v_x$','fontsize',14,'Interpreter','LaTex')
print(strcat('plots/vx_h_N',int2str(N)),'-depsc')

% Pressure
p_h=Fp(x',0.5.*ones(length(x),1));
figure;
plot(x,p_h,x_ref,p_h_ref);
legend('Program','Reference','Location','southeast')
xlabel('$x$','fontsize',14,'Interpreter','LaTex')
ylabel('$p$','fontsize',14,'Interpreter','LaTex')
print(strcat('plots/p_h_N',int2str(N)),'-depsc')

% Vorticity
xi_h=Fxi(x',0.5.*ones(length(x),1));
figure;
plot(x,xi_h,x_ref,xi_h_ref);
legend('Program','Reference','Location','southeast')
xlabel('$x$','fontsize',14,'Interpreter','LaTex')
ylabel('$\xi$','fontsize',14,'Interpreter','LaTex')
print(strcat('plots/xi_h_N',int2str(N)),'-depsc')

[tX,tY]=meshgrid(tx,tx);
Fpsi = scatteredInterpolant(reshape(tX,[(N+1)^2 1]),reshape(tY,[(N+1)^2 1]),psi_martin');
Fpsi(0.136,0.1118);

xi_N=[];
for i=[16 32 48 56 64]
   xi=dlmread(strcat('xi_N=',int2str(i),'.dat'));
   xipres=dlmread(strcat('xipres_N=',int2str(i),'.dat')); 
   xi_N=[xi_N sum(xi-xipres)];
end


figure
plot([16 32 48 56 64],xi_N)
xlabel('N','fontsize',14)
ylabel('$\int_\Omega \xi d\Omega$','fontsize',14,'Interpreter','LaTex')

xi_Re=[];
for i=[250 500 750 1000 1250 1500 1750 2000 2250 2500 3125 3750 4375 5000 5625 6250 6875 7500 8125 9375 10000]
    xi=dlmread(strcat('xi_N=',int2str(16),'Re=',int2str(i),'.dat'));
    length(xi);
    xipres=dlmread(strcat('xipres_N=',int2str(16),'.dat'));
    length(xipres);
    xi-xipres;
    xi_Re=[xi_Re sum(xi-xipres)];
end 

figure
plot([250 500 750 1000 1250 1500 1750 2000 2250 2500 3125 3750 4375 5000 5625 6250 6875 7500 8125 9375 10000],xi_Re)
xlabel('Re','fontsize',14)
ylabel('$\int_\Omega \xi d\Omega$','fontsize',14,'Interpreter','LaTex')