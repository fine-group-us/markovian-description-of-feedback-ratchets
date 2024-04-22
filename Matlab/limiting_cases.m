%% P(X,t) IN LIMITING CASES
% This script computes the theoretical marginal probability of the particle
% position and load the same data from simulations in the limiting cases
% \Delta t_m tending to 0 and infinity, and plot both.
%% Parameters
L=1;
F=5;
beta=1;
deltax=0.001;
Nbins=10000;
xvec=linspace(-L,L,Nbins);
%% Theoretical predection
% Protocol
xdis=linspace(0,2,1000);
theta_1_x = @(x) interp1(xdis-1,Protocol(deltax),x,"linear");
theta_0_x= @(x) 1-theta_1_x(x);
% Conditional probabilities
P1 = @(x) exp(- beta*F*abs(x));
norm1=integral(@(x) exp(- beta*F*abs(x)),-L,L);
P1 = @(x) P1(x)/norm1;
P0 = @(x) x./(2*x);

%% For \delta t_m tending to infinity
% Computing Bcc coefficients
b10= integral( @(x) P1(x).*theta_0_x(x) ,-L,L);
b01= integral( @(x) P0(x).*theta_1_x(x),-L,L);
b00=integral( @(x) P0(x).*theta_0_x(x),-L,L);
b11=integral( @(x) P1(x).*theta_1_x(x),-L,L);

% Computing p1=p(C=1)
p1=b10/(1-b11+b10);
p0=1-p1;

% Computing fsinfty=\rho^s(x)
fsinfty=@(x) p0*1/(2*L)+p1*exp(- beta*F*abs(x))/norm1;

%% For \delta t_m tending to zero

% Computing the effective force Feff
Feff= @(x) -sign(x)*F.*theta_1_x (x);
Veff= @(x) -integral(@(y) Feff(y),0,x);
Veffvec=zeros(1,Nbins);
for i=1:Nbins
    Veffvec(i)=Veff(xvec(i));
end
Veff= @(x) interp1(xvec,Veffvec,x,"linear");

% Ploting the effective potential
figure(1)
hold on
plot(xvec, Veffvec,'-', 'LineWidth',2)
shg
box on
ylabel('$V_{\\text{eff}}$','fontsize',20,'interpreter','latex')
xlabel('$x$','fontsize',20,'interpreter','latex')
set(gca,'FontSize',20)
set(0,'DefaultAxesFontName', 'Times New Roman')
hold off

% Computing fs0=\rho^s(x)
fs0=@(x) exp(- beta*Veff(x))+(exp(-beta*Veff(-L))-exp(-beta*Veff(L)))./(integral(@(x) exp( -beta*(Veff(L)-Veff(x))),-L,L))*integral(@(y) exp(- beta*(Veff(x)-Veff(y))),-L,x);
valor1=zeros(1,Nbins);
for i=1:Nbins
    z=xvec(i);
    valor1(i)=fs0(z);
end
fs0=@(x) interp1(xvec,valor1,x,"linear");
norm0=integral(fs0,-L,L);
fs0= @(x) fs0(x)/norm0;


%% Simutaliton resuslts

% Loading the data:

% For \Delta t_m tending to 0
nummed=5000;
resol=10;
nt_0=nummed*resol;
PC1x_0=load('datos_de_P_C_1_0.dat');
PC0x_0=load('datos_de_P_C_0_0.dat');
dimensions_0=size(PC1x_0);
factor_0= dimensions_0(1)/nt_0;
Binsim_0=dimensions_0(2);

% For \Delta t_m tending to infinity
nummed=5;
resol=100;
nt_infty=nummed*resol;
PC1x_infty=load('datos_de_P_C_1_infty.dat');
PC0x_infty=load('datos_de_P_C_0_infty.dat');
dimensions_infty=size(PC1x_infty);
factor_infty= dimensions_infty(1)/nt_infty;
Binsim_infty=dimensions_infty(2);

% Ploting the limiting distribution for x-variable
figure(2)
yvalue0=(PC0x_0(nt_0*factor_0,:)+PC1x_0(nt_0*factor_0,:))*(dimensions_0(2)/(2*L));
xvalue0=linspace(-L,L,dimensions_0(2));
yvalue0=yvalue0(1:2:dimensions_0(2));
xvalue0=xvalue0(1:2:dimensions_0(2));
box on
hold on
scatter(xvalue0,yvalue0,70,'s', 'MarkerEdgeColor', 'b','MarkerFaceColor', 'none',LineWidth=0.75);
hold on
yvalueinfty=(PC0x_infty(nt_infty*factor_temporal_infty,:)+PC1x_infty(nt_infty*factor_temporal_infty,:))*(dimensiones_infty(2)/(2*L));
xvalueinfty=linspace(-L,L,dimensiones_infty(2));
yvalueinfty=yvalueinfty(1:4:dimensiones_0(2));
xvalueinfty=xvalueinfty(1:4:dimensiones_0(2));
scatter(xvalueinfty,yvalueinfty,70,'o', 'filled', 'MarkerFaceColor', 'b');
hold on
plot(xvec,fs0(xvec),'k-', 'LineWidth',2)
plot(xvec,fsinfty(xvec),'k-','LineWidth',2)
ylabel('$ P(x,t_k^-)$','fontsize',20,'interpreter','latex')
xlabel('$x$','fontsize',20,'interpreter','latex')
set(gca,'FontSize',20)
set(0,'DefaultAxesFontName', 'Times New Roman')
