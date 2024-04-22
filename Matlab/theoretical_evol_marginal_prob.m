%% TIME EVOLUTION OF THE THEORETICAL MARGINAL PROBABILITIES
% This script computes the theoretical marginal probabilities a for
% different times (right before and after a measurement and at the
% measurement time) and plot it
%Parameters
L=1;
F=5;
beta=1;
deltatm=10^(-1);
value_t=[399,400,410,499];% vector selects some times of the whole evolution
deltax=0.0001;

%Initial condition: the potential is initially switched off
P_x_0= @(x) 1/(2*L)*x./x;
P_x_1= @(x) 0*x./x;

% Protocol
xdis=linspace(0,2,1000);
theta_1_x = @(x) interp1(xdis-1,Protocol(deltax),x,"linear");
theta_0_x= @(x) 1-theta_1_x(x);

%Vectors
nummed=5;
resol=100;
nt=nummed*resol;
deltat=deltatm/resol;
Nbins=300;
xvec=linspace(-L,L,Nbins);
measure=0;
pc=zeros(1,nummed);

%Evolution
for it=1:nt
    P_x=@(x) P_x_0(x)+P_x_1(x);
    if (mod(deltat*it,deltatm)==0)%the measure operator acts
        measure=measure+1;
        P_x_1=@(x_p) theta_1_x(x_p).*P_x(x_p);
        P_x_0=@(x_p) theta_0_x(x_p).*P_x(x_p);
        pc(measure)= integral(P_x_1,-L,L);
    else % Fokker-Planck evolution
        P_x_1_v=zeros(1,Nbins);
        P_x_0_v=zeros(1,Nbins);
        for in=1:Nbins
            x=xvec(in);
            k_c_1= @(x_p) K_C_1(x,x_p,deltat,F);
            P_x_1_v(in)=integral(@(x_p) k_c_1(x_p).*P_x_1(x_p),-L,L);
            k_c_0= @(x_p) K_C_0(x,x_p,deltat);
            P_x_0_v(in)=integral(@(x_p) k_c_0(x_p).*P_x_0(x_p),-L,L);
        end
        P_x_1 = @(x_p)interp1(xvec,P_x_1_v,x_p,"spline");
        P_x_0 = @(x_p)interp1(xvec,P_x_0_v,x_p,"spline");
    end

    % Ploting the results:

    %Computing P^k(x)
    xvec2=linspace(-L,L,100);
    for k=1:length(value_t)
        if (it==value_t(k))
            hold on
            box on
            ax = gca;
            ax.LineWidth = 2;
            figure(k)
            plot(xvec,P_x_0(xvec)+P_x_1(xvec),'b-','LineWidth',2)
            hold on
            plot(xvec2,P_x_1(xvec2),'g--','LineWidth',2)
            hold on
            plot(xvec2,P_x_0(xvec2),'r.','LineWidth',1)
            hold on
            ylabel('$P$','fontsize',35,'interpreter','latex')
            xlabel('$x$','fontsize',35,'interpreter','latex')
            set(gca,'FontSize',35)
            set(0,'DefaultAxesFontName', 'Times New Roman')
            axis([-1,1,0,1.4])
            shg
        end
    end
end

%Computing P^k(C)
figure(k+1)
box on
ax = gca;
ax.LineWidth = 2;
figure(k+1)
hold on
plot(pc,'o','LineWidth',2,'MarkerSize',10)
ylabel('$P^{k}_c$','fontsize',30,'interpreter','latex')
xlabel('$k$','fontsize',30,'interpreter','latex')
set(gca,'FontSize',30)
set(0,'DefaultAxesFontName', 'Times New Roman')



