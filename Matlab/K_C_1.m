function k_x_1 = K_C_1(x,x_p,t,F)
% This function computes the Green function for V-shape potential
L=1;
beta=1;
Nsum=20;
phi0_p=F*beta/(2*(1-exp(-F*beta*L)))*exp(-F*beta*abs(x_p));
phi0=F*beta/(2*(1-exp(-F*beta*L)))*exp(-F*beta*abs(x));
k_x_1=phi0;
for n=1:Nsum
    lambdan=F^2*beta/4+n^2*pi^2/(L^2*beta);
    % Antisymmetric eigenfunction
    phi1n=sqrt(phi0/L).*sin(n*pi*x/L);
    phi1n_p=sqrt(phi0_p/L).*sin(n*pi*x_p/L);
    % Symmetric eigenfunction
    phi2n=1/2*sqrt(phi0/(L*beta*lambdan)).*(2*n*pi/L*cos(n*pi*x/L)-F*beta*sin(n*pi*abs(x)/L));
    phi2n_p=1/2*sqrt(phi0_p/(L*beta*lambdan)).*(2*n*pi/L*cos(n*pi*x_p/L)-F*beta*sin(n*pi*abs(x_p)/L));
    k_x_1= k_x_1+ exp(-lambdan*t)*1./(phi0_p).*(phi1n.*phi1n_p+phi2n.*phi2n_p);
end
end