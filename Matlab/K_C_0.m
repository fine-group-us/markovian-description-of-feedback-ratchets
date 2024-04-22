function k_x_0 = K_C_0(x,x_p,t)
% This function computes the Green function for V(x)=0
L=1;
beta=1;
Nsum=20;
phi0_p=1/(2*L);
phi0=1/(2*L);
k_x_0=phi0;
for n=1:Nsum
    lambdan=1*n^2*pi^2/(beta*L^2);
    % Antisymmetric eigenfunction
    phi1n=sqrt(1/2)*1/L*sin(n*pi*x/L);
    phi1n_p=sqrt(1/2)*1/L*sin(n*pi*x_p/L);
    % Symmetric eigenfunction
    phi2n=sqrt(1/2)*1/L*cos(n*pi*x/L);
    phi2n_p=sqrt(1/2)*1/L*cos(n*pi*x_p/L);
    k_x_0= k_x_0+ exp(-lambdan*t)*1./(phi0_p).*(phi1n.*phi1n_p+phi2n.*phi2n_p);
end
end