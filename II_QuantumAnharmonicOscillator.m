%% Constants & Building the grid
clc
clear
z0 = 1;
zmax = 10*z0;
N = 2^10;   %Grid size (# of points)
dz = zmax/N;
v = 0:1:N-1;
z = -zmax/2 + v*dz;
kmax = 2*pi/dz;
dk = kmax/N;
k = -kmax/2 + v*dk;
wt=0:0.2:600;
zeta = 0.001;   %Damping Amplitude
omega = 0.01;   %Anharmonicity term

%%%%% Initial wavefunction %%%%%
    psi0 = exp(-0.5*((z-2*z0)./z0).^2);
    psi0 = psi0./sqrt(sum(abs(psi0).^2));
%%%%%%%%%%%%%%%%%%%

%% Building Hamiltonian matrix h_m,v & find eigenvalues and eigenvectors
hmv = zeros(N);

for m = 1:N
for v = 1:N  
        if m == v
        hmv(m,v) = (1-1i*zeta)*(      (z0/dz)^2 + 0.5*(z(m)/z0)^2      );
        elseif m+1 == v || m-1 == v
        hmv (m,v) =(1-1i*zeta)* (   -0.5*(z0/dz)^2  );
        end     
end
end

hmv = hmv  + omega*hmv*hmv;

[V,D] = eig(hmv);

%% Time Evolution of the wavefunction by linear superposition of of the time-dependent expansion coofficients
psi=V';
c = psi0*conj(psi');

psit = zeros(N,length(wt));
psit2 = zeros(N,length(wt));

for tt = 1:length(wt)
    for n = 1:length(c)
        psit(:,tt) = psit(:,tt)+ c(n).*psi(n,:)'.*exp(-1i*wt(tt)*D(n,n)); 
    end
    psit(:,tt) = psit(:,tt) ./ sqrt((sum(abs(psit(:,tt)).^2)));
    psit2(:,tt) = abs(psit(:,tt)).^2;
end
%% Fourier transform to k-space for the wavefunction at each time delay
psik = zeros(N,length(wt));

for tt = 1:length(wt) 
psik(:,tt) = ftxtop(psit(:,tt),dz,1);
 psik(:,tt) = psik(:,tt) ./ sqrt((sum(abs(psik(:,tt)).^2)));
end

%% Calculating expectation values <z>,<k> , standard deviations Sigma_z , Sigma_k and uncetainty relation 
Z = zeros(1,length(wt));
SigmaZ = zeros(1,length(wt));
P = zeros(1,length(wt));
K = zeros(1,length(wt));
SigmaP = zeros(1,length(wt));
SigmaK = zeros(1,length(wt));
Unc = zeros(1,length(wt));

for tt = 1:length(wt)
    
Z(tt) = z*(abs(psit(:,tt)).^2);
SigmaZ(tt) = sqrt(((z-Z(tt)).^2)*(abs(psit(:,tt)).^2));
K(tt) = k*(abs(psik(:,tt)).^2);
SigmaK(tt) = sqrt((((k-K(tt)).^2))*(abs(psik(:,tt)).^2));
Unc(tt) = SigmaZ(tt) * SigmaK(tt);
end

%% Plotting

 % Position expectation value <z> vs wt
figure(1)
plot(wt,Z);
set(gca,'FontSize',15);
        xlabel('$wt$','interpreter','latex')
        ylabel('${<Z(t)>}$','interpreter','latex')
        xline(157,'--','${wT/2 = \pi/2\Omega =157}$','FontSize',25,'interpreter','latex');
        xline(314,'--','${wT = \pi/\Omega = 314}$','FontSize',25,'interpreter','latex');
         title('$\zeta=0.01 \ , \ d=2z_0$','interpreter','latex')
         
 % Standard deviation Sigma_z vs wt
 figure(2)
 plot(wt, SigmaZ);
 set(gca,'FontSize',15);
        xlabel('$wt$','interpreter','latex')
        ylabel('${\sigma_z (t)}$','interpreter','latex')
          title('$\zeta=0.01 \ , \ d=2z_0$','interpreter','latex')
        xline(157,'--','${wT/2 = \pi/2\Omega =157}$','FontSize',25,'interpreter','latex');
        xline(314,'--','${wT = \pi/\Omega = 314}$','FontSize',25,'interpreter','latex');

%   Probability density (|Psi(z,t)|^2)
figure(3)
imagesc(wt,z,psit2);
 xlabel('$wt$','interpreter','latex')
        ylabel('${<Z(t)>}$','interpreter','latex')
        colorbar
        set(gca,'FontSize',15);
           title('$\zeta=0.01 \ , \ d=2z_0$','interpreter','latex')

% Wigner Function W(z,t) at wt =157
ind = find(wt == 157);
y = psit(:,ind)';
rho = y'*y;
rho_rotated =imrotate(rho,45,'bicubic','crop');
wigner = real(ftxtop(rho_rotated,dz,1));
wigner = wigner/sum(wigner,'all');

figure(4)
imagesc(z,k,wigner);
colorbar
xlabel('$<Z(t)>$','interpreter','latex')
        ylabel('${<K(t)>}$','interpreter','latex')
set(gca,'FontSize',15);
title('$Wigner \ Function \ of \ a \ cat \ state \ \zeta=0.001 \ , \ d=2z_0$','interpreter','latex')

