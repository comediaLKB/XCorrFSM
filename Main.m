%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is used to generate a thick forward scattering
% diffuser and calculate the cross-correlation product for
% two wavelengths lambda_1 and lambda_2 which pass through
% the diffuser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code is written by Payvand Arjmand & Marc Guillon
% Analytical models from Zhu et al. Optica 7(4), 338 (2020)
% https://doi.org/10.1364/OPTICA.382209
% 
% Contact info: 
% payvand.arjmand@u-paris.fr
% marc.guillon@u-paris.fr

clear

%% Parameters
%physical parameters 
lambda_1=0.8; %wavelength µm
lambda_2=0.81; %wavelength µm
k_1=2*pi/lambda_1; %wavenumber
k_2=2*pi/lambda_2; %wavenumber2
ls=100;%scattering mean free path (µm)
L=1000;%diffuser thickness (µm)
g=1-1E-3; %Anisotropy factor

%numerical parameters
N=512; %matrix dimension (square)
Nz=32; %number of z computed values (longitudinal coordinate)
N_diff=100; %number of diffuser layers
px=1; %pixel size (µm)

d=L/(N_diff-1); %distance between diffuser layers
Theta_0=sqrt(L*(1-g)/ls); % Theoretical output scattering angle
disp(['Output scattering angle (theoretical) : ',num2str(Theta_0),' rad'])

%% scattering angle
[kx,ky]=meshgrid([-N/2:N/2-1]);
kr=sqrt( kx.^2 + ky.^2 )/(N*px);

%% Fourier transform functions
FFT=@(x) fftshift(fft2(fftshift(x)));
IFFT=@(x) ifftshift(ifft2(ifftshift(x)));

%% Thick diffuseur generation
delta=ThickDiffuser(L, ls, g, lambda_1,  N, N_diff);

%% Propagation of lambda_1 & lambda_2 in the Thick Diffuser
%fields propagation
E_s_old=ones(N,N);%plane wave illumination
E_s_old_2=E_s_old;
for n=1:N_diff
    %   lambda 1
    E_s_new=E_s_old.*exp(i*k_1*delta(:,:,n));
    E_s_old=fresnel(E_s_new,lambda_1,d,px);
    
    %   lambda 2
    E_s_new_2=E_s_old_2.*exp(i*k_2*delta(:,:,n));
    E_s_old_2=fresnel(E_s_new_2,lambda_2,d,px);
end

E_inf=FFT(E_s_new);
I_inf=abs(E_inf).^2;
g_1=sum(I_inf(:).*cos(kr(:)*lambda_1))/sum(I_inf(:));
ang=sqrt((1-g_1));
disp(['Output scattering angle (computed @ lambda_1) : ',num2str(ang),' rad'])

E_inf_2=FFT(E_s_new_2);
I_inf_2=abs(E_inf_2).^2;
g_2=sum(I_inf_2(:).*cos(kr(:)*lambda_2))/sum(I_inf_2(:));
ang2=sqrt((1-g_2));
disp(['Output scattering angle (computed @ lambda_2) : ',num2str(ang2),' rad'])

z=L*(-Nz/2:Nz/2)*2/Nz;
for j=1:length(z)
     E_1=fresnel(E_s_new,lambda_1,z(j),px);
     E_2=fresnel(E_s_new_2,lambda_2,z(j),px);
     
     I_1=abs(E_1).^2;
     I_1=I_1-mean(I_1(:));%substracting the mean value
     En1=sum(I_1(:).^2);
     I_2=abs(E_2).^2;
     I_2=I_2-mean(I_2(:));%substracting the mean value
     En2=sum(I_2(:).^2);
     corr_product=real(IFFT( conj(FFT(I_1)) .* FFT(I_2) ));
     cross_corr(j)=corr_product(N/2+1,N/2+1)/sqrt(En1*En2);%normalization
end

% Analytical model from Zhu et al. Optica 7(4), 338 (2020)
k0=(k_1/k_2)*(k_1-k_2); %spectral detuning
C=1./( 1 + (k0*L*Theta_0^2)^2/18 + (Theta_0^2*k_1/k_2)^2 * ( k_1*(z+L/3)-k_2*(z+L/3) ).^2 ); %cross-correlation product

%% Figure
figure
rectangle('Position',[-L,0,L,1],'FaceColor',[1 1 0.75])%,'EdgeColor','k','LineWidth',0.5)
hold on
plot(z,C,'g','LineWidth',2)
plot(z,cross_corr,'+','LineWidth',2)
plot(L *[-1,-1],[0,1],'k','LineWidth',2)
plot(-L/3*[1,1],[0,1],'r','LineWidth',2)
plot(0*[1,1],[0,1],'k','LineWidth',2)
xlabel('axial coordinate (\mu m)')
title('spectro-axial cross-correlation product')
yticks([0 0.2 0.4 0.6 0.8 1])
xlim([-L-100 L])
xticks([-L 0 L])
ylim([0 1])
yticks([0 0.2 0.4 0.6 0.8 1])
legend('Analytical model','Numerical','Diffuser boundaries','-L/3','location','southeast') 
