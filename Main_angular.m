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
lambda=0.8; %wavelength µm
k=2*pi/lambda; %wavenumber
ls=100;%scattering mean free path (µm)
L=500;%diffuser thickness (µm)
g=0.985;%0.975;% %Anisotropy factor

Dteta=0.01;% incident angle (rad) (N_teta+1)*px/(2*z+L)

%numerical parameters
N=1024; %matrix dimension (square)
N_theta=32; %number of angular computed values
N_diff=30; %number of diffuser layers
px=0.5; %pixel size (µm)% WARNING: N px/L must be larger than Theta_0!!!

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
delta=ThickDiffuser(L, ls, g, lambda, px, N, N_diff);

%% reference field computation (normal incidence)
E_s_old_r=1;
for n=1:N_diff
        E_s_new_r=E_s_old_r.*exp(i*k*delta(:,:,n));
        E_s_old_r=fresnel(E_s_new_r,lambda,d,px);
end

E_r=fresnel(E_s_new_r,lambda,-L/2,px);
I_r=abs(E_r).^2;
I_r=I_r-mean(I_r(:));
En_r=sum(I_r(:).^2);

%% angular scanning of the beam
theta=Dteta*[0:N_theta]/N_theta; %theta values (rad)

h=waitbar(0,'please wait...');
for m=1:N_theta+1
    waitbar(m/(N_theta+1));
    E_s_old=exp(i*k*px*kx*theta(m));
    for j=1:N_diff
        E_s_new=E_s_old.*exp(i*k*delta(:,:,j));
        E_s_old=fresnel(E_s_new,lambda,d,px);
    end

    %field at z=-L/2
    E_tilt=fresnel(E_s_new,lambda,-L/2,px);
    I_tilt=abs(E_tilt).^2;
    I_tilt=I_tilt-mean(I_tilt(:));
    En_tilt=sum(I_tilt(:).^2);
    
    %correlation product
    corr_product=real(IFFT( conj(FFT(I_r)) .* FFT(I_tilt) ))/sqrt(En_r*En_tilt);
    
    figure(100)
    imagesc(corr_product)
    axis('equal')
    caxis([0 1])
    zoom(20)

    cross_corr(m,1)=corr_product(N/2+1,N/2+1);%u;
end
close(100)
close(h)

%% comparaison with model  

Cmax=exp(-theta.^2.*(k*L*Theta_0).^2/12); %cross-correlation product: model from optica Zhu et al.(2020)


%% 
figure
plot(theta,Cmax,'LineWidth',2)
hold on
plot(theta,cross_corr,'b+','LineWidth',2)
ax2 = gca;
ax2.FontSize=14;
xlabel('\theta (rad)','FontSize',14,'FontWeight','bold');
ylabel('cross-correlation product','FontSize',14,'FontWeight','bold');
lgd4=legend('model','computed');
lgd4.FontSize = 10;
ylim([0 1.1]);


