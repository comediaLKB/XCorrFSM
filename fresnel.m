%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Fresnel propagation in vacuum               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Marc Guillon (marc.guillon@u-paris.fr)
% 
% function B = fresnel(A,lambda,z,px)
%
% INPUTs:
% A: ( N x N matrix) wavefield at z=0
% lambda: wavelength
% z: axial shift
% px: pixel size
% 
% Nota Bene: lambda, z and px should be in the same unit.
% 
% OUTPUTs:
% B: ( N x N matrix) wavefield at z
%


function B=fresnel(A,lambda,z,px)

NA=1;

OS=lambda/(2*NA*px);
N=size(A,1);

[x,y]=meshgrid([-N/2:N/2-1]);
x=x/N;
y=y/N;
r2=x.^2+y.^2;

defocus=exp( 2*i*pi* z*NA^2/(2*lambda) *4* r2 * OS^2);

B=ifftshift(ifft2( ifftshift( defocus ) .* fft2(fftshift(A)) ));

return






