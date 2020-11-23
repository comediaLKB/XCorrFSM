% Generation of a smooth surfacic diffuser
% 
% © Marc Guillon et al. Opt. Express 25, 12640-12652 (2017)
% https://doi.org/10.1364/OE.25.012640
%
% function delta=ThinDiffuser(e,w,N)
%
% INPUTs:
% e: mean squared roughness (any unit)
% w: correlation-width of the phase mask (in pixels)
% N: matrix dimension (square)
% 
% OUTPUT:
% delta: optical path difference (same unit as input "e")

function delta=ThinDiffuser(e,w,N)

FFT=@(x) fftshift(fft2(fftshift(x)));
IFFT=@(x) ifftshift(ifft2(ifftshift(x)));

[x,y]=meshgrid([-N/2:N/2-1]);
r2=x.^2+y.^2;

G=w*e*2*sqrt(3/pi);
g=G/(2*pi*w^2)*exp(-0.5*r2/w^2);


delta=2*pi*(rand(N,N)-0.5);

delta= real( IFFT( FFT(g).* FFT( delta )) ); 

return

