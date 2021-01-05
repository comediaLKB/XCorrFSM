% Generation of a thick forward scattering diffuser
% made of a stack of phase plates
% 
% This code is written by Payvand Arjmand & Marc Guillon
% Contact info: 
% payvand.arjmand@u-paris.fr
% marc.guillon@u-paris.fr
% 
% function [e,w]=ThickDiffuser(L, ls, g, lambda, N_diff)
% 
% INPUTs:
% L: Diffuser Thickness
% ls: scattering mean free path
% g: anisotropy factor
% lambda: wavelength
% N: matrix size (square)
% N_diff: Number of phase plates
%
% OUTPUTs:
% delta: (NxNxN_diff) mean squared roughness
% Theta_0: (rad) output scattering angle

function [delta, Theta_0]=ThickDiffuser(L, ls, g, lambda,px, N, N_diff)

    if nargin==5
        N_diff=10;
    end
    
    k=2*pi/lambda;
    Theta_0=sqrt(L*(1-g)/ls);
    
    if Theta_0>pi/4
       warning('Output scattering angle larger than pi/4:Limit of the model') 
    end
    
    e = sqrt( L / (N_diff*ls*k^2) );
    w = sqrt(N_diff/2)*e / (px*Theta_0);
    
    for j=1:N_diff
        delta(:,:,j)=ThinDiffuser(e,w,N);
    end

end