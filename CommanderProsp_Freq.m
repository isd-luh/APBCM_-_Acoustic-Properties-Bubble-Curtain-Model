% APBCM - Acoustic Properties Bubble Curtain Model
% Copyright C 2026 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

%% Function to determine the effective sound velocity based Commander and Prosperetti (1989)
% There are some mistakes in the original publication from Commander and
% Prosperetti (1989). The final model equations has been taken from Kargl
% (2002)
function [c_mProspReal,c_mProspImag, Phi_vec] = CommanderProsp_Freq(Frequenz,Blasenradius_vec,Anzahl_Blasen_vec,Druck_stat,Temp,c_l,rho_l,sigma,gamma) 

%% Parameter
freq    = Frequenz;                     % Frequency  
p_stat  = Druck_stat;                   % static and equilibrium pressure in water
a_vec   = Blasenradius_vec';            % Bubble radius
N_a     = Anzahl_Blasen_vec';           % Number of bubbles per radius in unit volume
Temp    = Temp;                         % Temperature
c_l     = c_l;                          % Sound velocity water [ms^-1]
rho_l   = rho_l;                        % Density water [kgm^-3]
sigma   = sigma;                        % Surface tension
gamma   = gamma;                        % Adiabatic coefficient

%% Material values
mu_Stutz    = [1792*10^(-6);890*10^(-6)];           % Supporting points dynamic viscosity water from VDI Waermeatlas (2026)        
Temp_Stutz  = [0 ; 25];                             % Supporting points temperature  
mu          = interp1(Temp_Stutz,mu_Stutz,Temp);    

D_Stutz         = [1.89e-5 ; 0.378e-5];             % Supporting points thermal conductivity of air at 0 degree Celsius from VDI Waermeatlas (2026)
p_stat_Stutz    = [100000 ; 500000];                % Supporting points static pressure 
D               = interp1(p_stat_Stutz,D_Stutz,p_stat,'linear','extrap');

%% Model Commander and Prosperetti (1989)
a_vec   = reshape(a_vec,1,[]);
freq    = reshape(freq,[],1);
N_a     = reshape(N_a,[],1); 

[A,FREQ] = meshgrid(a_vec,freq);

omega   = 2*pi*FREQ;    
k_l     = omega/c_l;

p_0_vec     = p_stat + 2*sigma ./ A;                     
chi_vec     = D ./ (omega .* A.^2);                     
Phi_vec     = 3*gamma ./ (1 - 1i*3*(gamma-1) .* chi_vec .* ( (1i./chi_vec).^0.5 .* coth((1i./chi_vec).^0.5)-1 ));     
omega_0_vec = (1./(rho_l*A.^2).* ( p_0_vec.*real(Phi_vec)-(2*sigma)./A ) ).^0.5;                                        

b_viskos    = 4*mu*omega ./ (rho_l.*A.^2);
b_thermal   = p_0_vec .* imag(Phi_vec) ./ (rho_l .* A .^ 2);
b_rad       = omega.^2 .* k_l .* A;

b_vec       = b_viskos + b_thermal + b_rad;             

G = sum( (A ./ (omega_0_vec.^2 - omega.^2 +1i .* b_vec) ) * N_a, 2);             

k_mProsp = (k_l(:,1).^2 + 4*pi.*omega(:,1).^2.*G).^0.5;            
c_mProsp = omega(:,1) ./ k_mProsp;                           

c_mProspReal = real(c_mProsp);
c_mProspImag = imag(c_mProsp);    

end

% Reference
% Commander, K. W., and Prosperetti, A. (1989). "Linear pressure waves in
% bubbly liquids: Comparison between theory and experiments," J. Acoust. Soc. Am. 85.
% Kargl, S. G. (2002). "Effective medium approach to linear acoustics in
% bubbly liquids". en. In: The Journal of the Acoustical Society of America 111.1, p. 168. doi: 10.1121/1.1427356.
% VDI-Waermeatlas 2006
