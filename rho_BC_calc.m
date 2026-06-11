% APBCM - Acoustic Properties Bubble Curtain Model
% Copyright C 2026 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

%% Function to determine the local density in the bubble curtain

function rho_BC = rho_BC_calc(rho_w,air_frac_zx,rho_a_zx)

%% Parameter
rho_w       = rho_w;        % Density water
rho_a_zx    = rho_a_zx;     % Density air
air_frac_zx = air_frac_zx;  % Local air fraction   

%% Berechnung
rho_BC = rho_w .* (1-air_frac_zx) + rho_a_zx .* air_frac_zx;   % Blasenschleierdichte

end