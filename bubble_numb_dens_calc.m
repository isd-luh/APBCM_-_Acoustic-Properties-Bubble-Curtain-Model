% APBCM - Acoustic Properties Bubble Curtain Model
% Copyright C 2026 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

%% Function to determine the local bubble number density
% number of bubbles in size class

function n_bub_vec = bubble_numb_dens_calc(a,df,air_frac)

%% Parameter
a           = a;            % Class center of bubble radius
df          = df;           % Discrete probability of class
air_frac    = air_frac;     % Air fraction

%% Determine bubble number density
v_a         = 4/3*pi*a.^3;              % Creates array with bubble volume      
v_mean      = dot(df,v_a);              % Determine mean bubble volume
n_bub_ges   = air_frac / v_mean ;       % Determine absoulte number of all bubbles 
n_bub_vec   = df * n_bub_ges ;          % bubble number density 

end