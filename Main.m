% APBCM - Acoustic Properties Bubble Curtain Model
% Copyright C 2026 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

%% Script to determine the effective sound velocity in the bubble curtain

%% Parameter
freq                    = 100:100:500;
T                       = 40;                               % water depth (m)
B                       = 8;                                % maximal width of the bubble curtain (m)
n_z                     = 100;                              % number of supporting points (-)
n_x                     = 100;                              % number of supporting points (-)
a_char_1                = -5.41;                            % characterstic bubble size 1; depends on the choice of dstrbtn type
a_char_2                = 0.34;                             % characterstic bubble size 2; depends on the choice of dstrbtn type
a_max                   = 0.02;                             % maximum bubble size (m)
q_atm                   = 0.002;                            % air volume flow at atmospheric pressure per 1m nozzle hose (m^2/s)
dstrbn_type             = 'modeled';                        % type of bubble size distribution ('Gauss', 'modeled', 'lognormal', 'uniform')
n_class                 = 500;                              % number of size classes (-)
d_Orifice               = 0.005;                            % nozzle diameter (m)
dy_Orifice              = 0.3; 	                            % nozzle distance (m)
koal_cond               = [1,1];                            % coalescence conditions: [1,1]: Coalescence and Break up , [1,0]: Coalescence and no Break up, [0,1]: Break up and no Coalescence 
rho_l                   = 1000;                             % water density (kg/m^3)
p_atm                   = 100000;                           % atmospheric pressure (Pa)
c_l                     = 1480;                             % water sound velocity (m/s)
temp                    = 12;                               % temperature water (degree Celsius)
sigma                   = 0.072;                            % surface tension (N/m)
gamma                   = 1.4;                              % adiabatic coefficient (-)
z_eval                  = 2;                                % height at which axial bubble flow is evaluated and the bubble size distribution is extracted (m)

%% Include models
addpath(genpath('ABFM - Axial Bubble Formation Model'));
addpath(genpath('PBPM - Plane Bubble Plume Model'));

%% Generate bubble classes
[a,a_bnds] = class_create(0,a_max,n_class,'Boundary');

%% Assessing the bubble size distribution
if strcmp(dstrbn_type, 'modeled')
    q_atm_Orifice   = q_atm * dy_Orifice;
    [ ~, a, ~, ~, ~, ~, ~, ~, ~, f_a ] = bsd_bimodal( a, T, q_atm_Orifice, p_atm, rho_l, sigma, temp, d_Orifice, koal_cond, z_eval );
else
    f_a = zeros(length(a),1);
end

%trapz(a,f_a) % Test, should be close to 1 for dstrbn_type = 'modeled'

%% Build discrete bubble size distribution
df = dscrte_bsd_builder(a,a_bnds,f_a,dstrbn_type,a_char_1,a_char_2);

%sum(df) % Test, should be 1

%% Create mesh
[X_Stutz,~,Z_Stutz] = bc_mesh_create( T,B,1,n_x,1,n_z,'2-dimensional' );

%% Create uniform temperature field
temp_zx = ones(size(X_Stutz)) .* temp;

%% Determine static pressure field
p_stat_zx = p_stat_calc(Z_Stutz,rho_l,p_atm);          

%% Determine air density field
rho_g_zx = rho_Luft_calc(temp_zx,p_stat_zx);

%% Determine air fraction field
disp('Air fraction')

bc_plane_coeff.alpha    = 0.16;
bc_plane_coeff.lambda   = 0.2;
bc_plane_coeff.u_rel    = 0.4;
bc_plane_coeff.ampli    = 1.0;
bc_plane_coeff.h_0      = 0.0;
type                    = 'plane';
[ ~, ~, ~, ~, air_frac_zx, ~ ] = bc_Bohne( T, q_atm, p_atm, rho_l, sigma, temp, bc_plane_coeff , dy_Orifice, d_Orifice, X_Stutz, Z_Stutz, type );

%% Determine density field
rho_BC_zx = rho_BC_calc(rho_l,air_frac_zx,rho_g_zx);

%% Determine bubble number density field
n_bub_zx = cell(size(Z_Stutz));
for i = 1 : size(Z_Stutz,1)                         % z 
    for j = 1 : size(Z_Stutz,2)                     % x
        n_bub_zx{i,j} = bubble_numb_dens_calc(a,df,air_frac_zx(i,j));           
    end
end  

%% Determine effective sound velocity field
disp('Effective sound velocity')
c_real_Stutz_x_f = cell( size(Z_Stutz) );
c_imag_Stutz_x_f = cell( size(Z_Stutz) );
 
for i = 1 : size(Z_Stutz,1)                 % z
    for j = 1 : size(Z_Stutz,2)             % x

        if Z_Stutz(i,1) < -T
            c_real_Stutz_x_f{i,j} = ones( size(freq))*c_l;
            c_imag_Stutz_x_f{i,j} = zeros(size(freq));
        else
            [c_real_Stutz_x_f{i,j},c_imag_Stutz_x_f{i,j},~] = CommanderProsp_Freq(freq,a,n_bub_zx{i,j},p_stat_zx(i,j),temp_zx(i,j),c_l,rho_l,sigma,gamma);
        end

    end
end

%% Transform coordinates 
for i = 1 : size(Z_Stutz,1)                 % z
    for j = 1 : size(Z_Stutz,2)             % x
        for k = 1:length(c_real_Stutz_x_f{1,1})
            c_real_Stutz{k}(i,j) = c_real_Stutz_x_f{i,j}(k);
            c_imag_Stutz{k}(i,j) = c_imag_Stutz_x_f{i,j}(k);
        end
    end
end

%% Output
figure; contourf(X_Stutz,Z_Stutz,air_frac_zx); title('Air fraction')
figure; plot(a,df);  title('Bubble size distribution')
figure; contourf(X_Stutz,Z_Stutz, c_real_Stutz{end}); title('Effective sound velocity (Real)')
figure; contourf(X_Stutz,Z_Stutz, c_imag_Stutz{end}); title('Effective sound velocity (Imag)')


%% References
%Tobias Bohne, Tanja Griessmann , Raimund Rolfes, 2019. Modeling the noise mitigation of a bubble curtain. The Journal of the Acoustical Society of America 146 (4). doi: 10.1121/1.5126698
%Tobias Bohne, Tanja Griessmann, Raimund Rolfes, 2020. Development of and efficient buoyant jet integral model of a bubble plume coupled with a population dynamics model for bubble breakup and coalescence to predict the transmission loss of a bubble curtain. International Journal of Multiphase Flow 132, 103436. doi: 10.1016/j.ijmultiphaseflow.2020.103436