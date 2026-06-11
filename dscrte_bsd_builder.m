% APBCM - Acoustic Properties Bubble Curtain Model
% Copyright C 2026 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

%% Function for generating the discrete probability distribution df from a continous probability distribution f 
function  df = dscrte_bsd_builder(a, a_bnds, f_a, dstrbn_type, a_char_1, a_char_2) 

%% Parameter
a_bnds      = a_bnds;           % Bubble radius class boundaries
a           = a;                % Bubble radius class center
f_a         = f_a;              % Given distribution 
dstrbn_type = dstrbn_type;      % Type of bubble size distribution
a_char_1    = a_char_1;         % Characteristic bubble size 1
a_char_2    = a_char_2;         % Characteristic bubble size 2

%% Generate distribution
switch(dstrbn_type)
    case 'Gauss'      
        a_mu        = a_char_1;                                                        
        a_sigma     = a_char_2;                                                     
        f_a         = 1/((2*pi())^0.5*a_sigma)*exp(-0.5*((a-a_mu)./a_sigma).^2);
        f_bnds      = interp1(a,f_a,a_bnds,'pchip',0);
        flag_unfrm  = false;
        
    case 'modeled'             
        f_bnds      = interp1(a,f_a,a_bnds,'pchip',0);
        flag_unfrm  = false;
        
    case 'lognormal'
        a_mu        = a_char_1;                                                        
        a_sigma     = a_char_2;
        f_a         = lognpdf(a,a_mu,a_sigma);
        f_bnds      = interp1(a,f_a,a_bnds,'pchip',0);
        flag_unfrm  = false;
        
    case 'uniform'
        f_bnds      = zeros(1,length(a_bnds));
        [~,logical_cntr_char] = min(abs(a-a_char_1));
        flag_unfrm  = true;

    otherwise
        error('Unkown bubble size distribution')
end

%% Discrete probability distribution  (df)
df = zeros(1,length(a));
da = a_bnds(2:end)-a_bnds(1:end-1);

if flag_unfrm == true
    df(logical_cntr_char) = 1;
else
    f_bnds  = f_bnds./( trapz(a_bnds,f_bnds) );                       
    df      = ( f_bnds(1:end-1) + f_bnds(2:end) )/2 .* da;                
end

end

