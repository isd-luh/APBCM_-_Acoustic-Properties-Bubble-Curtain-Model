% APBCM - Acoustic Properties Bubble Curtain Model
% Copyright C 2026 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

%% Function to create bubble classes 
% Each class is characterized by its center and boundaries
function [cntr_vec, bnds_vec] = class_create(Min,Max,n_class,startingPoint)

% Grenze --> Mitte
if strcmp(startingPoint,'Boundary') == true
    bnds_vec  = linspace( Min , Max , n_class+1 );
    cntr_vec   = (bnds_vec(2:end) + bnds_vec(1:end-1)) / 2;                         

% Mitte --> Grenze
elseif strcmp(startingPoint,'Center') == true
    cntr_vec          = linspace(Min,Max,n_class);
    bnds_vec          = cntr_vec - (cntr_vec(3)-cntr_vec(2))/2;
    bnds_vec(end+1)   = cntr_vec(end) + (cntr_vec(3)-cntr_vec(2))/2;

else
    print('wrong startingPoint')

end