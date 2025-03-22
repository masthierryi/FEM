%% printando alguns paramentros para comparar codigos parecidos 
%{
% fprintf("\nU disp| %+e| %+e| %+e| %+e| \n",FEA.U_disp(3,3),FEA.U_disp(12,11),FEA.U_disp(5,18),FEA.U_disp(19,36));
% fprintf("U rot_| %+e| %+e| %+e| %+e| \n",FEA.U_rot(3,3),FEA.U_rot(12,11),FEA.U_rot(5,18),FEA.U_rot(19,36));
% fprintf("K_____| %+e| %+e| %+e| %+e| \n",FEA.matrices.CK(3,3),FEA.matrices.CK(25,23),FEA.matrices.CK(32,29),FEA.matrices.CK(32,30));
% fprintf("M_____| %+e| %+e| %+e| %+e| \n",FEA.matrices.CM(3,3),FEA.matrices.CM(25,23),FEA.matrices.CM(32,29),FEA.matrices.CM(32,30));
% fprintf("freq__| %+e| %+e| %+e| %+e| \n",FEA.result.natfreq(1),FEA.result.natfreq(2),FEA.result.natfreq(3),FEA.result.natfreq(5));
% fprintf("beam__| A = %+g| I = %+g| k = %+g| G = %+g| \n", FEA.Beam.A(1),FEA.Beam.I(1),FEA.Beam.k(1),FEA.Beam.G(1));
% fprintf("beam__| r_e = %+g| s_e = %+g| S = %+g| gama = %+g|",FEA.Beam.r_e(1),FEA.Beam.s_e(1),FEA.S,FEA.gama);
%}
%% TRUSS  em desenvolvimento, não adequado para a versão 2.0.0+
%{
% FEM parameters
% [element, initial node, final node , dof type*]
%*element type: 1- bar, 2- beam
data.t_el = [1 1 3 2;
             2 2 3 2];

% node coordinates
data.n_c = [0, 50; % node 1                 
            50, 50; % node 2
            100, 0];% node 3

% beam parameters
% for each node
data.rho = [7850];
data.E = [210e9];
data.poisson = [0.3];

% cross-section data
data.CS_form = [2]; % circle
data.CS_dim = [0.125,0.125]; 

% boundary conditions
% [node where it is applied, sheer stress parametrer, momentum parametrer];
% 1 = free; 0 = restrained
data.BC = [ 1 1 0 ;  % at node 1, BC is 
            2 1 0 ]; % x = L 
%}
%% input manual do exemplo 3 do artigo do zé 2019 modelagem de coluna humana lamec
%{
data(inp).rho = {[1000, 1  , 270],  [1830, 271, 275],  [1360, 276, 395],  [1830, 396, 400],...
                 [1000, 401, 670],  [1830, 671, 675],  [1360, 676, 795],  [1830, 796, 800],...
                 [1000, 801, 1070], [1830, 1071, 1075],[1360, 1076, 1195],[1830, 1196, 1200],...
                 [1000, 1201, 1470],[1830, 1471, 1475],[1360, 1476, 1595],[1830, 1596, 1600],...
                 [1000, 1601, 1870]}; %[kg/m^3] Specifc mass

data(inp).E   = {[500e6, 1  , 270],  [100e6, 271, 275],  [10e6, 276, 395],  [100e6, 396, 400],...
                 [500e6, 401, 670],  [100e6, 671, 675],  [10e6, 676, 795],  [100e6, 796, 800],...
                 [500e6, 801, 1070], [100e6, 1071, 1075],[10e6, 1076, 1195],[100e6, 1196, 1200],...
                 [500e6, 1201, 1470],[100e6, 1471, 1475],[10e6, 1476, 1595],[100e6, 1596, 1600],...
                 [500e6, 1601, 1870]}; %[Pa] Youngs modulus

data(inp).nu  = {[0.3, 1, lt]}; %[~] Poissons ratio

% cross-section data
data(inp).d =  {[0.032, 0, 1, 1  , 270],  [0.035, 0, 1, 271, 275],  [0.032, 0, 1, 276, 395],  [0.035, 0, 1, 396, 400],...
                [0.032, 0, 1, 401, 670],  [0.035, 0, 1, 671, 675],  [0.032, 0, 1, 676, 795],  [0.035, 0, 1, 796, 800],...
                [0.032, 0, 1, 801, 1070], [0.035, 0, 1, 1071, 1075],[0.032, 0, 1, 1076, 1195],[0.035, 0, 1, 1196, 1200],...
                [0.032, 0, 1, 1201, 1470],[0.035, 0, 1, 1471, 1475],[0.032, 0, 1, 1476, 1595],[0.035, 0, 1, 1596, 1600],...
                [0.032, 0, 1, 1601, 1870]}; 

% geo: 1 = common geometry; 2 = multilayer                             
% not working for timoshenko
data(inp).geo = 2; 

% multilayer data
% [radius, rho, E, nu, inital node, end node] 
% one line for each new layer
data(inp).layer = {[0.003, 1830, 5000e6, 0.3, 1  , 270];  [0.003, 1200, 175e6, 0.3, 276, 395];      ...
                   [0.003, 1830, 5000e6, 0.3, 401, 670];  [0.003, 1200, 175e6, 0.3, 676, 795];      ...
                   [0.003, 1830, 5000e6, 0.3, 801, 1070]; [0.003, 1200, 175e6, 0.3, 1076, 1195];    ...
                   [0.003, 1830, 5000e6, 0.3, 1201, 1470];[0.003, 1200, 175e6, 0.3, 1476, 1595];    ...
                   [0.003, 1830, 5000e6, 0.3, 1601, 1870]}; 
% end 

% boundary conditions
lno = data(inp).n_el+1;
data(inp).BC = [ 1      0 0 ;  
                 lt+1   0 0 ];
%}