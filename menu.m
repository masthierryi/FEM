
lear, clc, close all
% MEF crazy beam by masthierryi
% # version 3.2: 
% -[Geometry.common][Material]common geometry and material parameters for each 
% element using vectorization and a for loop;
% -[menu] support beams with different mesh, lengths...; 
% -[Coupling]{Coupled_mesh} couples the matrices and mesh
% -[ModeShapes] the mode shapes are plotted in function of the coupled mesh
% -{matrices}{result}{critical} response summarization with structs

tic()

% OUTPUT __________________________________________________________________
BT = 1; % Beam theories: 1 = EBT; 2 = RBT; % 3 = SBT; 4 = TBT
modes = 4; % number of displayed frequencies and modes
% _________________________________________________________________________

% INPUT ___________________________________________________________________

% Input 1
% -------------------------------------------------------------------------
inp = 1; % beam index, for each beam, add one on its index
% FEM parameters  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨  ¨
data(inp).n_el = 20; % number of elements
lt = data(inp).n_el; % number of the last elementa

% beam parameters ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% [value, node*]; *last node that has this value. one line for each change
data(inp).L = 1;  %[m] Length

data(inp).rho = {[7830, 1, lt]}; %[kg/m^3] Specifc mass

data(inp).E =   {[200e9, 1, lt]}; %[Pa] Youngs modulus

data(inp).nu =  {[0.3, 1, lt]}; %[~] Poissons ratio
% ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨

% cross-section data  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% [d1**, d2**, form*, inital node, end node] 
% *1 = circular; 2 = square; 3 = pipe; 4 = box; 20 = custom values;
%
% **form = 1 -> d1 is the radius
% **form = 2 -> d1 and d2 are horizontal and vertical measures, respectively
% **form = 3 -> d1 is the external radius, and d2, the internal radius
% **form = 4 -> d1 is the external heigth, and d2, the internal heigth
%
% p.s.: this type of data.d only works for cross section forms with
% two measures, for form = 1, the second measure is ignored.
%
% Use geo = 1 (common) for standard cross-sections: circle, retangle, pipe 
% and box). on future implementations it may have a geometry method for tapered.
% geo = 2 (multilayer) for a multilayer cross section, you must to use
% square or circle dimensions, the inertia and layer for retangular box is
% not present
%
% {[],[],[]...[]} and adjusting the nodes for stepped
data(inp).d = {[0.4, 0.39, 2, 1, lt/2],[0.2, 0.19, 2, lt/2+1, lt]}; %[0.4, 0.4, 1, 2, 1, lt]

% geo: 1 = common geometry; 2 = multilayer
data(inp).geo = 1; 
% ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨

% multilayer data ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% [radius, rho, E, nu, inital node, end node] 
% one line for each new layer
% data(inp).layer = {[0.04, 5000, 20e9, 0.3, 2,lt]}; % 0.03, 5000, 20e9, 0.3, 1,lt
% ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨

% boundary conditions ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% [node where it is applied, sheer stress parametrer, momentum parametrer];
% 1 = free; 0 = restrained %[N/m, N*m/rad]
lno = data(inp).n_el+1;
data(inp).BC = [ 1 1 1 ;  % at node 1, BC is 
               lno 0 1 ]; % 1 1 for coupling
% ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% -------------------------------------------------------------------------

% % Input 2
% % -------------------------------------------------------------------------
% inp = 2;
% FEM parameters  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% data(inp).n_el = 1;
% lt = data(inp).n_el; 
% 
% % beam parameters ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% data(inp).L = 1; 
% 
% data(inp).rho = [7850,lt];
% 
% data(inp).E = [210e9, lt];
% 
% data(inp).nu = [0.3, lt];
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% 
% % cross-section data  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ 
% data(inp).d = [0.5, 0.5, 1, lt];
% 
% data(inp).geo = 2;
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% 
% % multilayer data ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% data(inp).layer = [0.04, 5000, 20e9, 0.3, 1, lt];
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% 
% % boundary conditions ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % for coupling beams, the last node of actual beam must to be free 1 1, and 
% % the first node of the next beam must to be campled 0 0, but the code does
% % this automatically, so let it free 1 1
% lno = data(inp).n_el+1;
% data(inp).BC = [ 1 1 1 ; % at node 1, BC is 
%                lno 0 0 ]; % x = L 
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % -------------------------------------------------------------------------

% % TRUSS 
% % -------------------------------------------------------------------------
% % FEM parameters  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ 
% % [element, initial node, final node , dof type*]
% %*element type: 1- bar, 2- beam
% data.t_el = [1 1 3 2;
%              2 2 3 2];
% 
% % node coordinates
% data.n_c = [0, 50; % node 1                 
%             50, 50; % node 2
%             100, 0];% node 3
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% 
% % beam parameters ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % for each node
% data.rho = [7850];
% data.E = [210e9];
% data.poisson = [0.3];
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% 
% % cross-section data  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% data.CS_form = [2]; % circle
% data.CS_dim = [0.125,0.125]; 
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% 
% % boundary conditions ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % [node where it is applied, sheer stress parametrer, momentum parametrer];
% % 1 = free; 0 = restrained
% data.BC = [ 1 1 0 ;  % at node 1, BC is 
%             2 1 0 ]; % x = L 
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % -------------------------------------------------------------------------
% _________________________________________________________________________

% _________________________________________________________________________
% For each beam provided in the pages of the data (array inside the array), 
% its matrix are calculated individually. Then, these matrices are assembled 
% into a global matrix and the eigenvalue problem is solved for the entire system.
FEA = calculations(data,BT,modes); 

% basic analysis
% FEA = calculations.Draw(FEA); %2.1s for 100 el
FEA = ShapeModes(FEA,BT); % Mode shape for the chosen theory
% FEA = TheoriesModeShape(FEA); % Theories shape modes comparing

% parameter analysis
type = 3;
% type = the desired relationship to plot
% type = 1: Euler-Bernoulli and Rayleight
% type = 2: Euler-Bernoulli and Shear
% type = 3: Euler-Bernoulli and Timoshenko

% FEA = Rho_natFreq(FEA,modes,data,type);
% FEA = E_natFreq(FEA,modes,data,type);
% FEA = slendernessR_natFreq(FEA,modes,data,type); 

% _________________________________________________________________________

toc()

fprintf("\nU disp| %+e| %+e| %+e| %+e| \n",FEA.U_disp(3,3),FEA.U_disp(12,11),FEA.U_disp(5,18),FEA.U_disp(19,36));
fprintf("U rot_| %+e| %+e| %+e| %+e| \n",FEA.U_rot(3,3),FEA.U_rot(12,11),FEA.U_rot(5,18),FEA.U_rot(19,36));
fprintf("K_____| %+e| %+e| %+e| %+e| \n",FEA.matrices.CK(3,3),FEA.matrices.CK(25,23),FEA.matrices.CK(32,29),FEA.matrices.CK(32,30));
fprintf("M_____| %+e| %+e| %+e| %+e| \n",FEA.matrices.CM(3,3),FEA.matrices.CM(25,23),FEA.matrices.CM(32,29),FEA.matrices.CM(32,30));
fprintf("freq__| %+e| %+e| %+e| %+e| \n",FEA.result.natfreq(1),FEA.result.natfreq(2),FEA.result.natfreq(3),FEA.result.natfreq(5));
fprintf("beam__| A = %+g| I = %+g| k = %+g| G = %+g| \n", FEA.Beam.A(1),FEA.Beam.I(1),FEA.Beam.k(1),FEA.Beam.G(1));
fprintf("beam__| r_e = %+g| s_e = %+g| S = %+g| gama = %+g|",FEA.Beam.r_e(1),FEA.Beam.s_e(1),FEA.S,FEA.gama);

% 
% time chart --------------------------------------------------------------
% only computting the frequencies: 0.07s
% [calculations.Draw]: % 1.2075174s for 100 elements
% [ShapeModes]: % 0.382775s
% [TheoriesModeShape]: % 0.5978928s
% [Rho_natFreq]: 22.3442085s for nmr = 320 repetitions with 2 theories
% [E_natFreq]: 22.472827s for nmr = 320 repetitions with 2 theories
% [slendernessR_natFreq]: 106.303478s for nmr = 1600 repetitions
%
% testing pc: xeon e5 2640 v4 2.2~2.5GHz CPU; 16GB RAM; 250GB SSD nvme 3rd gen;
% radeon gx580 sp
% -------------------------------------------------------------------------

% 
% improvements to be made: ------------------------------------------------
% * r, s, A and I parameters for freq calculation for higher
% theories: what is these parameters on beams that have no constant
% material and geomtric parameters? talvez seja realmente o primeiro
% parametro na viga, usado pra calcular a freq e tals, ja q a mudaça
% relacionada a não constancia desses parametros é refletida nas matrizes
% de massa e rigidez. mas pq a freq da diferente entre vigas normais e
% modificadas?
% * tirar o U_disp e U_rot da resposta.
% * validate Rho_natFreq, E_natFreq; review slendernessR_natFreq, is not
% touching the x_axis even with nmr = 1600: ask slepton.
% * o Draw deve ser relacionar o form tbm, pros hollow muda mto
% * implement Inertia moment and calculations for rectangular box cross
% sections, then will be possible to make rectangular beam, same for elipse
% 
% -------------------------------------------------------------------------