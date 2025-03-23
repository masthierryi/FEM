clear, clc, close all
% MEF crazy beam by masthierryi - LAMEC  - UFPI
tic()

% OUTPUT __________________________________________________________________
BT = 1; % Beam theories: 1 = EBT; 2 = RBT; % 3 = SBT; 4 = TBT
modes = 3; % number of displayed frequencies and modes
% _________________________________________________________________________

% INPUT ___________________________________________________________________

% Input type 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
inp = 1; % beam index, for each beam, add one on its index
% parameters  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ 

data(inp).n_c = 10; % número de células unitárias

alpha_artigo = 1-(4/8);
LS2 = (alpha_artigo/2);
LS1 = 1-2*LS2;

data(inp).segments = [
  % [rho, E, nu,        d1, d2, form,       section_length, n_elements]

    % 7850, 205e9, 0.3,   0.01, 0.0250, 2,    LS2, 3; % simpa stepped
    % 7850, 205e9, 0.3,   0.01, 0.0375, 2,    LS1, 5;  
    % 7850, 205e9, 0.3,   0.01, 0.0250, 2,    LS2, 2;  

    2700,  69e9, 0.3,   0.01, 0.0250, 2,    LS2/10, 3 %1 simpa bimat
    7850, 205e9, 0.3,   0.01, 0.0250, 2,    LS1/10, 5;  
    2700,  69e9, 0.3,   0.01, 0.0250, 2,    LS2/10, 2;  
];   
    
data(inp).L_c   = sum(data(inp).segments(:,8)); % Soma do número de elementos por célula
data(inp).L = sum(data(inp).segments(:,7)) * data(inp).n_c;
% ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
data(inp).noty = 2; % input by length, then calculate nodes
data(inp).geo = 1; % lt = data(inp).n_el;

% multilayer data ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% [radius, rho, E, nu, inital node, end node] 
% one vector (line .: ";") for each new layer
lt = sum(data(inp).segments(:,end)) * data(inp).n_c;
data(inp).layer = {[0.01, 2710.3, 68.73e9, 0.3, 1,lt]}; 
% ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨

% boundary conditions ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
data(inp).BC = [ 1      0 1 ;  % at node 1, BC is 
                 lt+1   0 1 ]; % 1 1 for coupling
% ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% -----------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------------

% % Input type 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % -------------------------------------------------------------------------
% inp = 1; % beam index, for each beam, add one on its index
% % FEM parameters  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨  ¨
% data(inp).n_el = 100; % number of elements
% lt = data(inp).n_el; % number of the last elementa
% 
% % beam parameters ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % [value, node*]; *last node that has this value. one line for each change
% data(inp).L   = 0.32122; %[m] Length
% 
% fa = 1; % o valor do interno é x vezes o externo
% data(inp).rho = {[fa*2710.3, 1, lt]}; %[kg/m^3] Specifc mass
% 
% data(inp).E   = {[fa*68.73e9, 1, lt]}; %[Pa] Youngs modulus
% 
% data(inp).nu  = {[0.3, 1, lt]}; %[~] Poissons ratio
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% 
% % cross-section data  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % [d1**, d2**, form*, inital node, end node] 
% % *1 = circular; 2 = square; 3 = pipe; 4 = box; 20 = custom values;
% %
% % **form = 1 -> d1 is the radius
% % **form = 2 -> d1 and d2 are horizontal and vertical measures, respectively
% % **form = 3 -> d1 is the external radius, and d2, the internal radius
% % **form = 4 -> d1 is the external heigth, and d2, the internal heigth
% %
% % p.s.: this type of data.d only works for cross section forms with
% % two measures, for form = 1, the second measure is ignored.
% %
% % Use geo = 1 (common) for standard cross-sections: circle, retangle, pipe 
% % and box). on future implementations it may have a geometry method for tapered.
% % geo = 2 (multilayer) for a multilayer cross section, you must to use
% % square or circle dimensions, the inertia and layer for retangular box is
% % not present
% %
% % {[],[],[]...[]} and adjusting the nodes for stepped
% data(inp).d = {[0.010002, 0, 1, 1, lt]}; 
% 
% % geo: 1 = common geometry; 2 = multilayer                                 
% % not working for timoshenko
% data(inp).geo = 2; data(inp).noty = 1; 
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% 
% % multilayer data ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % [radius, rho, E, nu, inital node, end node] 
% % one line for each new layer
% data(inp).layer = {[0.01, 2710.3, 68.73e9, 0.3, 1,lt]}; 
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% 
% % boundary conditions ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % [node where it is applied, sheer stress parametrer, momentum parametrer];
% % 1 = free; 0 = restrained %[N/m, N*m/rad]
% lno = data(inp).n_el+1;
% data(inp).BC = [ 1      1 1 ;  % at node 1, BC is 
%                  lt+1   1 1 ]; % 1 1 for coupling
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % -------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------

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
