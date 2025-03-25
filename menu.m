clear, clc, close all
% MEF crazy beam by masthierryi - LAMEC  - UFPI
tic()

% OUTPUT __________________________________________________________________
BT = 1; % Beam theories: 1 = EBT; 2 = RBT; % 3 = SBT; 4 = TBT
modes = 1:4; % number of displayed frequencies and modes
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
    % 1000, 500e6, 0.3,   32e-3, 0.0250, 1,     27e-3, 6;
    % 1830, 100e6, 0.3,   35e-3, 0.0250, 1,    0.5e-3, 2;
    % 1360,  10e6, 0.3,   32e-3, 0.0250, 1,     12e-3, 3;
    % 1830, 100e6, 0.3,   35e-3, 0.0250, 1,    0.5e-3, 2;
% 
    7850, 205e9, 0.3,   0.01, 0.0250, 2,    LS2/10, 3; % simpa stepped
    7850, 205e9, 0.3,   0.01, 0.0375, 2,    LS1/10, 5;  
    7850, 205e9, 0.3,   0.01, 0.0250, 2,    LS2/10, 2;  
%
    % 2700,  69e9, 0.3,   0.01, 0.0250, 2,    LS2/10, 3 %1 simpa bimat
    % 7850, 205e9, 0.3,   0.01, 0.0250, 2,    LS1/10, 5;  
    % 2700,  69e9, 0.3,   0.01, 0.0250, 2,    LS2/10, 2;  
];   
    
data(inp).L_c   = sum(data(inp).segments(:,8)); % Soma do número de elementos por célula
data(inp).L = sum(data(inp).segments(:,7)) * data(inp).n_c;
% ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
data(inp).noty = 2; % input by length, then calculate nodes
data(inp).geo = 1; % lt = data(inp).n_el;

% multilayer data ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
lt = sum(data(inp).segments(:,end)) * data(inp).n_c;
% data(inp).n_c_l = 10;
% 
% % [radius, rho, E, nu, inital node, end node]
% data(inp).layer = [3e-3, 1830, 5000e6, 0.3, 1, 6;
%                                 0, 0, 0, 0, 7, 8;
%                    3e-3, 1200, 175e6, 0.3, 9, 11;
%                               0, 0, 0, 0, 12, 13;
% ]; 
% ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨

% boundary conditions ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
data(inp).BC = [ 1      0 1 ;  % at node 1, BC is 
                 lt+1   0 1 ]; % 1 1 for coupling
% ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% ----------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------

% % Input 2 ------- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------------------
% inp = 2; 
% % parameters  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ 
% data(inp).n_c = 1;
% data(inp).segments = [
%   1000, 500e6, 0.3,   32e-3, 0.0250, 1,     27e-3, 6;
% ];   
% 
% data(inp).L_c   = sum(data(inp).segments(:,8)); 
% data(inp).L = sum(data(inp).segments(:,7)) * data(inp).n_c;
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% data(inp).noty = 2; data(inp).geo = 2;
% 
% % multilayer data ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% lt = sum(data(inp).segments(:,end)) * data(inp).n_c;
% data(inp).n_c_l = 1;
% 
% data(inp).layer = [3e-3, 1830, 5000e6, 0.3, 1, 6;
% ]; 
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% 
% % boundary conditions ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% data(inp).BC = [ 1      1 1  ;  % 1 1 for coupling
%                  lt+1   0 0 ];  
% % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
% % ----------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------

% _________________________________________________________________________
FEA = calculations(data,BT,modes); 

% basic analysis
calculations.Draw(FEA); %2.1s for 100 el
ShapeModes(FEA,BT,modes); % Mode shape for the chosen theory
% TheoriesModeShape(FEA,modes); % Theories shape modes comparing

% parameter analysis
type = 3;
% type = the desired relationship to plot
% type = 1: Euler-Bernoulli and Rayleight
% type = 2: Euler-Bernoulli and Shear
% type = 3: Euler-Bernoulli and Timoshenko

% Rho_natFreq(FEA,modes,data,type);
% E_natFreq(FEA,modes,data,type);
% slendernessR_natFreq(FEA,modes,data,type); 
Freq_spectra(FEA,BT,22);
Periodic_Saturation(FEA,18:22,data); 

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
