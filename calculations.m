classdef calculations
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties 
        data
        mesh
        Beam

        GM
        GK
        U_disp
        U_rot
        
        critical
        eigenvalues                 %weg
        freq                        %rg
        natfreq                     %weg
        natfreqHz                   %weg

    end                                                             methods
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POO nonsense (thinking aloud, m sorry) ================================== 
    function self = calculations(data,BT)
        self.data = data;

        self = FEM_setup(self);
        self = parameters(self,BT);
        self = Eigenproblem(self);
        self = EigenSolve(self,BT);
    end

%% FEM setup ==============================================================
    function self  = FEM_setup(self)
            
            self.mesh.t_nel = size(self.data.t_el,2); % number of truss elements
            % self.mesh.n_el = self.data.n_el; % number of elements
            self.mesh.n_el = self.mesh.t_nel; % number of elements ########
            self.mesh.n_nodes = self.data.n_el+1;  % number of nodes
            self.mesh.l_el = (self.data.L/self.mesh.n_el)/2; % lenght of an isolated el
    
            % Coordinates increasing n_el times by n_el or L/n_el (adimensional)
            self.mesh.coordinates = linspace(0,self.data.L,self.mesh.n_nodes);
    end

%% geometry ===============================================================
    function self = parameters(self,BT)
                                                  % no = self.mesh.n_nodes;
                                                no = size(self.data.n_c,2);

        for j = 1:self.mesh.n_el % calculation for each truss element = = =
            % truss element length ----------------------------------------
            self.data.L_e = (1/2)*sqrt((self.data.n_c(self.data.t_el(i,3),1)-...
            self.data.n_c(self.data.t_el(i,2),2))^2 + ...
            (self.data.n_c(self.data.t_el(i,3),1)- ...
            self.data.n_c(self.data.t_el(i,2),2))^2);

            % geometry ----------------------------------------------------
            if self.data.CS_form == 1 % Circular
                A = pi/4*(self.data.CS_dim(:,1).^2);
                I = pi/64*(self.data.CS_dim(:,1).^4);
                K = 6*(1+self.data.poisson)./(7+6*self.data.poisson);
    
            elseif self.data.CS_form == 2 % Rectangular
                A = self.data.CS_dim(:,1).*self.data.CS_dim(:,2);
                I = (self.data.CS_dim(:,1).*self.data.CS_dim(:,2).^3)/12;
                K = 10*(1+self.data.poisson)./(12+11.*self.data.poisson);
            end
                                                     self.Beam(1:no,1) = A;
                                                     self.Beam(1:no,2) = I;
                                                     self.Beam(1:no,3) = K;
                                                     
            % material ----------------------------------------------------
                                           self.Beam(1:no,4) = self.data.E;
                                         self.Beam(1:no,5) = self.data.rho;
                                                            
                                     self.Beam(1:no,6) = self.data.poisson;
            
                                               % Transversal Elasticity (G)
               self.Beam(1:no,7) = self.data.E/(2*(1 + self.Beam(1:no,5)));

            % r and s parameters ------------------------------------------
            if BT == 1 % Euler
                r = 0; 
                s = 0; 
    
            elseif BT == 2 % Rayleight
                r = sqrt(self.Beam(1:no,2)./(self.Beam(1:no,1).*self.mesh.l_el^2));
                s = 0;
                r_t = sqrt(self.Beam(1:no,2)./(self.Beam(1:no,1).*self.data.L^2));
                s_t = 0;
    
            elseif BT == 3 % Shear
                r = 0;
                s = sqrt((self.Beam(1:no,4)./(self.Beam(1:no,3).* ...
                    self.Beam(1:no,7))).*r^2);
                s_t = sqrt((self.Beam(1:no,2)./(self.Beam(1:no,1).* ...
                    self.data.L)));
                r_t = 0;
    
            elseif BT == 4 % Timoshenko
                % r = \sqrt{\frac{I}{A L^2}}
                % s = \sqrt{\frac{E I}{k G A L^2}}
                r = sqrt(self.Beam(1:no,2)./(self.Beam(1:no,1).*self.mesh.l_el^2));
                s = sqrt((self.Beam(1:no,4)./(self.Beam(1:no,3).* ...
                self.Beam(1:no,7))).*r.^2);
                r_t = sqrt(self.Beam(1:no,2)./(self.Beam(1:no,1).*self.data.L^2));
                s_t = sqrt((self.Beam(1:no,4)./(self.Beam(1:no,3).* ...
                self.Beam(1:no,7))).*r_t.^2);
            end
                                                     self.Beam(1:no,8) = r;
                                                     self.Beam(1:no,9) = s;
                                                    self.data.r_t = r_t(1);
                                                    self.data.s_t = s_t(1);
        % -----------------------------------------------------------------
        % Explanation: an i x j matrix is ​​created, where i is equal to n_el, 
        % and each column contains, in order, the data:
        %
        % 1- Cross-Section Area; 2-Moment of Inertia; 3- Shear Correction Factor, 
        % 4- Young's Modulus, 5- Density, 6- Poisson's Modulus, 
        % 7- Transversal Elasticity Modulus, 8- Rotatory Inertia Parameter
        % 9- Shear deformation parameter
        %
        % This configuration allows allocation of all the data needed for 
        % modeling of conical, stepped, graduated beams and beams with other 
        % peculiarities.
        end % end truss elements parameters calculation = = = = = = = = = =

        % 1; 2; 3; 4; 5  ; 6 ; 7; 8; 9;
        % A; I; k; E; rho; nu; G; r; s;      # self.Beam()
    end % end of function parameters

%% mesh ===================================================================
    function self = Eigenproblem(self)
        
        % Compute the node(x) and weight(w) for integration ---------------
        [x,w] = LGQ(self,4,-1,1);

        % Pre-allocation --------------------------------------------------
        % Stiffness global matrix
        self.GK = zeros(2*self.mesh.n_nodes,2*self.mesh.n_nodes);
        % Mass global matrix
        self.GM = self.GK;
        % Shape functions
        Nd = zeros(4); Ns = Nd; Nd_d = Nd; Ns_d = Nd;
        % -----------------------------------------------------------------
        % L_e = self.mesh.l_el;

        for i = 1:self.mesh.n_el % Loop for each element  = = = = = = = = =
            %r e s pra cadael, isso pro mesh novo, aq so vetor catandos

            % Shape functions==============================================
            s = self.Beam(i,9);
            fc = 1/(4*(3*s^2+1)); % Constant for the function
            for j = 1:length(x)

                % functions for deflection --------------------------------
                Nd(1,j) = (2*(3*s^2+1)-3*(2*s^2+1)*x(j)+1*x(j)^3);
                Nd(2,j) = L_e*((3*s^2+1)-1*x(j)-(3*s^2+1)*x(j)^2+1*x(j)^3);
                Nd(3,j) = (2*(3*s^2+1)+3*(2*s^2+1)*x(j)-1*x(j)^3);
                Nd(4,j) = L_e*(-(3*s^2+1)-1*x(j)+(3*s^2+1)*x(j)^2+1*x(j)^3);
                Nd(:,j) = Nd(:,j)*fc;
                % functions for bending slope -----------------------------
                Ns(1,j) = (3*x(j)^2-3)/L_e;
                Ns(2,j) = (-1-2*(3*s^2+1)*x(j)+6*s^2+3*x(j)^2);
                Ns(3,j) = (3-3*x(j)^2)/L_e;
                Ns(4,j) = (-1+2*(3*s^2+1)*x(j)+6*s^2+3*x(j)^2);
                Ns(:,j) = Ns(:,j)*fc;
                % first derivate of Nd ------------------------------------
                Nd_d(1,j) = (-3*(2*s^2+1)+3*x(j)^2);
                Nd_d(2,j) = L_e*(-1-2*(3*s^2+1)*x(j)+3*x(j)^2);
                Nd_d(3,j) = (3*(2*s^2+1)-3*x(j)^2);
                Nd_d(4,j) = L_e*(-1+2*(3*s^2+1)*x(j)+3*x(j)^2);
                Nd_d(:,j) = Nd_d(:,j)*fc;
                % first derivate of Ns ------------------------------------
                Ns_d(1,j) = (6*x(j))/L_e;
                Ns_d(2,j) = (-2*(3*s^2+1)+6*x(j));
                Ns_d(3,j) = (-6*x(j))/L_e;
                Ns_d(4,j) = (2*(3*s^2+1)+6*x(j));
                Ns_d(:,j) = Ns_d(:,j)*fc;

            end % =========================================================

            % Prealocate for the sum --------------------------------------
            M_tra = zeros(4); % Translation element matrix
            M_rot = zeros(4); % Rotation element matrix
            K_ben = zeros(4); % Bending element matrix
            K_she = zeros(4); % Shear element matrix
            %--------------------------------------------------------------
        
            % Matrix constants --------------------------------------------
            % constant for translation mass matrix (\rho A L_e)
            C_tra = self.Beam(i,5).*self.Beam(i,1)*L_e;
            % constant for rotation mass matrix (r^2 \rho A L_e^3)
            C_rot = self.Beam(i,8).^2*self.Beam(i,5).*self.Beam(i,1).*L_e^3;
            % constant for stiffness bending matrix (\frac{EI}{A})
            C_ben = self.Beam(i,4).*self.Beam(i,2)./L_e;
            % constant for stiffness shear matrix (s^2\frac{kappa*G*A)^2 
            %   L^3}{E*I}) (EI / L_e s^2)
            C_she = (s^2.*(self.Beam(i,3).*self.Beam(i,7).*self.Beam(i,1)).^2* ...
                     L_e^3)./(self.Beam(i,4).*self.Beam(i,2)); 
            % -------------------------------------------------------------

            % Computing the element matrices
            for j=1:4
                for k=1:4
                    for d=1:length(x)
                        M_tra(j,k) = Nd(j,d)*Nd(k,d)*w(d)+M_tra(j,k);
                        M_rot(j,k) = Ns(j,d)*Ns(k,d)*w(d)+M_rot(j,k);
                        K_ben(j,k) = Ns_d(j,d)*Ns_d(k,d)*w(d)+K_ben(j,k);
                        K_she(j,k) = (Nd_d(j,d)/L_e-Ns(j,d))*(Nd_d(k,d)/L_e-Ns(k,d))*w(d)+K_she(j,k);
                    end
                end
            end

            % Local to global coordinates ---------------------------------
            % Complete mass element matrix %tinha um set dps das constantes
            M_m = M_tra.*C_tra + M_rot.*C_rot;
            % Complete stiffness element matrix
            M_k = K_ben.*C_ben + K_she.*C_she;

            % angle of elemnt to x axis
            ang = atan((self.data.n_c(self.data.t_el(i,3),1)- ...
            self.data.n_c(self.data.t_el(i,2),2))/ ...
            (self.data.n_c(self.data.t_el(i,3),1)-...
            self.data.n_c(self.data.t_el(i,2),2)));

            % rotation matrix
            T_m = [  cos(ang) sin(ang) 0; % longitudinal
                    -sin(ang) cos(ang) 0; % transversal
                     0        0        1];% rotation

            % adjust the matrices for element type
            M_gaux = zeros(6,6); %inicia matriz com todos dof 2D
            if self.data.t_el(i,4) == 1
                T_m = T_m(logical([1 1 0]), :); % Remove rows 
                T_m = T_m(:, logical([1 1 0])); % Remove columns
                M_m = M_gaux([1 4],[1 4])+M_m;
                M_k = M_gaux([1 4],[1 4])+M_k;
            elseif self.data.t_el(i,4) == 2
                T_m = T_m(logical([0 1 1]), :); % Remove rows 
                T_m = T_m(:, logical([0 1 1])); % Remove columns
                M_m = M_gaux([2 3 5 6],[2 3 5 6])+M_m;
                M_k = M_gaux([2 3 5 6],[2 3 5 6])+M_k;
            end
            T_m = blkdiag(T_m,T_m); % duplicate for 2 nodes

            % transforming the element matrices to global coordinates
            M_m = T_m'.*M_m.*T_m;
            M_k = T_m'.*M_k.*T_m;
            % -------------------------------------------------------------

            % Local to global coordinates ---------------------------------

            % se é barra, é 4x4, pelo menos uma viga é 6x6, se os el
            % tem msm direção é 4x4, se são iguais, repete matriz
            % if self.beam é td igual, repete as matrizes

            % adjust the tranform matrix for element type
            % if size(unique(self.data.t_el(:,4)),2) ~= 1
            % end

            for n=1:self.mesh.n_el
                for k=((2*n)-1):(2*(n+1)) 
                    for j=((2*n)-1):(2*(n+1))
                        self.GM(k,j) = self.GM(k,j) + M_m((k-(2*(n-1))),(j-(2*(n-1))));
                        self.GK(k,j) = self.GK(k,j) + M_k((k-(2*(n-1))),(j-(2*(n-1))));
                    end
                end
            end 
            % end positionating on global ---------------------------------
            
        end % end element loop  = = = = = = = = = = = = = = = = = = = = = = 

        % Boundary Conditions ---------------------------------------------
        % auxiliary vector that transforms BC conditions into a column array
        BC_c = reshape([self.data.BC(:, 2), self.data.BC(:, 3)].', [], 1);
    
        % defines the position of the dofs where the conditions are applied
        BC_pos = sort([2*self.data.BC(:,1);2*self.data.BC(:,1)-1]);
    
        % pre-allocation of vector with all dof boundary condition
        BC_v(1:2*self.mesh.n_nodes,1) = 1;
    
        % adjust each dof condition 
        for l = 1:size(BC_c)
            BC_v(BC_pos(l)) = BC_c(l);
        end
    
        % use logical to recreate the global matrices with only the unconstrained dofs
        self.GM = self.GM(logical(BC_v), :); % Remove rows 
        self.GM = self.GM(:, logical(BC_v)); % Remove columns
        self.GK = self.GK(logical(BC_v), :); % Remove rows
        self.GK = self.GK(:, logical(BC_v)); % Remove columns
        self.mesh.BC_v = BC_v; % vector with boundary conditions

        % ID matrix
        self.mesh.ID = ones(2,self.mesh.n_nodes);
        for l = 1:size(self.data.BC)
            self.mesh.ID(1:2,self.data.BC(l,1)) = [self.data.BC(l,2:3)]';
        end

        p = 0; % aux variable for ID matrix
        for l = 1:self.mesh.n_nodes
            for m = 1:2
                if self.mesh.ID (m,l) >= 1
                    p = p+1;
                    self.mesh.ID(m,l) = self.mesh.ID(m,l)+p-1;
                end
            end
        end
        % -----------------------------------------------------------------

    end % end the matrix
%% Solving the eigen problem ==============================================
    function self = EigenSolve(self,BT)
    
        r = self.data.r_t; s = self.data.s_t; %

        % Relationship of the matrix
        % Mat = self.GM\self.GK;
    
        % Computing the eigenvectors and eigenvalues
        [a_vet,a_val] = eig(self.GM\self.GK); %Mat
    
        % Natural frequencies
        omega = sqrt(diag(a_val));
    
        % Arranging the frequencies in ascending order
        [self.freq,Ind] = sort(omega);
    
        % Natural frequencie (hz)
        self.natfreqHz = self.freq./(2*pi);
    
        % Critical parameters ---------------------------------------------
        % critical b eigenvalue
        self.critical(3) = 108; % 1/(r*s)

        % Critical frequency
        self.critical(1) = sqrt((self.Beam(1,3).*self.Beam(1,7).* ...
                        self.Beam(1,1))./(self.Beam(1,5).*self.Beam(1,2)));

        % Critical wavenumber = Critical frequency [Hz]
        self.critical(2) = (1/(2*pi))*sqrt((self.Beam(1,3).*self.Beam(1,7).* ...
                        self.Beam(1,1))./(self.Beam(1,5).*self.Beam(1,2)));
        % -------------------------------------------------------------

        % Solution eigenvalues and Dispersion relations ===================
        B = sqrt((self.Beam(1,5)*self.Beam(1,1)*self.data.L^4)/ ...
                           (self.Beam(1,4)*self.Beam(1,2)));
        b = self.freq*B; 

        if BT == 1 % Euler-Bernoulli --------------------------------------
             % Beta
            self.eigenvalues(:,1) = sqrt(b);

            self.natfreq(:,1) = B*(b./(2*pi));

        elseif BT == 2 % Rayleight Beam -----------------------------------
            % Beta
            self.eigenvalues(:,1) = sqrt((b.^2/2).*(r.^2 + sqrt(r.^4 + 4./b.^2)));
            % Alpha
            self.eigenvalues(:,2) = sqrt((b.^2/2).*(-r.^2 + sqrt(r.^4 + 4./b.^2)));

            self.natfreq(:,1) = B*(sqrt(self.eigenvalues(:,1).^2 - self.eigenvalues(:,2).^2)/2*pi*r);

        elseif BT == 3 % Shear Beam ---------------------------------------
            % Beta
            self.eigenvalues(:,1) = sqrt((b.^2/2).*(s.^2 + sqrt(s.^4 + 4./b.^2)));
            % Alpha
            self.eigenvalues(:,2) = sqrt((b.^2/2).*(-s.^2 + sqrt(s.^4+ 4./b.^2)));

            self.natfreq(:,1) = B*(sqrt(self.eigenvalues(:,1).^2 - self.eigenvalues(:,2).^2)/2*pi*s);

        elseif BT == 4 % Timoshenko Beam ----------------------------------
            % Beta
            self.eigenvalues(:,1) = sqrt((b.^2/2).*((r.^2+s.^2) + sqrt((r.^2-s.^2)^2 + 4./b.^2)));
            % Alpha before critical frequency
            self.eigenvalues(:,2) = sqrt((b.^2/2).*(-(r.^2+s.^2) + sqrt((r.^2-s.^2)^2 + 4./b.^2)));
            % Alpha after the critical frequency
            self.eigenvalues(:,3) = sqrt((b.^2/2).*((r.^2+s.^2) - sqrt((r.^2-s.^2)^2 + 4./b.^2)));

            for j = 1:size(b)
                if b(j) < self.critical(3)
                    % Natural frequency before critical frequency
                    self.natfreq(j,1) = B*(sqrt(self.eigenvalues(j,1).^2 - self.eigenvalues(j,2).^2)/ ...
                        (2*pi*sqrt(r.^2+s.^2)));
                elseif b(j) > self.critical(3)
                    % Natural frequency after critical frequency
                    self.natfreq(j,2) = B*(sqrt(self.eigenvalues(j,1).^2 + self.eigenvalues(j,3).^2)/ ...
                        (2*pi*sqrt(r.^2+s.^2)));
                end 
            end

        end % =============================================================

        % %Dispersion relationship (wave number) ========================
        % 
        % if BT == 1 %Euler-Bernoulli -----------------------------------
        %     self.wavenumber = sqrt(self.freq.*self.data.L.*sqrt(self.Beam(5)/...
        %         self.Beam(4))*self.S(1));
        % 
        % elseif BT == 2 %Rayleight Beam --------------------------------
        %     B = sqrt(self.Beam(5)*self.Beam(1)/(self.Beam(4)*self.Beam(2)))*...
        %         self.freq*self.data.L^2;
        % 
        %     % Alpha
        %     self.wavenumber(:,1) = (B./sqrt(2)).*sqrt(-(1/self.S(1))^2 + sqrt((1/self.S(1))^4 + 4./B.^2));
        %     % Beta
        %     self.wavenumber(:,2) = (B./sqrt(2)).*sqrt((1/self.S(1))^2 + sqrt((1/self.S(1))^4 + 4./B.^2));
        % 
        % elseif BT == 3 %Shear Beam ------------------------------------
        %     B = sqrt(self.Beam(5)*self.Beam(1)/(self.Beam(4)*self.Beam(2)))*...
        %         self.freq*self.data.L^2;
        % 
        %     % Alpha
        %     self.wavenumber(:,1) = (B./sqrt(2)).*sqrt(-(1/self.S(1)*self.gama(1))^2+...
        %         sqrt((1/self.S(1)*self.gama(1))^4 + 4./B.^2));
        %     % Beta
        %     self.wavenumber(:,2) = (B./sqrt(2)).*sqrt((1/self.S(1)*self.gama(1))^2+...
        %         sqrt((1/self.S(1)*self.gama(1))^4 + 4./B.^2));
        % 
        % 
        %             % 1- Cross-Section Area; 2-Moment of Inertia; 3- Shear Correction Factor, 
        % % 4- Young's Modulus, 5- Density, 6- Poisson's Modulus, 
        % % 7- Transversal Elasticity Modulus, 8- Rotatory Inertia Parameter
        % % 9- Shear deformation parameter self.Beam()
        % 
        % elseif BT == 4 %Timoshenko Beam -------------------------------
        %     B = sqrt(self.Beam(5)*self.Beam(1)/(self.Beam(4)*self.Beam(2)))*...
        %         self.freq*self.data.L^2;
        % 
        %     % Alpha before critical frequencie
        %     self.wavenumber(:,3) = (B./sqrt(2)).*sqrt(((1/self.S(1))^2 +...
        %         ((1/self.S(1))*self.gama(1))^2) - sqrt( ((1/self.S(1))^2 -...
        %         ((1/self.S(1))*self.gama(1))^2)^2 + 4./B.^2));
        % 
        %     % Alpha after the critical frequencie
        %     self.wavenumber(:,2) = (B./sqrt(2)).*sqrt(-((1/self.S(1))^2 +...
        %         ((1/self.S(1))*self.gama(1))^2) + sqrt( ((1/self.S(1))^2 -...
        %         ((1/self.S(1))*self.gama(1))^2)^2 + 4./B.^2));
        % 
        %     % Beta
        %     self.wavenumber(:,1) = (B./sqrt(2)).*sqrt(((1/self.S(1))^2 +...
        %         ((1/self.S(1))*self.gama(1))^ 2) + sqrt( ((1/self.S(1))^2 -...
        %         ((1/self.S(1))*self.gama(1))^2)^2 + 4./B.^2));
        % 
        % end % =========================================================

        % Pre-allocation --------------------------------------------------
        A_vet = zeros(length(Ind)); % eigenvectors ordened
        % -----------------------------------------------------------------

        % arranging the eigenvectors in ascending order
        for i=1:length(Ind)
            A_vet(:,i)=a_vet(:,Ind(i));
        end

        % Adding the contribution of the loop constants
        cont_disp = sum(self.data.BC(:,2));
        cont_rot = sum(self.data.BC(:,3));

        % Displacement ----------------------------------------------------
        disp = nonzeros(self.mesh.ID(1,:));
        self.U_disp = zeros(length(disp)+cont_disp,length(A_vet));
        % Re-ordenate eigenvectors rows of displacement
        for i=1:length(disp)
            self.U_disp(i+self.mesh.BC_v(1),:) = A_vet(disp(i),:);
        end
        % Displacement normalization
        self.U_disp = self.U_disp./max(abs((self.U_disp)));
        for  j = 1:size(A_vet)
            if  self.U_disp(2,j) <= 0
                self.U_disp(:,j) = -self.U_disp(:,j);
            end
        end
        %------------------------------------------------------------------

        % Rotation --------------------------------------------------------
        rot = nonzeros(self.mesh.ID(2,:));
        self.U_rot = zeros(length(rot)+...
            cont_rot,length(A_vet));
        % Re-ordenate eigenvectors rows of rotation
        for i=1:length(rot)
            self.U_rot(i+self.mesh.BC_v(2),:) = A_vet(rot(i),:);
        end
        % Rotation normalization
        self.U_rot = self.U_rot./max(abs((self.U_rot)));
        for  j = 1:size(A_vet)
            if  self.U_rot(2,j) <= 0
                self.U_rot(:,j) = -self.U_rot(:,j);
            end
        end
        %------------------------------------------------------------------
    
    end % end FrequenciesComp


%% Legendre-Gauss Quadrature by Greg von Winckel (2004)
    function [x,w] = LGQ(~,N,a,b)
    
        % This script is for computing definite integrals using Legendre-Gauss
        % Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
        % [a,b] with truncation order N
        %
        % Suppose you have a continuous function f(x) which is defined on [a,b]
        % which you can evaluate at any x in [a,b]. Simply evaluate it at all of
        % the values contained in the x vector to obtain a vector f. Then compute
        % the definite integral using sum(f.*w);
        %
        % Written by Greg von Winckel in 02/25/2004
    
        N=N-1;
        N1=N+1; N2=N+2;
        xu=linspace(-1,1,N1)';
        % Initial guess
        y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
        % Legendre-Gauss Vandermonde Matrix
        L=zeros(N1,N2);
        % Derivative of LGVM
        Lp=zeros(N1,N2);
        % Compute the zeros of the N+1 Legendre Polynomial
        % using the recursion relation and the Newton-Raphson method
        y0=2;
        % Iterate until new points are uniformly within epsilon of old points
        while max(abs(y-y0))>eps
    
            % Removed Lp for optimization, by slepton
            L(:,1)=1;
            L(:,2)=y;
    
            for k=2:N1
                L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
            end
    
            Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
    
            y0=y;
            y=y0-L(:,N2)./Lp;
    
        end
        % Linear map from[-1,1] to [a,b]
        x=(a*(1-y)+b*(1+y))/2;
        % Compute the weights
        w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
    
    end % end Legendre-Gauss Quadrature

        % %% Ploting the mode shapes - displacement -------------------------
        % function self = PlotStructuralDisplacements(self,modes,BT,material)
        % 
        %     % Removing rigid body movement
        %     RBM = length(self.omega_n(self.omega_n<4));
        % 
        %     % Data display
        %     frequencies(1,:) = self.wave_number(1+RBM:modes+RBM);
        %     display(frequencies)
        %     % fprintf('The number of rigid body moviment is: %d\n\n',RBM);
        %     fprintf('The factor ri is: %i and si is: %i\n\n',...
        %         1/self.S(1),self.gama(1)/self.S(1))
        %     fprintf('The factor r is: %i and s is: %i\n\n',...
        %         self.r(1),self.s(1))
        % 
        %     if BT == 1
        %         model = 'Euler-Bernoulli';
        %     elseif BT == 2
        %         model = 'Rayleigh';
        %     elseif BT == 3
        %         model = 'Shear';
        %     elseif BT == 4
        %         model = 'Timoshenko';
        %     end
        %     % 
        %     % % Normal mode shape -------------------------------------------
        %     % figure
        %     % hold on
        %     % box on
        %     % grid on
        %     % title(sprintf('First %d mode of a %s beam:',modes,model));
        %     % xlabel('\xi')
        %     % ylabel('V(\xi)')
        %     % axis([0 max(self.mesh.coordinates)/self.mesh.Ly -1.5 1.5])
        %     % for  mode_cnt = 1:modes
        %     %     plot(self.mesh.coordinates/self.mesh.Ly,...
        %     %         self.U_disp(:,mode_cnt + RBM),...
        %     %     "DisplayName",sprintf('%dº mode',mode_cnt));
        %     %     legend('-DynamicLegend');
        %     % end
        %     % legend ('show');
        %     % legend('Location','southwest')
        %     % hold off
        %     % %--------------------------------------------------------------
        %     % 
        %     % % Rotation mode shape -----------------------------------------
        %     % figure
        %     % hold on
        %     % box on
        %     % grid on
        %     % title(sprintf('First %d rotation modes of a %s beam:',modes,model));
        %     % xlabel('\xi')
        %     % ylabel('V(\xi)')
        %     % axis([0 max(self.mesh.coordinates)/self.mesh.Ly -1.5 1.5])
        %     % for  mode_cnt = 1:modes
        %     %     plot(self.mesh.coordinates/self.mesh.Ly,...
        %     %         self.U_rot(:,mode_cnt + RBM),...
        %     %         "DisplayName",sprintf('%dº mode',mode_cnt))
        %     %     legend('-DynamicLegend');
        %     % end
        %     % legend ('show');
        %     % legend('Location','southwest')
        %     % hold off
        %     % %--------------------------------------------------------------
        % 
        %     % % Ploting all theories for this boundary condition ------------
        %     % figure;
        %     % box on
        %     % grid on
        %     % title(sprintf('First %d Mode Shapes of all theories',modes));
        %     % xlabel('\xi')
        %     % ylabel('V(\xi)')
        %     % for i = 1:4
        %     %     BT = i;
        %     % 
        %     %     % Reruning the code for each beam theory
        %     %     self = AddSolidMaterial(self,BT,material);
        %     %     self = AddMefMatrix(self,BT);
        %     % 
        %     %     % Ploting the value for this theory
        %     %     for  mode_cnt = 1:modes
        %     %         subplot(modes,1,mode_cnt);
        %     %         hold on
        %     %         p = plot(self.mesh.coordinates/self.mesh.Ly,self.U_disp(:,mode_cnt + RBM));
        %     %         title(sprintf('mode %d',mode_cnt));
        %     %         axis([0 max(self.mesh.coordinates)/self.mesh.Ly -1.5 1.5])
        %     % 
        %     %         if i == 1
        %     %             p.LineStyle = '-';
        %     %             p.DisplayName = 'Euler';
        %     %         elseif i == 2
        %     %             p.LineStyle = '-.';
        %     %             p.DisplayName = 'Rayleigh';
        %     %         elseif i == 3
        %     %             p.LineStyle = ':';
        %     %             p.DisplayName = 'Shear';
        %     %         elseif i == 4
        %     %             p.LineStyle = '--';
        %     %             p.DisplayName = 'Timoshenko';
        %     %         end
        %     % 
        %     %     end
        %     % 
        %     % end
        %     % legend('Location','southoutside','Orientation','horizontal')
        %     % hold off
        % 
        %     % % % EBT x TBT for various frequencies ------------
        %     % % figure;
        %     % % box on
        %     % % grid on
        %     % % 
        %     % % xlabel('\xi')
        %     % % ylabel('V(\xi)')
        %     % % for i = [1,4]
        %     % %     BT = i;
        %     % % 
        %     % %     % Reruning the code for each beam theory
        %     % %     self = AddSolidMaterial(self,BT,material);
        %     % %     self = AddMefMatrix(self,BT);
        %     % % 
        %     % %     % Ploting the value for this theory
        %     % %     for  mode_cnt = 1:modes
        %     % %         subplot(modes,1,mode_cnt);
        %     % %         hold on
        %     % %         p = plot(self.mesh.coordinates/self.mesh.Ly,self.U_disp(:,mode_cnt + RBM));
        %     % %         title(sprintf('mode %d',mode_cnt));
        %     % %         axis([0 max(self.mesh.coordinates)/self.mesh.Ly -1.5 1.5])
        %     % % 
        %     % %     if i == 1
        %     % %         p.LineStyle = '-';
        %     % %         p.DisplayName = 'Euler';
        %     % %     elseif i == 4
        %     % %         p.LineStyle = '--';
        %     % %         p.DisplayName = 'Timoshenko';
        %     % %     end
        %     % % 
        %     % %     end
        %     % % 
        %     % % end
        %     % % legend('Location','southoutside','Orientation','horizontal')
        %     % % hold off
        % 
        %     % % EBT x TBT for various frequencies (Rotation) ------------
        %     % figure;
        %     % box on
        %     % grid on
        %     % title(sprintf('First %d Mode Shapes of EBT and TBT theories (Rotation)',modes));
        %     % xlabel('\xi')
        %     % ylabel('V(\xi)')
        %     % for i = [1,4]
        %     %     BT = i;
        %     % 
        %     %     % Reruning the code for each beam theory
        %     %     self = AddSolidMaterial(self,BT,material);
        %     %     self = AddMefMatrix(self,BT);
        %     % 
        %     %     % Ploting the value for this theory
        %     %     for  mode_cnt = 1:modes
        %     %         subplot(modes,1,mode_cnt);
        %     %         hold on
        %     %         p = plot(self.mesh.coordinates/self.mesh.Ly,self.U_rot(:,mode_cnt + RBM));
        %     %         title(sprintf('mode %d',mode_cnt));
        %     %         axis([0 max(self.mesh.coordinates)/self.mesh.Ly -1.5 1.5])
        %     % 
        %     %     if i == 1
        %     %         p.LineStyle = '-';
        %     %         p.DisplayName = 'Euler';
        %     %     elseif i == 4
        %     %         p.LineStyle = '--';
        %     %         p.DisplayName = 'Timoshenko';
        %     %     end
        %     % 
        %     %     end
        %     % 
        %     % end
        %     % legend('Location','southoutside','Orientation','horizontal')
        %     % hold off
        % 
        % end % End of the plotting

    end % end methods =====================================================
end % end classdef
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%