classdef calculations % v3.3
    properties                                                             
    data                        % beam data (struct)for each beam
    mesh                        % FEM mesh (struct) for each beam
    Beam                        % beam parameters for each beam
    layer                       % multilayer data
    structure                   % final structure
    coupled_mesh                % coupled mesh for plottings
    modes                       % number of shape modes to be ploted
    
    r                           % rotatory inertia parameter for the entire beam
    s                           % shear deformation parameter for the entire beam
    S                           % slenderness ratio = 1/r 
    gama                        % \frac{E}{\kappa G}
    L                           % sum of lenght of all coupled beams

    matrices
    % .K                        % global FEM stiffiness
    % .M                        % global FEM mass matrix
    % .CM                       % coupled mass matrix for all beams
    % .CK                       % coupled stiffness matrices
    
    U_disp                      % displacement matrix for ploting
    U_rot                       % rotation matrix for plotin
    
    result
    % .eigenvalues              % adimensional wave number
    % .natfreq                  % omega_n
    % .natfreqHz                % omega_n in [Hz]
    % .freq                     % dimensional omega_n 
    critical                    % critical parameters

    end 
    % ============================== FEM ================================ %
                                                                    methods
    %% CONSTRUCTOR                                                         
    function self = calculations(data,BT,modes)
        self.modes = modes;
        
        % Initialize self.data
        self.data = data;

        self.L = sum([data.L]); % sum of lenght of all coupled beams

        for input = 1:size(data,2)
            self.data(input) = data(input);
            
            self = FEM_setup(self,input);
            self = parameters(self,BT,input);
        end

        % for input = 1:size(data,2)
        %     self.data(input) = data(input);
        % 
        %     self = FEM_setup(self,input);
        %     self = parameters(self,BT,input);
        % end

        self = Eigenproblem(self); % setting the matrices
        self = Coupling(self); % coupling all beam matrices
        self = EigenSolve(self,BT); % solving for the frequencies
    end

    %% FEM SETTING                                                         
    function self = FEM_setup(self,input)

        self.mesh(input).n_el = self.data(input).n_el; % number of elements
        self.mesh(input).n_nodes = self.data(input).n_el+1;  % number of nodes
        self.mesh(input).l_el = (self.data(input).L/self.mesh(input).n_el)/2; % lenght of an isolated el

        % Coordinates increasing n_el times by n_el or L/n_el 
        self.mesh(input).coordinates = [linspace(0,self.data(input).L,self.mesh(input).n_nodes)]';
    end

    %% PARAMETERS #####                                                    
    function self = parameters(self,BT,input)

        self = self.geometry(self,input);
        self = self.material(self,input);

        % calculate each parameters for each element
        for el = 1:self.mesh(input).n_el
            E = self.Beam(input).E(el);
            A = self.Beam(input).A(el,1);
            I = self.Beam(input).I(el,1);
            k = self.Beam(input).k(el,1);
            G = self.Beam(input).G(el);
            l_e = self.mesh(input).l_el;
            % _____________________________________________________________
            
            % r and s parameters 
            % -------------------------------------------------------------
            if BT == 1 % Euler BT ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
                r_el = 0; % r for the element
                s_el = 0; % s for the element
            % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
    
            elseif BT == 2 % Rayleight BT ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
                r_el = sqrt(I./ (A.* l_e^2));
                s_el = 0;
            % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
    
            elseif BT == 3 % Shear BT ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
                r_el = 0;
                s_el = sqrt((E./ (k.* G)).*(sqrt(I./ (A.* l_e^2))).^2);
            % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
    
            elseif BT == 4 % Timoshenko BT  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
                % r_el = sqrt(I./ (A.* l_e^2));
                % s_el = sqrt((E./ (k.* G)).* r_el.^2);
                r_el = sqrt(I./ (A.* l_e^2));
                s_el = sqrt((E./(k.*G)).*r_el.^2);
                
            end
            % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
            % -------------------------------------------------------------
                                         self.Beam(input).r_e(el,1) = r_el;
                                         self.Beam(input).s_e(el,1) = s_el;
            % _____________________________________________________________
                                                                          %
            % This configuration allows allocation of all the data needed for 
            % modeling of conical, stepped, graduated beams and beams with other 
            % peculiarities.
        end

        % _________________________________________________________________
        E = self.Beam(input).E(1);
        A = self.Beam(input).A(1);
        I = self.Beam(input).I(1);
        k = self.Beam(input).k(1);

                                  % r and s for the entire beam - - - - - -
                            self.r(input) = sqrt(I(1)./ (A(1).* self.L^2));
           self.s(input) = sqrt((E(1)./ (k(1).* G(1))).* self.r(input).^2);
                                  % S and gamma - - - - - - - - - - - - - -
                                               self.S(input) = 1/self.r(1);
                             self.gama(input) = sqrt(E(1)./ (k(1).* G(1)));
                                                % self.gama(input) = 2.205;
                                  % - - - - - - - - - - - - - - - - - - - -
        % _________________________________________________________________
                                                  
        % Final structure _________________________________________________
        if input == 1 % pre-allocation 
            self.structure = self.Beam(input);
        end

        if self.data(input).geo == 1  % if normal geometry
            self.structure(input) = self.Beam(input);
        else % add layer ###### tirar o geo
            beam = self.Beam(input); 
            self.structure(input) = struct(...
                'k', beam.k, ...
                'G', beam.G, ...
                'r_e', beam.r_e, ...
                's_e', beam.s_e, ...
                'nu', beam.nu, ...
                'A', beam.A + sum(self.layer(input).A, 2), ...
                'I', beam.I + sum(self.layer(input).I, 2), ...
                'rho', beam.rho + sum(self.layer(input).rho, 2), ...
                'E', beam.E + sum(self.layer(input).E, 2) );
        end
        % _________________________________________________________________

    end % end of function parameters

    %% EIGENPROBLEM                                                        
    function self = Eigenproblem(self)

        % Compute (x) and weight(w) for integration  ----------------------
        [x,w] = self.LGQ(self,4,-1,1);

        % Shape functions pre-allocation ----------------------------------
        Nd = zeros(4); Ns = Nd; Nd_d = Nd; Ns_d = Nd;
        % Matrix cells pre-allocation
        self.matrices.K = cell(1,size(self.structure,2));
        self.matrices.M = self.matrices.K;

        % Loop for each different beam
        % _________________________________________________________________
        for inp = 1:length(self.data)

            % Stiffness global matrix pre-allocation
            GK = zeros(2*self.mesh(inp).n_nodes, 2*self.mesh(inp).n_nodes);
            % Mass global matrix
            GM = GK;
            % aux for allocating the elementar matrices on global
            satus = 1;

            L_e = self.mesh(inp).l_el;
            for i = 1:self.mesh(inp).n_el % Loop for each element
            % =============================================================

                % Shape functions 
                % ---------------------------------------------------------
                s_e = self.structure(inp).s_e(i);
                fc = 1/(4*(3*s_e^2+1)); % Constant for the function
                for j = 1:length(x)
    
                    % functions for deflection  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
                    Nd(1,j) = (2*(3*s_e^2+1)-3*(2*s_e^2+1)*x(j)+1*x(j)^3);
                    Nd(2,j) = L_e*((3*s_e^2+1)-1*x(j)-(3*s_e^2+1)*x(j)^2+1*x(j)^3);
                    Nd(3,j) = (2*(3*s_e^2+1)+3*(2*s_e^2+1)*x(j)-1*x(j)^3);
                    Nd(4,j) = L_e*(-(3*s_e^2+1)-1*x(j)+(3*s_e^2+1)*x(j)^2+1*x(j)^3);
                    Nd(:,j) = Nd(:,j)*fc;
                    % functions for bending ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ 
                    Ns(1,j) = (3*x(j)^2-3)/L_e;
                    Ns(2,j) = (-1-2*(3*s_e^2+1)*x(j)+6*s_e^2+3*x(j)^2);
                    Ns(3,j) = (3-3*x(j)^2)/L_e;
                    Ns(4,j) = (-1+2*(3*s_e^2+1)*x(j)+6*s_e^2+3*x(j)^2);
                    Ns(:,j) = Ns(:,j)*fc;
                    % first derivate of Nd  ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ 
                    Nd_d(1,j) = (-3*(2*s_e^2+1)+3*x(j)^2);
                    Nd_d(2,j) = L_e*(-1-2*(3*s_e^2+1)*x(j)+3*x(j)^2);
                    Nd_d(3,j) = (3*(2*s_e^2+1)-3*x(j)^2);
                    Nd_d(4,j) = L_e*(-1+2*(3*s_e^2+1)*x(j)+3*x(j)^2);
                    Nd_d(:,j) = Nd_d(:,j)*fc;
                    % first derivate of Ns¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
                    Ns_d(1,j) = (6*x(j))/L_e;
                    Ns_d(2,j) = (-2*(3*s_e^2+1)+6*x(j));
                    Ns_d(3,j) = (-6*x(j))/L_e;
                    Ns_d(4,j) = (2*(3*s_e^2+1)+6*x(j));
                    Ns_d(:,j) = Ns_d(:,j)*fc;
                    % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
                end 
                % ---------------------------------------------------------
            
                % Matrix constants 
                % ---------------------------------------------------------
                A = self.structure(inp).A(i);
                I = self.structure(inp).I(i);
                k = self.structure(inp).k(i);
                E = self.structure(inp).E(i);
                rho = self.structure(inp).rho(i);
                G = self.structure(inp).G(i);
                r_e = self.structure(inp).r_e(i);

                % constant for translation mass matrix
                C_tra = rho.* A* L_e;
                % constant for rotation mass matrix 
                C_rot = r_e.^2* rho.* A.* L_e^3;
                % constant for stiffness bending matrix (\frac{EI}{A})
                C_ben = (E.* I)./ L_e; 
                % constant for stiffness shear matrix
                C_she = (s_e^2.* (k.* G.* A).^2* L_e^3)./(E.* I); 
                % ---------------------------------------------------------
    
                % Computing the element matrix 
                % ---------------------------------------------------------
                M_tra = Nd * diag(w) * Nd';
                M_rot = Ns * diag(w) * Ns';
                K_ben = Ns_d * diag(w) * Ns_d';
                K_she = (Nd_d/self.mesh.l_el-Ns) * diag(w) * (Nd_d/self.mesh.l_el-Ns)';
                % (Nd_d/self.Le-Ns)
    
                % Complete mass element matrix
                M_m = M_tra.*C_tra + M_rot.*C_rot;
                % Complete stiffness element matrix
                M_k = K_ben.*C_ben + K_she.*C_she;
                % ---------------------------------------------------------
    
                % alocating the element matrices on the global
                % ---------------------------------------------------------
                finis = satus + 3;
                GM(satus:finis,satus:finis) = GM(satus:finis,satus:finis) + M_m(:,:);
                GK(satus:finis,satus:finis) = GK(satus:finis,satus:finis) + M_k(:,:);
                satus = satus + 2;
                % ---------------------------------------------------------

            end % end element loop
            % =============================================================
            clear A C_ben C_rot C_she C_tra E G I K_ben K_she L_e M_k M_m M_rot M_tra ...
            Nd Nd_d Ns Ns_d fc i j k r_e rho s_e BC_v % also reset BC_v

            % Boundary Conditions 
            % -------------------------------------------------------------
            % auxiliary vector that transforms BC conditions into a column array
            BC_c = reshape([self.data(inp).BC(:, 2), self.data(inp).BC(:, 3)].', [], 1);
        
            % defines the position of the dofs where the conditions are applied
            BC_pos = sort([2*self.data(inp).BC(:,1);2*self.data(inp).BC(:,1)-1]);
        
            % pre-allocation of vector with all dof boundary condition
            BC_v(1:2*self.mesh(inp).n_nodes,1) = 1;
        
            % adjust each dof condition 
            for l = 1:size(BC_c)
                BC_v(BC_pos(l)) = BC_c(l);
            end
        
            % recreate the global matrices with only unconstrained dofs
            % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
            GM = GM(logical(BC_v), :,:); % Remove rows 
            GM = GM(:, logical(BC_v),:); % Remove columns
            GK = GK(logical(BC_v), :,:); % Remove rows
            GK = GK(:, logical(BC_v),:); % Remove columns
            self.mesh(inp).BC_v = BC_v; % vector with boundary conditions
            % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨

            self.matrices.M{inp} = GM;
            self.matrices.K{inp} = GK;
            % -------------------------------------------------------------

            % ID matrix
            % -------------------------------------------------------------
            self.mesh(inp).ID = ones(2, self.mesh(inp).n_nodes);
            
            % Apply boundary conditions directly using vectorized approach
            bc_indices = self.data(inp).BC(:, 1);      % Nodes where BC applies
            bc_values = self.data(inp).BC(:, 2:3)';    % Corresponding BC values
            self.mesh(inp).ID(:, bc_indices) = bc_values;
            
            % Generate sequential IDs for free degrees of freedom
            free_dofs = (self.mesh(inp).ID == 1);      % Logical array for free DOFs
            self.mesh(inp).ID(free_dofs) = 1:nnz(free_dofs); % Assign sequential IDs
            % -------------------------------------------------------------


        end % end loop for different bems
        % _________________________________________________________________

    end % end the matrix
    
    %% COUPLING #####                                                      
    function self = Coupling(self)
        
        % get the size of matrices inside self.matrices.K
        N_sysn(:) = cellfun(@length, self.matrices.K);
        % pre-allocate the system global matrix
        self.matrices.CM = zeros(N_sysn(1) + sum(N_sysn(2:end)-2), N_sysn(1) + sum(N_sysn(2:end)-2));
        % self.matrices.CM = zeros(N_sysn(1) + sum(N_sysn(2:end)), N_sysn(1) + sum(N_sysn(2:end)));
        self.matrices.CK = self.matrices.CM;

        satus = 1; % aux for loccating the local matrices on the global
        for inp = 1:length(self.matrices.K)
            finis = satus + N_sysn(inp)-1; % same as satus

            % couples the mass matrices
            self.matrices.CM(satus:finis,satus:finis) = ...
                self.matrices.CM(satus:finis,satus:finis) + self.matrices.M{inp}(:,:);

            % couples the stifffness matrices
            self.matrices.CK(satus:finis,satus:finis) = ...
                self.matrices.CK(satus:finis,satus:finis) + self.matrices.K{inp}(:,:);

            satus = finis - 1;
        end
        clear satus finis inp N_sysn


        % coupling the mesh data for plotting the mode shapes
        % -----------------------------------------------------------------
        self.coupled_mesh.n_el = sum([self.mesh(:).n_el]); % number of elements
        self.coupled_mesh.n_nodes = self.coupled_mesh.n_el + 1;  % number of nodes
        self.coupled_mesh.l_el = (self.L/self.coupled_mesh.n_el)/2; % lenght of an isolated el

        % coordinates 
        self.coupled_mesh.coordinates = zeros(self.coupled_mesh.n_nodes, 1);
        
        % Flatten coordinates and apply offsets
        offset = 0; % Initialize offset
        current_index = 1; % Initialize insertion index
        
        for inp = 1:length(self.mesh)
            coords = self.mesh(inp).coordinates; % Current mesh coordinates
            start_idx = 1 + (inp > 1); % Skip the first point if not the first mesh
            n_points = numel(coords) - (inp > 1); % Adjust number of points accordingly
            self.coupled_mesh.coordinates(current_index:current_index + n_points - 1) = coords(start_idx:end) + offset; % Adjust and add
            offset = offset + coords(end); % Update offset with the last point of current mesh
            current_index = current_index + n_points; % Update insertion index
        end
        % -----------------------------------------------------------------

        % Coupled BC
        % ----------------------------------------------------------------
        % preallocate BC_c as a cell array
        BC_c = cell(length(self.data), 1);

        % fill BC_beam with boundary conditions
        for i = 1:length(self.data)
            if i == 1
                % first beam
                BC_c{i} = self.data(i).BC; % Último nó livre
            else
                % other beams
                BC_c{i} = [self.data(i).BC(2:end, 1) + self.data(i-1).BC(2:end, 1)-1,...
                      self.data(i).BC(2:end, 2:3)]; % Ignorar o primeiro nó
            end
        end

        % combine all boundary conditions into the global BC matrix
        self.coupled_mesh.BC = vertcat(BC_c{:});
        % -----------------------------------------------------------------

        % Coupled ID matrix
        % -----------------------------------------------------------------
        self.coupled_mesh.ID = ones(2, self.coupled_mesh.n_nodes);
        
        % Apply boundary conditions directly using vectorized approach
        bc_indices = self.coupled_mesh.BC(:, 1);       % Nodes where BC applies
        bc_values = self.coupled_mesh.BC(:, 2:3)';    % Corresponding BC values
        self.coupled_mesh.ID(:, bc_indices) = bc_values;
        
        % Generate sequential IDs for free degrees of freedom
        free_dofs = (self.coupled_mesh.ID == 1);      % Logical array for free DOFs
        self.coupled_mesh.ID(free_dofs) = 1:nnz(free_dofs); % Assign sequential IDs
        % -----------------------------------------------------------------

        % Coupled boundary conditions vector 
        % -----------------------------------------------------------------
        BC_c = reshape([self.coupled_mesh.BC(:, 2), self.coupled_mesh.BC(:, 3)].', [], 1);

        % defines the position of the dofs where the conditions are applied
        BC_pos = sort([2*self.coupled_mesh.BC(:,1);2*self.coupled_mesh.BC(:,1)-1]);

        % pre-allocation of vector with all dof boundary condition
        self.coupled_mesh.BC_v(1:2*self.coupled_mesh.n_nodes,1) = 1;

        % adjust each dof condition
        for l = 1:size(BC_c)
            self.coupled_mesh.BC_v(BC_pos(l)) = BC_c(l);
        end
        % -----------------------------------------------------------------
        
    end

    %% EIGENSOLVING                                                        

    function self = EigenSolve(self,BT)
    
        r_t = self.r(1); s_t = self.s(1); 
        A = self.structure(1).A(1);
        I = self.structure(1).I(1);
        k = self.structure(1).k(1);
        E = self.structure(1).E(1);
        rho = self.structure(1).rho(1);
        G = self.structure(1).G(1);
               
        % Computing the eigenvectors and eigenvalues
        [a_vet,a_val] = eig(self.matrices.CM\self.matrices.CK); 
    
        % Natural frequencies
        omega = sqrt(diag(a_val));
    
        % Arranging the frequencies in ascending order
        [self.result.natfreq,Ind] = sort(omega);
    
        % Natural frequencie (hz)
        self.result.natfreqHz = self.result.natfreq./(2*pi);
    
        % Critical parameters ---------------------------------------------
        % critical b eigenvalue
        self.critical.b = 1/(r_t*s_t);

        % Critical frequency
        self.critical.f = sqrt((k.*I.*A)./(rho.*I));

        % Critical wavenumber = Critical frequency [Hz]
        self.critical.omega = (1/(2*pi))*sqrt((k.*G.* A)./(rho.*I));
        % -------------------------------------------------------------

        % Solution eigenvalues and Dispersion relations ===================
        B = sqrt((rho* A* self.L^4)/ (E* I));
        b = self.result.natfreq*B; 

        % Solution eigenvalues and Dispersion relations ===================
        if BT == 1 % Euler-Bernoulli --------------------------------------
             % Beta
            self.result.dispersion(:,1) = sqrt(b);
            self.result.freq(:,1) = B*(b./(2*pi));

        elseif BT == 2 % Rayleight Beam -----------------------------------
            % Beta
            self.result.dispersion(:,1) = sqrt((b.^2/2).*(r_t.^2 + sqrt(r_t.^4 + 4./b.^2)));
            % Alpha
            self.result.dispersion(:,2) = sqrt((b.^2/2).*(-r_t.^2 + sqrt(r_t.^4 + 4./b.^2)));

            self.result.freq(:,1) = B*(sqrt(self.result.eigenvalues(:,1).^2 ...
                - self.result.dispersion(:,2).^2)/2*pi*r_t);

        elseif BT == 3 % Shear Beam ---------------------------------------
            % Beta
            self.result.dispersion(:,1) = sqrt((b.^2/2).*(s_t.^2 + sqrt(s_t.^4 + 4./b.^2)));
            % Alpha
            self.result.dispersion(:,2) = sqrt((b.^2/2).*(-s_t.^2 + sqrt(s_t.^4+ 4./b.^2)));

            self.result.freq(:,1) = B*(sqrt(self.result.eigenvalues(:,1).^2 ...
                - self.result.dispersion(:,2).^2)/2*pi*s_t);

        elseif BT == 4 % Timoshenko Beam ----------------------------------
            % Beta
            self.result.dispersion(:,1) = sqrt((b.^2/2).*((r_t.^2+s_t.^2) ...
                + sqrt((r_t.^2-s_t.^2)^2 + 4./b.^2)));
            % Alpha before critical frequency
            self.result.dispersion(:,2) = sqrt((b.^2/2).*(-(r_t.^2+s_t.^2) ...
                + sqrt((r_t.^2-s_t.^2)^2 + 4./b.^2)));
            % Alpha after the critical frequency
            self.result.dispersion(:,3) = sqrt((b.^2/2).*((r_t.^2+s_t.^2) ...
                - sqrt((r_t.^2-s_t.^2)^2 + 4./b.^2)));

            for j = 1:size(b)
                if b(j) < self.critical.b
                    % Natural frequency before critical frequency
                    self.result.freq(j,1) = B*(sqrt(self.result.dispersion(j,1).^2 ...
                        - self.result.dispersion(j,2).^2)/ ...
                        (2*pi*sqrt(r_t.^2+s_t.^2)));
                elseif b(j) > self.critical.b
                    % Natural frequency after critical frequency
                    self.result.freq(j,2) = B*(sqrt(self.result.dispersion(j,1).^2 ...
                        + self.result.dispersion(j,3).^2)/ ...
                        (2*pi*sqrt(r_t.^2+s_t.^2)));
                end 
            end

        end % =============================================================

        % %Dispersion relationship (wave number) ============================
        % 
        % if BT == 1 %Euler-Bernoulli ---------------------------------------
        %     self.result.eigenvalues = sqrt(self.result.natfreq.*self.L.*sqrt(rho/ E)*self.S(1));
        % 
        % elseif BT == 2 %Rayleight Beam ------------------------------------
        %     % Alpha
        %     self.result.eigenvalues(:,1) = (b./sqrt(2)).*sqrt(-(1/self.S(1))^2 + sqrt((1/self.S(1))^4 + 4./b.^2));
        %     % Beta
        %     self.result.eigenvalues(:,2) = (b./sqrt(2)).*sqrt((1/self.S(1))^2 + sqrt((1/self.S(1))^4 + 4./b.^2));
        % 
        % elseif BT == 3 %Shear Beam -----------------------------------------
        %     % Alpha
        %     self.result.eigenvalues(:,1) = (b./sqrt(2)).*sqrt(-(1/self.S(1)*self.gama(1))^2+...
        %         sqrt((1/self.S(1)*self.gama(1))^4 + 4./b.^2));
        %     % Beta
        %     self.result.eigenvalues(:,2) = (b./sqrt(2)).*sqrt((1/self.S(1)*self.gama(1))^2+...
        %         sqrt((1/self.S(1)*self.gama(1))^4 + 4./b.^2));
        % 
        % elseif BT == 4 %Timoshenko Beam -----------------------------------
        %     % Alpha before critical frequencie
        %     self.result.eigenvalues(:,3) = (b./sqrt(2)).*sqrt(((1/self.S(1))^2 +...
        %         ((1/self.S(1))*self.gama(1))^2) - sqrt( ((1/self.S(1))^2 -...
        %         ((1/self.S(1))*self.gama(1))^2)^2 + 4./b.^2));
        % 
        %     % Alpha after the critical frequencie
        %     self.result.eigenvalues(:,2) = (b./sqrt(2)).*sqrt(-((1/self.S(1))^2 +...
        %         ((1/self.S(1))*self.gama(1))^2) + sqrt( ((1/self.S(1))^2 -...
        %         ((1/self.S(1))*self.gama(1))^2)^2 + 4./b.^2));
        % 
        %     % Beta
        %     self.result.eigenvalues(:,1) = (b./sqrt(2)).*sqrt(((1/self.S(1))^2 +...
        %         ((1/self.S(1))*self.gama(1))^ 2) + sqrt( ((1/self.S(1))^2 -...
        %         ((1/self.S(1))*self.gama(1))^2)^2 + 4./b.^2));
        % 
        % end % =============================================================

        % rigid body movement
        % Rbm = length(self.result.natfreq(self.result.natfreq<8));
        %
        % % FREQ DISPLAY
        % frequencies(1,:) = self.result.eigenvalues(1+Rbm:self.modes+Rbm);
        % display(frequencies)

        % clear

        % Pre-allocation --------------------------------------------------
        A_vet = zeros(length(Ind)); % eigenvectors ordened
        % -----------------------------------------------------------------

        % arranging the eigenvectors in ascending order
        for i=1:length(Ind)
            A_vet(:,i)=a_vet(:,Ind(i));
        end

        % Adding the contribution of the loop constants
        BC_all = vertcat(self.coupled_mesh.BC);
        cont_disp = sum(1 - BC_all(:,2));
        cont_rot = sum(1 - BC_all(:,3));

        % Displacement ----------------------------------------------------
        disp = nonzeros(self.coupled_mesh.ID(1,:));
        self.U_disp = zeros(length(disp) + cont_disp,length(A_vet));
        % Re-ordenate eigenvectors rows of displacement
        for i=1:length(disp)
            self.U_disp(i+1-self.coupled_mesh.BC_v(1),:) = A_vet(disp(i),:);
        end
        % Displacement normalization
        self.U_disp = self.U_disp./max(abs(self.U_disp));%max()
        for  j = 1:size(A_vet)
            if  self.U_disp(2,j) <= 0
                self.U_disp(:,j) = -self.U_disp(:,j);
            end
        end
        % self.U_disp = self.U_disp.*(-1); for inverting   the plots
        % vertically. same on U_rot.
        %------------------------------------------------------------------

        % Rotation --------------------------------------------------------
        rot = nonzeros(self.coupled_mesh.ID(2,:));
        self.U_rot = zeros(length(rot) + cont_rot,length(A_vet));
        % Re-ordenate eigenvectors rows of rotation
        for i=1:length(rot)
            self.U_rot(i+1-self.coupled_mesh.BC_v(2),:) = A_vet(rot(i),:);
        end
        % Rotation normalization
        self.U_rot = self.U_rot./max(max(abs(self.U_rot)));
        for  j = 1:size(A_vet)
            if  self.U_rot(2,j) <= 0
                self.U_rot(:,j) = -self.U_rot(:,j);
            end
        end
        % self.U_rot = self.U_rot.*(-1);
        %------------------------------------------------------------------
    
    end % end FrequenciesComp

    end 
    % ============================ ANALYSIS ============================= %
                                                                    methods
    %% MODE SHAPES                                                         
    function self = ShapeModes(self,BT)
    
        % Removing rigid body movement
        Rbm = length(self.result.natfreq(self.result.natfreq<8));

        % % Data display
        % frequencies(1,:) = self.result.eigenvalues(1+Rbm:self.modes+Rbm);
        % display(frequencies)
        
        if BT == 1
            model = 'Euler-Bernoulli';
        elseif BT == 2
            model = 'Rayleigh';
        elseif BT == 3
            model = 'Shear';
        elseif BT == 4
            model = 'Timoshenko';
        end

        maxdisp = max(max(self.U_disp(:,1+ Rbm:self.modes+ + Rbm)))+0.05; 
        maxrot = max(max(self.U_rot(:,1+ Rbm:self.modes+ + Rbm)))+0.05; 
        mindisp = min(min(self.U_disp(:,1+ Rbm:self.modes+ + Rbm)))-0.05; 
        minrot = min(min(self.U_rot(:,1+ Rbm:self.modes+ + Rbm)))-0.05; 


        
        % Rotation mode shape -----------------------------------------
        figure("Position",[460 40 500 300])
        hold on
        box on
        grid on
        title(sprintf('First %d rotation modes of a %s beam:',self.modes,model));
        xlabel('\xi')
        ylabel('V(\xi)')
        axis([0 max(self.coupled_mesh.coordinates)/self.L minrot maxrot])
        for  mode_cnt = 1:self.modes
            plot(self.coupled_mesh.coordinates/self.L, self.U_rot(:,mode_cnt + Rbm),...
                "DisplayName",sprintf('%dº mode',mode_cnt))
            legend('-DynamicLegend');
        end
        legend ('show');
        legend('Location','southwest')
        hold off
        %--------------------------------------------------------------

                % Normal mode shape -------------------------------------------
        figure("Position",[460 342 500 300])
        hold on
        box on
        grid on
        title(sprintf('First %d mode of a %s beam:',self.modes,model));
        xlabel('\xi')
        ylabel('V(\xi)')
        axis([0 max(self.coupled_mesh.coordinates)/self.L mindisp maxdisp])
        for  mode_cnt = 1:self.modes
            plot(self.coupled_mesh.coordinates/self.L,  self.U_disp(:,mode_cnt + Rbm),  ...
                "DisplayName",sprintf('%dº mode',mode_cnt));
            legend('-DynamicLegend');
        end
        legend ('show');
        legend('Location','southwest')
        hold off
        % -------------------------------------------------------------

    end % End of the plotting

    %% ALL THEORIES SHAPE MODES COMPARING                                  
    function self = TheoriesModeShape(self)

        % Removing rigid body movement
        Rbm = length(self.result.natfreq(self.result.natfreq<8));

        % Ploting all theories for this boundary condition ------------
        figure;
        box on
        grid on
        title(sprintf('First %d Mode Shapes of all theories',self.modes));
        xlabel('\xi')
        ylabel('V(\xi)')
        for i = 1:4
            BT = i;

            % Reruning the code for each beam theory
            % for input = 1:size(self.data,2)
            %     self = parameters(self,BT,input);
            %     self = Eigenproblem(self);
            % end

            for input = 1:size(self.data,2)
                self = parameters(self,BT,input);
            end

            self = Eigenproblem(self); % setting the matrices
            self = Coupling(self); % coupling all beam matrices
            self = EigenSolve(self,BT); % solving for the frequencies

            % Ploting the value for this theory
            for  mode_cnt = 1:self.modes
                subplot(self.modes,1,mode_cnt);
                hold on
                p = plot(self.coupled_mesh.coordinates/self.L,self.U_disp(:,mode_cnt + Rbm));
                title(sprintf('mode %d',mode_cnt));
                axis([0 max(self.coupled_mesh.coordinates)/self.L -1.5 1.5])

                if i == 1
                    p.LineStyle = '-';
                    p.DisplayName = 'Euler';
                elseif i == 2
                    p.LineStyle = '-.';
                    p.DisplayName = 'Rayleigh';
                elseif i == 3
                    p.LineStyle = ':';
                    p.DisplayName = 'Shear';
                elseif i == 4
                    p.LineStyle = '--';
                    p.DisplayName = 'Timoshenko';
                end

            end

        end
        legend('Location','southoutside','Orientation','horizontal')
        hold off % end all theories

        % % ############
        % % EBT x TBT for various frequencies ---------------------------
        % figure;
        % box on
        % grid on
        % 
        % xlabel('\xi')
        % ylabel('V(\xi)')
        % for i = [1,4]
        %     BT = i;
        % 
        %     % Reruning the code for each beam theory
        %     self = parameters(self,BT);
        %     self = Eigenproblem(self);
        % 
        %     % Ploting the value for this theory
        %     for  mode_cnt = 1:self.modes
        %         subplot(self.modes,1,mode_cnt);
        %         hold on
        %         p = plot(self.mesh.coordinates/self.data.L,self.U_disp(:,mode_cnt + RBM));
        %         title(sprintf('mode %d',mode_cnt));
        %         axis([0 max(self.mesh.coordinates)/self.data.L -1.5 1.5])
        % 
        %         if i == 1
        %             p.LineStyle = '-';
        %             p.DisplayName = 'Euler';
        %         elseif i == 4
        %             p.LineStyle = '--';
        %             p.DisplayName = 'Timoshenko';
        %         end
        % 
        %     end
        % 
        % end
        % legend('Location','southoutside','Orientation','horizontal')
        % hold off % end RBT x TBT various freq rotation
        % 
        % % ############
        % % EBT x TBT for various frequencies (Rotation) ----------------
        % figure;
        % box on
        % grid on
        % title(sprintf('First %d Mode Shapes of EBT and TBT theories (Rotation)',self.modes));
        % xlabel('\xi')
        % ylabel('V(\xi)')
        % for i = [1,4]
        %     BT = i;
        % 
        %     % Reruning the code for each beam theory
        %     self = parameters(self,BT);
        %     self = Eigenproblem(self);
        % 
        %     % Ploting the value for this theory
        %     for  mode_cnt = 1:self.modes
        %         subplot(self.modes,1,mode_cnt);
        %         hold on
        %         p = plot(self.mesh.coordinates/self.data.L,self.U_rot(:,mode_cnt + RBM));
        %         title(sprintf('mode %d',mode_cnt));
        %         axis([0 max(self.mesh.coordinates)/self.data.L -1.5 1.5])
        % 
        %         if i == 1
        %             p.LineStyle = '-';
        %             p.DisplayName = 'Euler';
        %         elseif i == 4
        %             p.LineStyle = '--';
        %             p.DisplayName = 'Timoshenko';
        %         end
        % 
        %     end
        % 
        % end
        % legend('Location','southoutside','Orientation','horizontal')
        % hold off % end RBT x TBT various freq rotation
        
    end % end allthoeries function

    %% Rho x natFreq #####                                                 
    function self = Rho_natFreq(self,modes,data,type)

            % Suport variable
            second_spec = zeros(1,modes); % Is the change for alphas 
            Rbm = length(self.result.natfreq(self.result.natfreq<8));
            
            % Number of increases necessary for 1/S goes from 0 to 0.4
            nmr = 320; %320
            % Prealocating the frequencies and changes matrix
            frequencies = ones(nmr,modes*5);

            % Reseting the diameter
            rhov = linspace (1000,12000,nmr);

            for j = 1:nmr
                BT = 1;
        
                % ---------------------------------------------------------
                for input = 1:size(data,2)
                    self.data(input).rho{1}(1) = rhov(j);
                    self = parameters(self,BT,input);
                end
        
                self = Eigenproblem(self); % setting the matrices
                self = Coupling(self); % coupling all beam matrices
                self = EigenSolve(self,BT); % solving for the frequencies

                wv_n = self.result.eigenvalues; % wave number
                natF_Hz = self.result.natfreqHz;
                % ---------------------------------------------------------

                frequencies(j,1:modes) = wv_n(1+Rbm:modes+Rbm,1);

                frequencies(j,modes+1:2*modes) = natF_Hz(1+Rbm:modes+Rbm,1);

                % Calculating the frequencies for the next Beam theory ----
                BT = 1+type;

                % ---------------------------------------------------------
                for input = 1:size(data,2)
                    self.data(input).rho{1}(1) = rhov(j);
                    self = parameters(self,BT,input);
                end

                self = Eigenproblem(self); % setting the matrices
                self = Coupling(self); % coupling all beam matrices
                self = EigenSolve(self,BT); % solving for the frequencies

                wv_n = self.result.eigenvalues; % wave number
                natF = self.result.natfreq;
                natF_Hz = self.result.natfreqHz;
                omega_crit = self.critical.omega;
                % ---------------------------------------------------------

                A = self.structure(1).A(1);
                I = self.structure(1).I(1);
                E = self.structure(1).E(1);
                rho = self.structure(1).rho(1);

                % Alpha relation for Rayleight and shear beam, and Beta for
                % Timoshenko. Change made because Alpha of tbt change
                frequencies(j,2*modes+1:3*modes) = wv_n(1+Rbm:modes+Rbm,1);

                if type == 1
                    % Beta relation
                    frequencies(j,3*modes+1:4*modes) = wv_n(1+Rbm:modes+Rbm,2);

                    % Relation for frequency curves
                    frequencies(j,4*modes+1:5*modes) = natF_Hz(1+Rbm:modes+Rbm,1);

                elseif type == 2
                    % Beta relation
                    frequencies(j,3*modes+1:4*modes) = wv_n(1+Rbm:modes+Rbm,2);

                    % Relation for frequency curves
                    frequencies(j,4*modes+1:5*modes) = natF_Hz(1+Rbm:modes+Rbm,1);

                elseif type == 3
                    B = sqrt(rho*A/(E*I))*natF*self.L^2;

                    for i = 1:modes
                        if B(i+Rbm) > ((self.S(1))^2/self.gama(1))
                            second_spec(i) = 1;
                        end
                        % Alpha relation
                        frequencies(j,3*modes+i) = wv_n(i+Rbm,2+second_spec(i));
                    end

                    % Critical frequency
                    frequencies(j,5*modes+1) = omega_crit(1);

                    % Relation for frequency curves
                    frequencies(j,4*modes+1:5*modes) = natF_Hz(1+Rbm:modes+Rbm,1);

                end

            end

            % Ploting the frequencies -------------------------------------
            figure
            hold on
            box on
            if type == 1
                el2 = plot(rhov(:),frequencies(:,modes+1:2*modes),'-',...
                    'DisplayName','Euler','color','#EDB120','HandleVisibility', 'off');
                sl =  plot(rhov(:),frequencies(:,4*modes+1:5*modes),'--',...
                    'DisplayName','Rayleight','color','b','HandleVisibility', 'off');
            elseif type == 2
                el2 = plot(rhov(:),frequencies(:,modes+1:2*modes),'-',...
                    'DisplayName','Euler','color','#EDB120','HandleVisibility', 'off');
                sl =  plot(rhov(:),frequencies(:,4*modes+1:5*modes),'--',...
                    'DisplayName','Shear','color','b','HandleVisibility', 'off');
            elseif type == 3 
                el2 = plot(rhov(:),frequencies(:,modes+1:2*modes),'-',...
                    'DisplayName','Euler','color','#EDB120','HandleVisibility', 'off');
                sl =  plot(rhov(:),frequencies(:,4*modes+1:5*modes),'--',...
                    'DisplayName','Timoshenko','color','b','HandleVisibility', 'off');
            end

            % Legend ploting
            el2(1).HandleVisibility = 'on';
            sl(1).HandleVisibility = 'on';

            xlabel('Rho (Kg/m^3))')
            ylabel('\omega (Hz)')
            axis('tight')
            title('Natural Frequency x Specific Weight') 
            legend('show','Location','northeast') 

            hold off
            % -------------------------------------------------------------
    end
    %% E x natFreq #####                                                   
    function self = E_natFreq(self,modes,data,type)
            % Suport variable
            second_spec = zeros(1,modes); % Is the change for alphas 
            Rbm = length(self.result.natfreq(self.result.natfreq<8));
            
            nmr = 320; % number of repetitions
            % Prealocating the frequencies and changes matrix
            frequencies = ones(nmr,modes*5);

            % Reseting the diameter
            youngv = linspace (110e9,207e9,nmr);

            for j = 1:nmr
                BT = 1;

                % ---------------------------------------------------------
                for input = 1:size(data,2)
                    self.data(input).E{1}(1) = youngv(j);
                    self = parameters(self,BT,input);
                end
        
                self = Eigenproblem(self); % setting the matrices
                self = Coupling(self); % coupling all beam matrices
                self = EigenSolve(self,BT); % solving for the frequencies

                wv_n = self.result.eigenvalues; % wave number
                natF_Hz = self.result.natfreqHz;
                % ---------------------------------------------------------

                frequencies(j,1:modes) = wv_n(1+Rbm:modes+Rbm,1);

                % Relation for frequency curves
                if type == 1 || type == 3

                    frequencies(j,modes+1:2*modes) = natF_Hz(1+Rbm:modes+Rbm,1);
                    
                elseif type == 2

                    frequencies(j,modes+1:2*modes) = natF_Hz(1+Rbm:modes+Rbm,1);

                end

                % Calculating the frequencies for the next Beam theory ----
                BT = 1+type;

               % ---------------------------------------------------------
                for input = 1:size(data,2)
                    self.data(input).E{1}(1) = youngv(j);
                    self = parameters(self,BT,input);
                end

                self = Eigenproblem(self); % setting the matrices
                self = Coupling(self); % coupling all beam matrices
                self = EigenSolve(self,BT); % solving for the frequencies

                wv_n = self.result.eigenvalues; % wave number
                natF = self.result.natfreq; %#########
                natF_Hz = self.result.natfreqHz;
                omega_crit = self.critical.omega;

                A = self.structure(1).A(1);
                I = self.structure(1).I(1);
                E = self.structure(1).E(1);
                rho = self.structure(1).rho(1);
                % ---------------------------------------------------------

                % Alpha relation for Rayleight and shear beam, and Beta for
                % Timoshenko. Change made because Alpha of tbt change
                frequencies(j,2*modes+1:3*modes) = wv_n(1+Rbm:modes+Rbm,1);

                if type == 1
                    % Beta relation
                    frequencies(j,3*modes+1:4*modes) = wv_n(1+Rbm:modes+Rbm,2);

                    % Relation for frequency curves
                    frequencies(j,4*modes+1:5*modes) = natF_Hz(1+Rbm:modes+Rbm,1);

                elseif type == 2
                    % Beta relation
                    frequencies(j,3*modes+1:4*modes) = wv_n(1+Rbm:modes+Rbm,2);

                    % Relation for frequency curves
                    frequencies(j,4*modes+1:5*modes) = natF_Hz(1+Rbm:modes+Rbm,1);

                elseif type == 3
                    B = sqrt(rho*A/(E*I))*natF*self.L^2;

                    for i = 1:modes
                        if B(i+Rbm) > ((self.S(1))^2/self.gama(1))
                            second_spec(i) = 1;
                        end
                        % Alpha relation
                        frequencies(j,3*modes+i) = wv_n(i+Rbm,2+second_spec(i));
                    end

                    % Critical frequency
                    frequencies(j,5*modes+1) = omega_crit(1);

                    % Relation for frequency curves
                    frequencies(j,4*modes+1:5*modes) = natF_Hz(1+Rbm:modes+Rbm,1);

                end
            end

            % Ploting the frequencies -------------------------------------
            figure
            hold on
            box on
            if type == 1
                el2 = plot(youngv(:),frequencies(:,modes+1:2*modes),'-',...
                    'DisplayName','Euler','color','#EDB120','HandleVisibility', 'off');
                sl =  plot(youngv(:),frequencies(:,4*modes+1:5*modes),'--',...
                    'DisplayName','Rayleight','color','b','HandleVisibility', 'off');
            elseif type == 2
                el2 = plot(youngv(:),frequencies(:,modes+1:2*modes),'-',...
                    'DisplayName','Euler','color','#EDB120','HandleVisibility', 'off');
                sl =  plot(youngv(:),frequencies(:,4*modes+1:5*modes),'--',...
                    'DisplayName','Shear','color','b','HandleVisibility', 'off');
            elseif type == 3
                el2 = plot(youngv(:),frequencies(:,modes+1:2*modes),'-',...
                    'DisplayName','Euler','color','#EDB120','HandleVisibility', 'off');
                sl =  plot(youngv(:),frequencies(:,4*modes+1:5*modes),'--',...
                    'DisplayName','Timoshenko','color','b','HandleVisibility', 'off');
            end

            % Legend ploting
            el2(1).HandleVisibility = 'on';
            sl(1).HandleVisibility = 'on';

            xlabel('Young`s modulus (GPa)')
            ylabel('\omega (Hz)')
            axis('tight')
            title('Natural Frequency x Young´s modulus')
            legend('show','Location','northeast')

            hold off
            % -------------------------------------------------------------
    end
    %% Slenderness Ratio x natFreq BY SLEPTON (2023) #####                 
    function self = slendernessR_natFreq(self,modes,data,type)
            % Suport variable
            second_spec = zeros(1,modes); % Is the change for alphas 
            Rbm = length(self.result.natfreq(self.result.natfreq<8));
            
            % Number of increases necessary for 1/S goes from 0 to 0.4
            nmr = 1600; %1600
            % Prealocating the frequencies and changes matrix
            frequencies = ones(nmr,modes*5);
            change = ones(nmr,3);

            % Reseting the diameter
            d_change = zeros(length(self.data),2);

            for j = 1:nmr
                BT = 1;

                % Recalculating the beam with new diameter
                d_change = d_change(1,1)+0.001; %[m] diameter;

                % ---------------------------------------------------------
                for input = 1:size(data,2)
                    % self.data(input).d{1} = [d_change, d_change, 1, ...
                    %     self.data(input).d{1}(1,4), self.data(input).d{1}(1,5)];
                    self.data(input).d{1}(1:2) = d_change;
                    % will not work for self.data(input).d(1,3), form =
                    % hollow geometries
                    self = parameters(self,BT,input);
                end
        
                self = Eigenproblem(self); % setting the matrices
                self = Coupling(self); % coupling all beam matrices
                self = EigenSolve(self,BT); % solving for the frequencies

                wv_n = self.result.eigenvalues; % wave number
                natF = self.result.natfreq;

                A = self.structure(1).A(1);
                I = self.structure(1).I(1);
                E = self.structure(1).E(1);
                rho = self.structure(1).rho(1);
                % ---------------------------------------------------------

                frequencies(j,1:modes) = wv_n(1+Rbm:modes+Rbm,1);
                % Relation for frequency curves
                if type == 1 || type == 3

                    frequencies(j,modes+1:2*modes) = natF(1+Rbm:modes+Rbm,1)*...
                        (1/self.S(1))*sqrt(rho*A*self.L^4/(E*I));
                    
                elseif type == 2

                    frequencies(j,modes+1:2*modes) = natF(1+Rbm:modes+Rbm,1)*...
                        self.gama(1)*(1/self.S(1))*sqrt(rho*A*self.L^4/(E*I));

                end

                % Calculating the frequencies for the next Beam theory
                BT = 1+type;

                % ---------------------------------------------------------
                for input = 1:size(data,2)
                    % self.data(input).d{1} = [d_change, d_change, 1, ...
                    %     self.data(input).d{1}(1,4), self.data(input).d{1}(1,5)];
                    self.data(input).d{1}(1:2) = d_change; 
                    self = parameters(self,BT,input);
                end

                self = Eigenproblem(self); % setting the matrices
                self = Coupling(self); % coupling all beam matrices
                self = EigenSolve(self,BT); % solving for the frequencies

                wv_n = self.result.eigenvalues; % wave number
                natF = self.result.natfreq;
                % natF_Hz = self.result.natfreqHz;
                omega_crit = self.critical.omega;


                A = self.structure(1).A(1);
                I = self.structure(1).I(1);
                E = self.structure(1).E(1);
                rho = self.structure(1).rho(1);
                % ---------------------------------------------------------

                % Alpha relation for Rayleight and shear beam, and Beta for
                % Timoshenko. Change made because Alpha of tbt change
                frequencies(j,2*modes+1:3*modes) = wv_n(1+Rbm:modes+Rbm,1);

                if type == 1
                    % Beta relation
                    frequencies(j,3*modes+1:4*modes) = wv_n(1+Rbm:modes+Rbm,2);

                    % Relation for frequency curves
                    frequencies(j,4*modes+1:5*modes) = natF(1+Rbm:modes+Rbm,1)*...
                        (1/self.S(1))*sqrt(rho*A*self.L^4/(E*I));

                elseif type == 2
                    % Beta relation
                    frequencies(j,3*modes+1:4*modes) = wv_n(1+Rbm:modes+Rbm,2);

                    % Relation for frequency curves
                    frequencies(j,4*modes+1:5*modes) = natF(1+Rbm:modes+Rbm,1)*...
                        self.gama(1)*(1/self.S(1))*sqrt(rho*A*self.L^4/(E*I));

                elseif type == 3
                    B = sqrt(rho*A/(E*I))*natF*self.L^2;

                    for i = 1:modes
                        if B(i+Rbm) > ((self.S(1))^2/self.gama(1))
                            second_spec(i) = 1;
                        end
                        % Alpha relation
                        frequencies(j,3*modes+i) = wv_n(i+Rbm,2+second_spec(i));
                    end

                    % Critical frequency
                    frequencies(j,5*modes+1) = omega_crit(1);

                    % Relation for frequency curves
                    frequencies(j,4*modes+1:5*modes) = natF(1+Rbm:modes+Rbm,1)*...
                        (1/self.S(1))*sqrt(rho*A*self.L^4/(E*I));
                end

                % X-axis
                change(j,1) = 1/self.S(1);
                change(j,2) = self.gama(1)/self.S(1); 

            end

            % wave number curves ------------------------------------------
            figure
            hold on
            box on
            el = plot(change(:,1),frequencies(:,1:modes),'-',...
                'DisplayName','Euler','color','#EDB120','HandleVisibility', 'off');
            if type == 1
                al = plot(change(:,1),frequencies(:,2*modes+1:3*modes),'--',...
                    'DisplayName','\alpha','color','#A2142F','HandleVisibility', 'off');
                bl = plot(change(:,1),frequencies(:,3*modes+1:4*modes),'-.',...
                    'DisplayName','\beta','color','b','HandleVisibility', 'off');
                xlabel('1/S')
                ylabel('Wave Number')
                axis([0 0.3 0 16])
            elseif type == 2
                al = plot(change(:,2),frequencies(:,2*modes+1:3*modes),'--',...
                    'DisplayName','\alpha','color','#A2142F','HandleVisibility', 'off');
                bl = plot(change(:,2),frequencies(:,3*modes+1:4*modes),'-.',...
                    'DisplayName','\beta','color','b','HandleVisibility', 'off');
                xlabel('\gamma/S')
                ylabel('Wave Number')
                axis([0 0.3 0 16])
            elseif type == 3
                ex = plot(change(:,1),frequencies(:,5*modes+1),'-',...
                    'DisplayName','\beta_{crit}','color','#109411','HandleVisibility', 'off');
                bl = plot(change(:,1),frequencies(:,2*modes+1:3*modes),'--',...
                    'DisplayName','\beta','color','b','HandleVisibility', 'off');
                for i = 1:modes-1 % For all timoshenko alpha
                    plot(change(:,1),frequencies(:,3*modes+i),'-.',...
                        'color','#A2142F','HandleVisibility', 'off'); % alpha
                end
                al = plot(change(:,1),frequencies(:,3*modes+modes),'-.',...
                    'DisplayName','\alpha','color','#A2142F','HandleVisibility', 'off');
                xlabel('1/S')
                ylabel('Wave Number')
                axis([0 0.4 0 16])
            end

            % Legend ploting
            el(1).HandleVisibility = 'on';
            al(1).HandleVisibility = 'on';
            bl(1).HandleVisibility = 'on';
            ex(1).HandleVisibility = 'on';
            legend('show','Location','northeast')

            hold off
            % -------------------------------------------------------------

            % Ploting the frequencies curves ------------------------------

            figure
            hold on
            box on
            if type == 1
                el2 = plot(change(:,1),frequencies(:,modes+1:2*modes),'-',...
                    'DisplayName','Euler','color','#EDB120','HandleVisibility', 'off');
                sl =  plot(change(:,1),frequencies(:,4*modes+1:5*modes),'--',...
                    'DisplayName','Rayleight','color','b','HandleVisibility', 'off');
                xlabel('1/S')
                ylabel('\omega*1/S*(\rho*A*L^4/E*I)^{1/2}')
                axis([0 0.3 0 16])
            elseif type == 2
                el2 = plot(change(:,2),frequencies(:,modes+1:2*modes),'-',...
                    'DisplayName','Euler','color','#EDB120','HandleVisibility', 'off');
                sl =  plot(change(:,2),frequencies(:,4*modes+1:5*modes),'--',...
                    'DisplayName','Shear','color','b','HandleVisibility', 'off');
                xlabel('\gamma/S')
                ylabel('\omega*\gamma/S*(\rho*A*L^4/E*I)^{1/2}')
                axis([0 0.3 0 16])
            elseif type == 3
                el2 = plot(change(:,1),frequencies(:,modes+1:2*modes),'-',...
                    'DisplayName','Euler','color','#EDB120','HandleVisibility', 'off');
                sl =  plot(change(:,1),frequencies(:,4*modes+1:5*modes),'--',...
                    'DisplayName','Timoshenko','color','b','HandleVisibility', 'off');
                xlabel('1/S')
                ylabel('\omega*1/S*(\rho*A*L^4/E*I)^{1/2}')
                % axis('tight')
                axis([0 0.4 0 6])
            end

            % Legend ploting
            el2(1).HandleVisibility = 'on';
            sl(1).HandleVisibility = 'on';
            legend('show','Location','southeast')

            hold off
    end

    end
    % ============================= STATIC ============================== %
                                                           methods (Static)
    %% LEGENDRE-GAUSS QUADRATURE BY GREG VON WINCKEL (2004)                
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
    %% MATERIAL #                                                          
    function self = material(self,input)

        % reorganizing values for each element ____________________________
        % depois fzer outr método só pras seções bugadas e pra tapered, mudando o
        % formato do input, deixar os tipos de input como comentario echamar os
        % metodos de geo em cada input, depois fzer o mesmo pra fgmb
        %
        %secoes bugadas colocar d1, d2..., dn
        %tapered: definir função com o handle e fzer o calculo de cada nó no metodo
        fields = {'E', 'rho'};
        for f = 1:numel(fields)
            masu = vertcat(self.data(input).(fields{f}){:});
            [values, intervals] = deal(masu(:, 1), masu(:, end-1:end));

            % calculate the size of the vector
            new_values = zeros(self.data(input).n_el, 1); % pre-allocate new_values

            % fill new_values with the values to its respective nodes
            for i = 1:numel(values)
                new_values(intervals(i,1) : intervals(i,2)) = values(i);
            end

            % store the vector on its due field
            self.Beam(input).(fields{f}) = new_values;
        end
        
        % Transversal Elasticity (G)
        self.Beam(input).G(:,1) = self.Beam(input).E./(2*(1 + self.Beam(input).nu));
        % self.Beam(input).G(:,1) = ones(1,size(self.Beam(input).E,1))*77.5e9;
        % _________________________________________________________________
    end
    %% GEOMETRY #                                                          
    function self = geometry(self,input)
        % reorganizing values for each element ____________________________
        fields = {'d1', 'd2','d3','nu'};% d3 is the form
        
        for f = 1:numel(fields)
           if startsWith(fields{f}, 'd')
                % get if its d1 d2 or d3
                col_idx = str2double(fields{f}(2)); 
                % values for each node set
                masu = vertcat(self.data(input).d{:});
                % nodes where its is applied
                [values, intervals] = deal(masu(:, col_idx), masu(:, end-1:end));
           else
                masu = vertcat(self.data(input).(fields{f}){:});
                [values, intervals] = deal(masu(:, 1), masu(:, end-1:end));
            end

            % calculate the size of the vector
            new_values = zeros(self.data(input).n_el, 1); % pre-allocate new_values

            % fill new_values with the values to its respective nodes
            for i = 1:size(intervals,1)
                new_values(intervals(i,1) : intervals(i,2)) = values(i);
            end

            % store the vector on its due field
            self.Beam(input).(fields{f}) = new_values;
        end
        clear idx fields f new_values intervals values i col_idx;
        % _________________________________________________________________

        % -----------------------------------------------------------------
        for el = 1:self.mesh(input).n_el

            % rewritting the parameters for convenience
            d1 = self.Beam(input).d1(el);
            d2 = self.Beam(input).d2(el);
            nu = self.Beam(input).nu(el);

            if self.Beam(input).d3(el) == 1 % Circular¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
                A = pi/4*(d1.^2);
                I = pi/64*(d1.^4);
                k = 6*(1+nu)./(7+6*nu);
    
            elseif self.Beam(input).d3(el) == 2 % Rectangular ¨ ¨ ¨ ¨ ¨ ¨ ¨
                A = d1.*d2;
                I = (d1.*d2.^3)/12;
                % k = 10*(1+nu)./(12+11.*nu);
                k = 5/6;

            elseif self.Beam(input).d3(el) == 3 % Hollow Circle ¨ ¨ ¨ ¨ ¨ ¨ 
                A = pi*(d1^2 - d2^2);
                I = (pi/2)*(d1^4 - d2^4);
                k = (2 *(1+nu))/(4 + 3*nu);

            elseif self.Beam(input).d3(el) == 4 %Thin-Walled Square Tube¨ ¨
                A = d1^2 - d2^2;
                I = ((d1^4)/12) - ((d2^4)/12);
                k = (20*(1 + nu))/(48 + 39*nu);

            elseif self.Beam(input).d3(el) == 20 %Thin-Walled Square Tube¨ ¨
                A = 0.0097389;
                I = 0.0001171;
                k =  0.53066;

            end
            % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
            % -------------------------------------------------------------
                                              self.Beam(input).A(el,1) = A;
                                              self.Beam(input).I(el,1) = I;
                                              self.Beam(input).k(el,1) = k;
        end
        % -----------------------------------------------------------------

        % % modifications ___________________________________________________
        % if self.data(input).geo ==  1
        % else
        %     self = self.multilayer(self,input);
        % end
        % % _________________________________________________________________
    end
    %% DRAW #####                                                          
    function self = Draw(self)

        x_ax = self.coupled_mesh.coordinates;
 
        if any([self.data.geo] == 2)
            % cellData = cell(size(self.data, 2), 1);
            % for k = 1:size(self.data, 2)
            %     cellData{k} = self.layer(k).dl1(:, 1); % Extrai a primeira coluna
            % end
            % d1 = vertcat(cellData{:});
            % d1 = [d1./2; d1(end)./2]; d2 = d1;
            % 
            % cellData = cell(size(self.data, 2), 1);
            % for k = 1:size(self.data, 2)
            %     cellData{k} = self.layer(k).dl2(:, 1); % Extrai a primeira coluna
            % end
            % dl1 = vertcat(cellData{:});
            % dl1 = [dl1./2; dl1(end)./2];
            d1 = [vertcat(self.Beam(:).d1)./2; self.Beam(end).d1(end)./2];
            d2 = d1;
        else
            d1 = [vertcat(self.Beam(:).d1)./2; self.Beam(end).d1(end)./2];
            d2 = [vertcat(self.Beam(:).d2)./2; self.Beam(end).d2(end)./2];
        end
        

        rho = [vertcat(self.Beam(:).rho); self.Beam(end).rho(end)];
        E = [vertcat(self.Beam(:).E); self.Beam(end).E(end)];
        nu = [vertcat(self.Beam(:).nu); self.Beam(end).nu(end)];


        scr_sz = get( groot, 'Screensize' );
        figure("Position",[0 0 500 scr_sz(end)])
        % -----------------------------------------------------------------
        subplot(3,1,1)
        hold on
        box on
        grid off
        % xlabel('L (m)')
        ylabel('heigth (m)')
        title('Side view of the beam and mesh')
        %
        % plot the vertical dimension line - - -
        stairs(x_ax,  [d1,-d1],"Color","k","LineWidth",1.5);

        % plot the nodes - - -
        scatter(x_ax,zeros(length(x_ax)),10,'MarkerEdgeColor',"k")

        % plot the central line - - -
        plot([0;max(x_ax)], [0;0],"Color","k","LineWidth",0.2)

        % % plot layers line - - -
        % stairs(x_ax,  [dl1,-dl1],"Color","b");

        axis([0 max(x_ax) -max(d1)-(max(d1)) max(d1)+max(d1)])
        legend ('off');
        hold off
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        subplot(3,1,2)
        hold on
        box on
        grid off
        xlabel('L (m)')
        ylabel('width (m)')
        title('Top view of the beam and mesh')
        
        % plot the horizontal dimension line - - -
        stairs(x_ax,  [d2,-d2],"Color","k","LineWidth",1.7);
        
        % % plot layers line - - -
        % stairs(x_ax,  [dl1,-dl1],"Color","b");

        axis([0 max(x_ax) -max(d2)-(max(d2)) max(d2)+max(d2)])
        legend ('off');
        hold off
        % -----------------------------------------------------------------

        properties = {rho, E, nu};
        unities = {'Pa', 'kg/m^3', '~'};
        prop_name = {'\rho', 'E', '\nu'};
        colors =[1, 0.824, 0.408;  % Azul (Density)
                 0.62, 0.886, 1;  % Vermelho (Young Modulus)
                 0.706, 1, 0.322]; % Verde (Poisson Ratio)
        y_offsets = [0, -0.2, -0.4]; 
        % -----------------------------------------------------------------
        subplot(3,1,3)
        hold on
        box off
        grid off
        title('Material change from light to dark')
        % Material change colormap
        for i = 1:numel(properties)
            values = properties{i};
            normalized = 1 - (0.7 * (values - min(values)) / (max(values) - min(values) + 1e-10)); 
            changes = [1; find(diff(values) ~= 0) + 1; length(values)]; % Indices of change

            % Plot segments
            for j = 1:length(changes)-1
                x_seg = x_ax(changes(j):changes(j+1));
                y_seg = y_offsets(i) * ones(size(x_seg));
                col = colors(i, :) *normalized(changes(j)); % Adjust color intensity
                plot(x_seg, y_seg, 'Color', col, 'LineWidth', 10);
                formatted_text = sprintf('%s [ %g, %g ] [ %s ]', prop_name{i}, min(values), max(values), unities{i});
                text(0.01, y_offsets(i) + 0.1, formatted_text);
            end
        end
        axis([0 max(x_ax) (min(y_offsets)-0.2) (max(y_offsets)+0.2)])
        xticklabels([])
        yticklabels([])
        legend ('off');
        hold off
        % -----------------------------------------------------------------
        % end figure         
    end

    %% MULTILAYER #####                                                    
    function self = multilayer(self,inp)

        % [radius, rho, E, nu, inital node, end node]  
        % fields = {'r','rho','E','nu','no_i','no_f'};
        % for f = 1:numel(fields)
        %         prop_name = fields{f};
        %         masu =  vertcat(self.data(inp).layer{:});
        %     self.layer(inp).data.(prop_name) = masu(:,1);
        % end

        fields = {'no _i','no_f'};
        % for f = 1:numel(fields)
        %     prop_name = fields{f};
        %     masu =  vertcat(self.data(inp).layer{:});
        %     self.layer(inp).data.(prop_name) = masu(:,1);
        % end

            masu =  vertcat(self.data(inp).(fields){:});
        self.layer(inp).data.r = masu(:,1);
        self.layer(inp).data.rho = self.data(inp).layer{:}(:,2);
        self.layer(inp).data.E = self.data(inp).layer{:}(:,3);
        self.layer(inp).data.nu = self.data(inp).layer{:}(:,4);
        self.layer(inp).data.no_i = self.data(inp).layer{:}(:,5);
        self.layer(inp).data.no_f = self.data(inp).layer{:}(:,6);

        % % reorganizing values for each element ____________________________
        % fields = {'r','rho','E','nu'};
        % 
        % for f = 1:numel(fields)
        %     % pre-allocate new_values
        %     new_values = zeros(self.data(inp).n_el, size(self.data(inp).layer,1));
        % 
        %     for j = 1:size(self.data(inp).layer,1)
        % 
        %         % Extrai valores e nós para os outros campos
        %         prop_name = fields{f};
        %         values = self.layer(inp).data.(prop_name)(j);
        % 
        %         % setting the values node intervals
        %         intervals = [self.layer(inp).data.no_f(j); self.layer(inp).data.no_i(j)];
        % 
        %         % fill new_values with the values to its respective nodes
        %         for i = 1:numel(values) %
        %             new_values(intervals(i+1):intervals(i),j) = values(i);
        %         end
        % 
        %         % store the vector on its due field
        %         self.layer(inp).(fields{f}) = new_values;
        %     end
        % end
        % clear idx fields f new_values intervals counts values total_size i col_idx;
        % % _________________________________________________________________
                % reorganizing values for each element ____________________________
        fields = {'r','rho','E','nu'};

        for f = 1:numel(fields)
            % pre-allocate new_values
            new_values = zeros(self.data(inp).n_el, size(self.data(inp).layer,1));

            % for each 
            for j = 1:size(self.data(inp).layer,1)

                % Extrai valores e nós para os outros campos
                prop_name = fields{f};
                values = self.layer(inp).data.(prop_name)(j);

                % setting the values node intervals
                intervals = [self.layer(inp).data.no_f(j); self.layer(inp).data.no_i(j)];

                % fill new_values with the values to its respective nodes
                for i = 1:numel(values) %
                    new_values(intervals(i+1):intervals(i),j) = values(i);
                end

                % store the vector on its due field
                self.layer(inp).(fields{f}) = new_values;
            end
        end
        clear idx fields f new_values intervals counts values total_size i col_idx;
        % _________________________________________________________________
        %         % reorganizing values for each element ____________________________
        % fields = {'d1', 'd2','d3','nu'};% d3 is the form
        % 
        % for f = 1:numel(fields)
        %    if startsWith(fields{f}, 'd')
        %         % get if its d1 d2 or d3
        %         col_idx = str2double(fields{f}(2)); 
        %         % values for each node set
        %         masu = vertcat(self.data(input).d{:});
        %         % nodes where its is applied
        %         [values, intervals] = deal(masu(:, col_idx), masu(:, end-1:end));
        %    else
        %         masu = vertcat(self.data(input).(fields{f}){:});
        %         [values, intervals] = deal(masu(:, 1), masu(:, end-1:end));
        %     end
        % 
        %     % calculate the size of the vector
        %     new_values = zeros(self.data(input).n_el, 1); % pre-allocate new_values
        % 
        %     % fill new_values with the values to its respective nodes
        %     for i = 1:size(intervals,1)
        %         new_values(intervals(i,1) : intervals(i,2)) = values(i);
        %     end
        % 
        %     % store the vector on its due field
        %     self.Beam(input).(fields{f}) = new_values;
        % end
        % clear idx fields f new_values intervals values i col_idx;
        % % _________________________________________________________________

        % rewritting ___________________________________________________
        form = self.Beam(inp).d3;
        db1 = self.Beam(inp).d1;
        db2 = self.Beam(inp).d2;
        r = self.layer(inp).r;
        % pre allocating 
        self.layer(inp).dl1 = zeros(self.data(inp).n_el, size(r, 2));
        self.layer(inp).dl2 = zeros(self.data(inp).n_el, size(r, 2));
        % each layers dimensions
        for j = 1:size(r, 2)
            for i = 1:self.data(inp).n_el
                % Selecionar valor inicial ou acumulado de dl1
                if j == 1
                    self.layer(inp).dl1(i, j) = (form(i) == 1 || form(i) == 2) * db1(i) + ...
                        (form(i) == 3) * (db1(i) + db2(i)) + ...
                        (form(i) == 4) * db2(i);
                else
                    self.layer(inp).dl1(i, j) = self.layer(inp).dl2(i, j-1);
                end

                % Calcular dl2 com base em dl1 e form(i)
                self.layer(inp).dl2(i, j) = self.layer(inp).dl1(i, j) + r(i, j) * ...
                    (form(i) == 1 || form(i) == 3) + 2 * r(i, j) * (form(i) == 2 || form(i) == 4);
            end
        end
        % _________________________________________________________________
        
        % _________________________________________________________________
        for el = 1:self.mesh(inp).n_el

            % rewritting the parameters for convenience
            d2 = self.layer(inp).dl1(el);
            d1 = self.layer(inp).dl2(el);
            nu = self.Beam(inp).nu(el);

            if form(el) == 1 || form(el) == 3 % pipe 
                A = pi*(d1^2 - d2^2);
                I = (pi/2)*(d1^4 - d2^4);
                k = (2 *(1+nu))/(4 + 3*nu);

            elseif form(el) == 2 || form(el) == 4 % box
                A = d1^2 - d2^2;
                I = ((d1^4)/12) - ((d2^4)/12);
                k = (20*(1 + nu))/(48 + 39*nu);

            end
            % ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨ ¨
            % -------------------------------------------------------------
                                            self.layer(inp).A(el,1) = A;
                                            self.layer(inp).I(el,1) = I;
                                            self.layer(inp).k(el,1) = k;
        end
        % _________________________________________________________________

    end % end multilayer

    end % static methods
    methods (Access = private)
        % function result = cellcatcol(field, input, column)
        %     % only for convenience
        %     res = vertcat(self.data(input).(field){:});
        %     result = res(:,column);
        % end
    end

    % ============================== END ================================ %
end % classdef