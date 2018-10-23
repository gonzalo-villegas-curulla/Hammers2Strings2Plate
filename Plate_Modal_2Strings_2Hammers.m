%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Gonzalo Villegas Curulla.
%
%      0. Preamble.
%      Collision 2 hammers - 2 strings.
%      Stiff strings (EB model). Simply supported Boundary Conditionws (BC).
%      Simple spring connection (K stiffness parameter) strings-to-plate.
%      Squared isotropic modal plate. Simply Supported Boundary Conditions (BC).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all, clc;

% Add path
% cd /Users/gonzalo/Documents/GitHub/Hammers2Strings2Plate
% addpath(genpath(pwd));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.Custom parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% EDIT HERE %%%%%%%%

% ---------------------------------------
% I/O

SR   = 44100.0;              % Sample rate (Hz)
Tf   = 3.5;                  % Duration simulation (s)
win_dur = 0.01;              % Fadeout time waveform (s)

k    = 1/SR;                 % time step (s)
Nf   = floor(SR*Tf);         % Time steps for simulation
% ---------------------------------------

% ---------------------------------------
% Flags
flag_sound        = 0;       % Listen to plate output
flag_collision    = 0;       % Visalise collision hammer string in new figure
refresh_collision = 10;      % Set frame rate refresh for collision plot
flag_energy       = 1;       % Hammer/string/coupling/plate/losses/ normalised total
flag_spectrum     = 1;       % FFT string and plate
flag_waveform     = 1;       % Waveform plate out
NumIt             = 20;      % Iterations Newton-Raphson solver (Hammer collision)
% ---------------------------------------

% ---------------------------------------
% Hammer1 setups
M_hammer     = 0.010;             % Hammer mass (kg)
u0_hammer    = -0.04;             % Initial displacement (m)
v0_hammer    = 80;                % Initial velocity (m/s)
xin_hammer   = 0.12;              % Hammer-in point (m) or Chabassier parameter
K_collision  = 4.5e9;             % Collision stiffness (Pa)
alpha        = 2.5;               % Collision exponent
% ---------------------------------------
% ---------------------------------------
% Hammer2 setups
M_hammer2     = 0.030;             % Hammer mass (kg)
u0_hammer2    = -0.04;             % Initial displacement (m)
v0_hammer2    = 100;               % Initial velocity (m/s)
xin_hammer2   = 0.12;              % Hammer-in point (m) or Chabassier parameter
K_collision2  = 4.5e9;             % Collision stiffness (Pa)
alpha2        = 2.5;               % Collision exponent
% ---------------------------------------

% ---------------------------------------
% STRING_1.  Setups

% String.Params.f0 = 65.56; %%% (Hz)
[String] = getParamsChabassier(1,'unwrapped');

E_st   = 2e11;                        % Young's modulus (Pa)
L_st   = String.Params.L;             % Length (m)
r_st   = String.Params.d/2/1000;      % Diameter(mm)->(m)
rho_st = String.Params.rho;           % Material density (kg/m^3)
T_st   = String.Params.T;             % Longitudinal tension (N)
xin_hammer = String.Params.x0;        % Collision point with hammer (m)
op_st      = 0.29;                    % Readout point (m)
sig0_st    = 0.2;                     % Manually set sigma decays
sig1_st    = 0.3;

% Connection_1 with plate
K_connection = 1e5;       % Stiffness1 coupling with plate (N)
% ---------------------------------------

% ---------------------------------------
% STRING_2.  Params

rho_st2 = 7860;       % Material density (steel string) (Kg/m^3)
T_st2   = 4188.2;     % Tension (N);
E_st2   = 2e11;       % Young's modulus (Pa)
r_st2   = 5e-4;       % Radius (m)
L_st2   = 1.32;       % Length (m)
op_st2  = 0.12;       % Read-out audio point (m)
sig0_st2  = 0.3;      % Manually set sigma decays
sig1_st2  = 0.5;

% Connection_2
K_connection2 = 1e5;  % Stiffness2 coupling with plate (N)
% ---------------------------------------

% ---------------------------------------
% PLATE.
Lx   = 2;                    % Plate X (m)
Ly   = 2;                    % Plate Y (m)
H    = 0.01;                 % Plate thickness (m)
mat_inx = 3;                 % Material selection index (1 = steel, 2 = gold, 3 = silver, 4 = aluminium)
maxFreqPlate = 1e4;          % Max frequency mode in plate (Hz)

% COUPLING POSITIONS.
coupling_point  = [0.2*Lx,0.2*Ly];   % Input point from string1 (w1w)
coupling_point2 = [0.7*Lx,0.7*Ly];   % With string2 (w2w)
out_point_L_pl  = [0.01*Lx, 0.5*Ly]; % Left channel out plate (pickup 1)
out_point_R_pl  = [0.7*Lx,0.7*Ly];   % Right channel out plate (pickup2)


% PLATE. Physical parameters.
mat_tab = [2e11  7860 0.30;  % Steel
    7.9e10 19300 0.40;       % Gold
    8.3e10 10490 0.37;       % Silver
    1.16e11 4500 0.30;       % Titanium
    7.0e10  2700 0.35];      % Aluminium

% PLATE. Losses per octaves.
mt =  1;
T60_01 = 5.50*mt;
T60_02 = 4.20*mt;
T60_03 = 3.89*mt;
T60_04 = 3.75*mt;
T60_05 = 3.24*mt;
T60_06 = 3.29*mt;
T60_07 = 3.15*mt;
T60_08 = 3.10*mt;
% ---------------------------------------

%%%%%%%% END OF EDITS %%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Derived parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------
% 2.1.1.PLATE. Physical parameters.
E_p        =  mat_tab(mat_inx,1);           % Young's modulues (Pa)
rho_p      =  mat_tab(mat_inx,2);           % Material density (kg/m3)
poisson_p  =  mat_tab(mat_inx,3);           % Poisson coefficient
D          =  E_p*H^3/12/(1-poisson_p^2);   % Flexural rigidity (Kg*s^2/m^2)

% ---------------------------------------
% 2.1.2.PLATE. Modes computation.

maxOmega = 2*pi*maxFreqPlate;       % Omega max (rad/s)
maxL     = max([Lx,Ly]);            % Max dimension length of plate (m)
DD       = floor(maxL/pi*sqrt((sqrt(4*maxOmega^2*rho_p*H*D))/(2*D))); % Max mode
omega_vec       = zeros(DD*DD,3);          % omegas vector
ind = 0;                            % Counter

for rx = 1 : DD                     % Modes index X
    for ry = 1 : DD                 % Modes index Y
        ind = ind + 1;
        omega_r_2 = D/rho_p/H*((rx*pi/Lx)^2+(ry*pi/Ly)^2)^2;
        omega_r = sqrt(omega_r_2);
        omega_vec(ind,:) = [omega_r,rx,ry];
    end
end

omega_vec = sortrows(omega_vec);            % Sort in ascending order frequency. Keep respective 'rx' and 'ry' indices
for n = 1:DD*DD
    if (omega_vec(n,1) > maxOmega)
        break
    end
end
DIM = n-1;                     % Pick max index mode ('lined up')
omega_vec  = omega_vec(1:DIM,:);             % Shorten omega vector to only those freqs within maxOmega = 14000*2*pi
ff  = omega_vec(:,1)/2/pi;            % All modal freqs in Hz sorted
disp('Plate modes computed.');

% ---------------------------------------
% 2.1.3.PLATE. Losses per octaves.

c = zeros(DIM,1);
for m = 1 : DIM
    if ff(m) < 160
        c(m) = 1/(T60_01);
    elseif ff(m) >= 160 && ff(m) < 320
        c(m) = 1/(T60_02);
    elseif ff(m) >= 320 && ff(m) < 640
        c(m) = 1/(T60_03);
    elseif ff(m) >= 640 && ff(m) < 1250
        c(m) = 1/(T60_04);
    elseif ff(m) >= 1250 && ff(m) < 2500
        c(m) = 1/(T60_05);
    elseif ff(m) >= 2500 && ff(m) < 5000
        c(m) = 1/(T60_06);
    elseif ff(m) >= 5000 && ff(m) < 10000
        c(m) = 1/(T60_07);
    elseif ff(m) >= 10000
        c(m) = 1/(T60_08);
    end
end
cf = 6*log(10);
c  = c*cf;
clear m
disp('Plate losses computed.');
% c = zeros(DIM,1);


% ---------------------------------------
% 2.1.5 STABILITY CONDITIONS. Plate.
den_stab = 0.25*max(max(omega_vec)) + 4*DIM*K_connection / (Lx*Ly*H*rho_p);
if gt(k^2, 1/den_stab)
  error('Condition plate on k^2 violated. Check dimensions, maxFreqPlate and SR');
end

maxK_connectionAllowed = 0.25*Lx*Ly*rho_p*H*(SR^2 - 0.25*max(max(omega_vec))^2);
maxOmegaAllowed = 2*sqrt(SR^2 - 4*DIM*K_connection/(Lx*Ly*H*rho_p));

% ---------------------------------------

% ---------------------------------------
% 2.2.1.STRING 1. Grid and stability
A_st     = pi*(r_st^2);           % Cross-section area (m^2)
I_st     = 0.25*pi*(r_st^4);      % Radius of gyration or area moment inertia

% Newton-Raphson solver for String1 grid-spacing
h0  = 0.1;           % Initial guess
dh  = 1;             % Diff consecutive guesses
thr = 1e-10;         % Threshold between two consecutive iterations.
while (dh > thr)
  f_h = (rho_st*A_st)*(h0^4) -K_connection*(k^2)*(h0^3) - T_st*(k^2)*(h0^2) - 4*E_st*I_st*(k^2);
  fp_h = 4*(rho_st*A_st)*(h0^3) - 3*K_connection*(k^2)*(h0^2) - 2*T_st*(k^2)*h0;
  h1 = h0 - f_h / fp_h;
  dh = abs(h1-h0);
  h0 = h1;
end
hmin = h0;

% Grid
N      = floor(L_st/hmin);                % N segments
h      = L_st/N;                          % Grid spacing (m)
li     = round(xin_hammer/h);            % Collision point
g      = zeros(N,1);                     % Distribution collision
g(li)  = 1/h;
lo_st  = round(op_st/h);             % Readout point string
disp('Stability condition on String1 computed.');
% ---------------------------------------

% ---------------------------------------
% 2.2.1.STRING 2. Grid and stability
A_st2     = pi*(r_st2^2);                     % Cross-section area (m^2)
I_st2     = 0.25*pi*(r_st2^4);                % Radius of gyration or area moment inertia

% Newton-Raphson solver
h0  = 0.1;           % Initial guess
dh  = 1;             % Diff consecutive guesses
thr = 1e-10;        % Threshold between two consecutive iterations.
cnt = 2;
while (dh > thr)
  f_h = (rho_st2*A_st2)*(h0^4) -K_connection2*(k^2)*(h0^3) - T_st2*(k^2)*(h0^2) - 4*E_st2*I_st2*(k^2);
  fp_h = 4*(rho_st2*A_st2)*(h0^3) - 3*K_connection2*(k^2)*(h0^2) - 2*T_st2*(k^2)*h0;
  h1 = h0 - f_h / fp_h;
  dh = abs(h1-h0);
  h0 = h1;
  cnt = cnt+1;

end
hmin2 = h0;

% Grid
N2      = floor(L_st2/hmin2);              % N2 segments
h2      = L_st2/N2;                        % Grid spacing2 (m)
li2     = round(xin_hammer2/h2);           % Collision point
g2      = zeros(N2,1);                     % Distribution collision
g2(li2) = 1/h2;
lo_st2  = round(op_st2/h2);                % Readout point string2
disp('Stability condition on String2 computed.');
% ---------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Error check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------
% Error check on stability condition
if omega_vec(end,1) > 2/k
    error('Condition on omega max violated');
end

% Plate sizes and points
if out_point_L_pl(1)>Lx | out_point_R_pl(1) > Lx | out_point_L_pl(2) > Ly | out_point_R_pl(2) >Ly | coupling_point(1)>Lx | coupling_point(2)>Ly
    warning('Check dimensions and points on plate');
end
% ---------------------------------------

% Sizes Lx, Ly, H, L_st, L_st2, readoutpoints,




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Initialisations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------
% 4.1.1 HAMMER
m = (k^2)*((h*g'*g)/((rho_st*A_st)*(1 + sig0_st*k)) +1/M_hammer); % Scheme factor

U2_hammer = u0_hammer;
U1_hammer = U2_hammer + k*v0_hammer;
U_hammer  = 0;

collision_history = zeros(Nf,1);
% ---------------------------------------


% ---------------------------------------
% 4.1.2 HAMMER2
mm = (k^2)*((h2*g2'*g2)/((rho_st2*A_st2)*(1 + sig0_st2*k)) +1/M_hammer2); % Scheme factor

U2_hammer2 = u0_hammer2;
U1_hammer2 = U2_hammer2 + k*v0_hammer2;
U_hammer2  = 0;
% ---------------------------------------



% ---------------------------------------
% 4.2.1 %%%%% STRING 1 %%%%% (hammer collision)
out_st = zeros(Nf,1);     % Out audio
w1w2   = zeros(N, 1);
w1w1   = zeros(N,1);
w1w    = zeros(N,1);

xax_st = (1:N)'*h;         % Axis along string in physical units (m)

% 4.2.2. Hammer-String1 Interpenetrations
eta1 = U1_hammer - w1w1(li);  % Interpenetration time step 'n' (m)
eta2 = U2_hammer - w1w2(li);
eta  = 0;

% 4.2.3. STRING 1. DIFFERENCE MATRICES.
v          = ones(N,1);
Dxplus     = 1/h*spdiags([-v, v],-1:0,N,N); % size: NxN

Dxx        = toeplitz([-2/(h^2) 1/(h^2) zeros(1,N-2)]);  % size: NxN
Dxx(end,:) = 0;       % Set according to BCond dxxW_N = 0
Dxx        = sparse(Dxx);

Dxxxx        = toeplitz([6/(h^4) -4/(h^4) 1/(h^4) zeros(1,N-3)]);    % Size: NxN
Dxxxx(1,1)   = 5/(h^4);                                              % Row 1
Dxxxx(end-1,end-1) = 5/(h^4);   Dxxxx(end-1, end) = -2/(h^4);        % Row N-1
Dxxxx(end,:) = 0;                                                    % Row N
Dxxxx(end,end-2 : end) = [1 , -(T_st*h^2/E_st/I_st + 2) , ( (h^3*K_connection + T_st*h^2)/E_st/I_st + 1) ]./(h^4); % BC coupling
Dxxxx        = sparse(Dxxxx);

% 4.2.4 MATRIX FORM Scheme
AA = (1 + k*sig0_st)*eye(N);
AA = sparse(AA);
BB = T_st*(k^2)/(rho_st*A_st)*Dxx - E_st*I_st*(k^2)/(rho_st*A_st)*Dxxxx + 2*(eye(N)+sig1_st*k*Dxx);
BB = sparse(BB);
CC = (1 - sig0_st*k)*eye(N) + 2*sig1_st*k*Dxx;
CC = sparse(CC); % Remember taking '-' sign into account
% ---------------------------------------



% ---------------------------------------
% 4.3. COUPLING 1.
u_xc  = 0;  % Plate displacement @ coupling point next time step
u_xc1 = 0;
u_xc2 = 0;
fvec      = zeros(N,1);                         % Against string 1
fvec(end) = k^2* K_connection / h / (rho_st*A_st) ; % Force value factor @ coupling point 1
fin       = 0;
% ---------------------------------------


% ---------------------------------------
% 4.4  %%%%% String 2 %%%%%
out_st2 = zeros(Nf,1);     % Out audio
w2w     = zeros(N2,1);
w2w1    = zeros(N2,1);
w2w2    = zeros(N2,1);

% 4.4.2. Hammer2-String2 Interpenetrations
etaa1 = U1_hammer2 - w2w1(li2);  % Interpenetration time step 'n' (m)
etaa2 = U2_hammer2 - w2w2(li2);
etaa  = 0;

% 4.4.3. STRING 2. DIFFERENCE MATRICES .
v2          = ones(N2,1);
Dxplus2     = 1/h2*spdiags([-v2, v2],-1:0,N2,N2); % size: N2xN2

Dxx2        = toeplitz([-2/(h2^2) 1/(h2^2) zeros(1,N2-2)]);  % size: N2xN2
Dxx2(end,:) = 0;       % Set according to BCond dxx2W_N2 = 0
Dxx2        = sparse(Dxx2);

Dxxxx2        = toeplitz([6/(h2^4) -4/(h2^4) 1/(h2^4) zeros(1,N2-3)]);    % Size: N2xN2
Dxxxx2(1,1)   = 5/(h2^4);                                              % Row 1
Dxxxx2(end-1,end-1) = 5/(h2^4);   Dxxxx2(end-1, end) = -2/(h2^4);        % Row N2-1
Dxxxx2(end,:) = 0;                                                    % Row N2
Dxxxx2(end,end-2 : end) = [1 , -(T_st2*h2^2/E_st2/I_st2 + 2) , ( (h2^3*K_connection2 + T_st2*h2^2)/E_st2/I_st2 + 1) ]./(h2^4); % BC coupling
Dxxxx2        = sparse(Dxxxx2);

% 4.4.4 MATRIX FORM SCHEME
AA2 = (1 + k*sig0_st2)*eye(N2);
AA2 = sparse(AA2);
BB2 = T_st2*(k^2)/(rho_st2*A_st2)*Dxx2 - E_st2*I_st2*(k^2)/(rho_st2*A_st2)*Dxxxx2 + 2*(eye(N2)+sig1_st2*k*Dxx2);
BB2 = sparse(BB2);
CC2 = (1 - sig0_st2*k)*eye(N2) + 2*sig1_st2*k*Dxx2;
CC2 = sparse(CC2); % Remember taking '-' sign into account
% ---------------------------------------


% ---------------------------------------
% 4.5 Coupling 2
u2_xc  = 0;  % Plate displacement @ coupling point next time step
u2_xc1 = 0;
u2_xc2 = 0;

fvec2      = zeros(N2,1);  % Against string 2.
fvec2(end) = k^2* K_connection2 / h2 / (rho_st2*A_st2) ; % Force value factor @ coupling point 2
fin2       = 0;
% ---------------------------------------



% ---------------------------------------
% 4.6. %%%%% PLATE %%%%%
% Updates.

fac = 1./(1/k^2 + c/k);
G1  = fac.*(2/k^2 - omega_vec(:,1).^2);
G2  = fac.*(c/k - 1/k^2);
P   = (4/Lx/Ly) * fac .* sin(coupling_point(1)*pi*omega_vec(:,2)/Lx) .* sin(coupling_point(2)*pi*omega_vec(:,3)/Ly);
P2   = (4/Lx/Ly) * fac .* sin(coupling_point2(1)*pi*omega_vec(:,2)/Lx) .* sin(coupling_point2(2)*pi*omega_vec(:,3)/Ly);


proj_vec = ( sin(coupling_point(1)*pi*omega_vec(:,2)/Lx).*sin(coupling_point(2)*pi*omega_vec(:,3)/Ly) )';
proj_vec2 = ( sin(coupling_point2(1)*pi*omega_vec(:,2)/Lx).*sin(coupling_point2(2)*pi*omega_vec(:,3)/Ly) )';

rpL =  sin(omega_vec(:,2)*pi*out_point_L_pl(1)/Lx) .* sin(omega_vec(:,3)*pi*out_point_L_pl(2)/Ly); % Modal shapes in x_out
rpL = rpL.';
rpR =  sin(omega_vec(:,2)*pi*out_point_R_pl(1)/Lx) .* sin(omega_vec(:,3)*pi*out_point_R_pl(2)/Ly); % Modal shapes in x_out
rpR = rpR.';

% Modal states at time steps nth and 'n-1'th
q2  = zeros(DIM,1);
q1  = zeros(DIM,1);
q   = zeros(DIM,1);
% ---------------------------------------

% ---------------------------------------
% 4.7.1 Energy vectors

H_string1  = zeros(Nf,1);      % Energy string potential (J)
H_string2  = zeros(Nf,1);      % Energy string kinetic (J)
H_hammer   = zeros(Nf,1);      % Energy hammer (J)
Q          = zeros(Nf,1);      % Losses string (J)
H_coupling = zeros(Nf,1);      % Energy coupling connection (J)
H_plate    = zeros(Nf,1);      % Energy plate (J)
Qp         = zeros(Nf,1);      % Losses plate (J)
AcumQ      = zeros(Nf,1);      % Integrated losses vector


% 4.7.2 Energy coefficients
Ecoef1 = 0.5*(rho_st*A_st)*h*(SR^2); % String1
Ecoef2 = 0.5*T_st*h;
Ecoef3 = 0.5*E_st*I_st*h;
Ecoef4 = 0.5*sig1_st*(rho_st*A_st)*h*SR;

Ecoef13 = 0.5*(rho_st2*A_st2)*h2*(SR^2); % String2
Ecoef14 = 0.5*T_st2*h2;
Ecoef15 = 0.5*E_st2*I_st2*h2;
Ecoef16 = 0.5*sig1_st2*(rho_st2*A_st2)*h2*SR;

Ecoef17 = 0.5*(rho_st2*A_st2)*sig0_st2*h2*(SR^2); % Losses string 2
Ecoef18 = 0.5*(rho_st2*A_st2)*sig1_st2*h2*(SR^2); % Losses string 1


Ecoef5  = 0.5*M_hammer*(SR^2); % Hammer
Ecoef6  = 0.5*(rho_st*A_st)*sig0_st*h*(SR^2); % Losses string 1
Ecoef7  = 0.5*(rho_st*A_st)*sig1_st*h*(SR^2); % Losses string 1
Ecoef8  = 0.25*Lx*Ly; %Plate
Ecoef9  = Ecoef8*0.5*(rho_p*H)*(SR^2);
Ecoef10 = Ecoef8*0.5*(rho_p*H).*(omega_vec(:,1).^2).';
Ecoef11 = 0.5*K_connection; % Coupling spring
Ecoef12 = Ecoef8*0.5*rho_p*H*(SR^2)*(c'); % Losses plate
% ---------------------------------------

% ---------------------------------------
% 4.8. Audio Readout vectors
outL  = zeros(Nf,1);
outR  = zeros(Nf,1);
% ---------------------------------------

% ---------------------------------------
% 4.9 Visualisations, axes, figures
if flag_collision
  fcoll = figure('Name', 'Collision', 'Units','Normalized','Position',[0.125 .125 .75 .75]);
  axcoll = axes(fcoll);
  axis_coll = [0 L_st -0.15 0.15];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Main Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
for n = 1 : Nf

    % Forces Coupling 1 ---------------------------------------
    f = K_connection * (w1w1(end) - u_xc1);
    fin =  f/rho_p/H;                     % Against plate
    f_connection = u_xc1*fvec;            % Against string

    % Forces Coupling 2 ---------------------------------------
    f2 = K_connection2 * (w2w1(end) - u2_xc1);
    fin2 =  f2/rho_p/H;                     % Against plate
    f_connection2 = u2_xc1*fvec2;

    % Collision 1---------------------------------------
    a = eta2;
    nu1 = (2/k/k)*(w1w1-w1w2) + (T_st/(rho_st*A_st))*Dxx*w1w1 - (E_st*I_st/(rho_st*A_st))*Dxxxx*w1w1 + (2*sig1_st/k)*Dxx*(w1w1-w1w2);
    b = -2*(U1_hammer - U2_hammer) + k^2*nu1(li)/(1+sig0_st*k);
    r0 = -b;
    if (r0+a) < 0 && a < 0 % No interpenetration
      r = r0;
    else
      for nwt = 1 : NumIt % Newton solver
        G_r0 = r0 + (m/r0)*(phi(K_collision,alpha,r0+a) - phi(K_collision,alpha,a)) + b; % G(r0) polynomial
        Gp_r0 = 1 + (m/r0/r0)*(r0*phiPrimed(K_collision,alpha,r0+a)-phi(K_collision,alpha,r0+a) +...
         phi(K_collision,alpha,a)); % G'(r0)
        r0 = r0 - G_r0 / Gp_r0;
      end
      r = r0;
      % n
    end
    eta = r + eta2;
    f_hammer = (phi(K_collision,alpha,eta) - phi(K_collision,alpha,eta2)) / (eta - eta2);
    collision_history(n) = f_hammer;
    U_hammer = 2*U1_hammer - U2_hammer - f_hammer*(k^2)/M_hammer;

    % Collision 2---------------------------------------
    aa = etaa2;
    nuu1 = (2/k/k)*(w2w1-w2w2) + (T_st2/(rho_st2*A_st2))*Dxx2*w2w1 - (E_st2*I_st2/(rho_st2*A_st2))*Dxxxx2*w2w1 + (2*sig1_st2/k)*Dxx2*(w2w1-w2w2);
    bb = -2*(U1_hammer2 - U2_hammer2) + k^2*nuu1(li2)/(1+sig0_st2*k);
    r00 = -bb;
    if (r00+aa) < 0 && aa < 0 % No interpenetration
      rr = r00;
    else
      for nwt = 1 : NumIt % Newton solver
        G_r00 = r00 + (mm/r00)*(phi(K_collision2,alpha2,r00+aa) - phi(K_collision2,alpha2,aa)) + bb; % G(r02) polynomial
        Gp_r00 = 1 + (mm/r00/r00)*(r00*phiPrimed(K_collision2,alpha2,r00+aa)-phi(K_collision2,alpha2,r00+aa) +...
         phi(K_collision2,alpha2,aa)); % G'(r0)
        r00 = r00 - G_r00 / Gp_r00;
      end
      rr = r00;
      % n
    end
    etaa = rr + etaa2;
    f_hammer2 = (phi(K_collision2,alpha2,eta2) - phi(K_collision2,alpha2,etaa2)) / (etaa - etaa2);
    % collision_history2(n) = f_hammer2;
    U_hammer2 = 2*U1_hammer2 - U2_hammer2 - f_hammer2*(k^2)/M_hammer2;

    % String 1% ---------------------------------------
    w1w = AA \ (BB*w1w1 - CC*w1w2 + f_connection + (k^2)/(rho_st*A_st)*g*f_hammer);
    % String 2% ---------------------------------------
    w2w = AA2 \ (BB2*w2w1 - CC2*w2w2 + f_connection2 + (k^2)/(rho_st2*A_st2)*g2*f_hammer2);

    % Plate % ---------------------------------------
    q = G1.*q1  + G2.*q2 + P*fin + P2*fin2;


    % Read plate displacement at connection points
    u_xc = proj_vec*q;
    u2_xc = proj_vec2*q;


    % Visualize hammer and string ---------------------------------------

    if flag_collision & ~mod(n,refresh_collision) %& (abs(U_hammer) < max(abs(axis_coll(3:4))))
        fcoll;
        plot(axcoll,xax_st, w1w,'k', 'linewidth',1.1);
        hold on
        scatter(axcoll,xin_hammer,U_hammer,'MarkerFaceColor','k','SizeData',80);
        hold off
        axis(axis_coll)
        drawnow
    end

    % Energy % ---------------------------------------

    if flag_energy

        H_string1(n) = Ecoef1*(w1w1-w1w2)'*(w1w1-w1w2) + Ecoef2*(Dxplus*w1w1)'*(Dxplus*w1w2) + ...
            Ecoef3*(Dxx*w1w1)'*(Dxx*w1w2) - Ecoef4*(Dxplus*(w1w1-w1w2))'*(Dxplus*(w1w1-w1w2));

        H_string2(n) = Ecoef13*(w2w1-w2w2)'*(w2w1-w2w2) + Ecoef14*(Dxplus2*w2w1)'*(Dxplus2*w2w2) + ...
            Ecoef15*(Dxx2*w2w1)'*(Dxx2*w2w2) - Ecoef16*(Dxplus2*(w2w1-w2w2))'*(Dxplus2*(w2w1-w2w2));


        H_plate(n)    =  Ecoef9*(q1-q2)'*(q1-q2) + Ecoef10* (q1.*q2)  ;
        H_coupling(n) = Ecoef11*(w1w1(end) - u_xc1)*(w1w2(end) - u_xc2);
        H_hammer(n)   = Ecoef5*(U1_hammer-U2_hammer)^2 + 0.5*(phi(K_collision,alpha,eta1) +...
              phi(K_collision,alpha,eta2));

        % Strings losses
        Q(n) = Ecoef6*(w1w1-w1w2)'*(w1w1-w1w2) + Ecoef7*(Dxplus*(w1w1-w1w2))'*(Dxplus*(w1w1-w1w2)) +...
               Ecoef17*(w2w1-w2w2)'*(w2w1-w2w2) + Ecoef18*(Dxplus2*(w2w1-w2w2))'*(Dxplus2*(w2w1-w2w2));

        % Plate loss
        Qp(n) = Ecoef12*((q-q2).^2);

        % Integrated losses
        if n == 1
          AcumQ(n) = -k*(Q(n)+Qp(n));
        else
          AcumQ(n) = AcumQ(n-1) - k*(Q(n)+Qp(n));
          % Inexpensive version of AcumQ(n) = sum(Q(1:n));
        end

    end

    % Read-outs % ---------------------------------------
    outL(n) = rpL*q;
    outR(n) = rpR*q;
    out_st(n) = w1w(lo_st);
    out_st2(n) = w2w(lo_st2);

    % Feed next time step % ---------------------------------------
    w1w2 = w1w1; w1w1 = w1w;    % String 1
    w2w2 = w2w1; w2w1 = w2w;    % String 2

    q2 = q1; q1 = q;                  % Modes plate

    u_xc2 = u_xc1; u_xc1 = u_xc;      % Plate displacement @ connection 1
    u2_xc2 = u2_xc1; u2_xc1 = u2_xc;  % Plate displacement @ connection 2

    U2_hammer =  U1_hammer; U1_hammer = U_hammer;     % Hammer displacement
    eta2 = eta1; eta1 = eta;          % Interpentrations hammer-string1


end % End main loop
toc

% Normalise outputs

outL    = outL / max(abs(outL)); % Plate
outR    = outR / max(abs(outR));
out     = [outL, outR];
out_st  = out_st/max(abs(out_st)); % String
out_st2 = out_st2/max(abs(out_st2)); % String2




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Outs, plots and visualizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tvec = k*[0:Nf-1]';


% 6.1 Audio out
if flag_sound
    fade_N        = floor(win_dur*SR);
    fade_half_win = (0.5 - 0.5 * cos(pi*(fade_N-1:-1:0)/fade_N))';
    out(end-fade_N+1:end,1) = out(end-fade_N+1:end,1).*fade_half_win;
    out(end-fade_N+1:end,2) = out(end-fade_N+1:end,2).*fade_half_win;
    out_st(end-fade_N+1:end) = out_st(end-fade_N+1:end).*fade_half_win;
    out_st(1:fade_N) = out_st(1:fade_N).*flipud(fade_half_win);
    out_st2(end-fade_N+1:end) = out_st2(end-fade_N+1:end).*fade_half_win;

    out(:,1) = [diff(out(:,1)); 0]; out(:,1) = out(:,1)/max(abs(out(:,1)));
    out(:,2) = [diff(out(:,2)); 0]; out(:,2) = out(:,2)/max(abs(out(:,2)));

    audiovec = [out_st out_st; out; out_st2 out_st2];
    soundsc(audiovec,SR);

end
% ---------------------------------------

% ---------------------------------------
% 6.2 Spectra
if flag_spectrum

    OUTL      = 20*log10(abs(fft(outL)));
    df        = SR/length(OUTL);
    f_axis    = (0:df:SR-df)';
    OUTR      = 20*log10(abs(fft(outR)));
    OUTS      = 20*log10(abs(fft(out_st)));
    OUTS2     = 20*log10(abs(fft(out_st2)));
    freq0_st  = sqrt(T_st/rho_st/A_st)/2/L_st;
    maxfplate = max(max(OUTL),max(OUTR));


    fig2 = figure(2);

    axspc1 = axes(fig2,'Units','Normalized','OuterPosition',[0. .66 1 .33]);     % [0.3*freq0_st 1.2*8*freq0_st]
    title3 = "Spectrum plate H = %.3f (m)";
    title3 = compose(title3,H);
    semilogx(axspc1,f_axis, OUTL-maxfplate, f_axis, OUTR-maxfplate);
    xlim(axspc1,[30 SR/2]);
    title(axspc1,title3);
    xlabel(axspc1, 'Frequency (Hz)');
    ylabel(axspc1, 'dBFS');
    hold on
    grid on
    line(axspc1,[freq0_st freq0_st],get(gca, 'ylim'),'linestyle',':','color','r', 'linewidth',1.3);
    hold off


    axspc2 = axes(fig2,'Units','Normalized','OuterPosition',[0. .33 1 .33]);% [0.3*freq0_st 1.2*8*freq0_st]
    title1 = "Spectrum string1 (%.1f Hz) Hammer1-excited";
    title1 = compose(title1, freq0_st);
    semilogx(axspc2,f_axis, OUTS-max(OUTS));
    xlim(axspc2,[30 SR/2]);
    title(axspc2,title1);
    xlabel(axspc2,'Frequency (Hz)');
    ylabel(axspc2,'dBFS');
    hold on
    grid on
    line(axspc2,[freq0_st freq0_st],get(gca, 'ylim'),'linestyle',':','color','r', 'linewidth',1.3);
    hold off


    axspc3 = axes(fig2,'Units','Normalized','OuterPosition',[0. .0 1. .33]);     % [0.3*freq0_st 1.2*8*freq0_st]
    freq0_st2 = 0.5*sqrt(T_st2/rho_st2/A_st2)/L_st2;
    title2 = "Spectrum string2 (%.1f Hz) Hammer2-excited";
    title2 = compose(title2, freq0_st2);
    semilogx(axspc3,f_axis, OUTS2-max(OUTS2));
    xlim(axspc3,[30 SR/2]);
    title(axspc3,title2);
    xlabel(axspc3,'Frequency (Hz)');
    ylabel(axspc3,'dBFS');
    hold on
    grid on
    line(axspc3,[freq0_st2 freq0_st2],get(gca, 'ylim'),'linestyle',':','color','r', 'linewidth',1.3);
    hold off

    set(fig2, 'Units', 'normalized', 'Position', [0.,0.,1.,1.]);

    % fig2.Children


end
% ---------------------------------------


% ---------------------------------------
% 6.3 Waveforms
if flag_waveform
  font_text = 12;

  fig3 = figure('Name', 'Waveforms', 'Units', 'Normalized', 'Outerposition', [.05 .0344 0.4139 0.45]);
  axwv1 = axes(fig3,'Units','Normalized','OuterPosition',[0. .66 1 .33]);
  plot(axwv1, (1:Nf)'*k,out); title(axwv1,'Waveform plate', 'FontSize',font_text);
  xlabel(axwv1,'Time (s)', 'FontSize',font_text);
  ylabel(axwv1,'Displ. (norm)', 'FontSize',font_text);

  axwv2 = axes(fig3,'Units','Normalized','OuterPosition',[0. .33 1 .33]);
  plot(axwv2, (1:Nf)'*k,out_st);
  title(axwv2, 'Waveform string 1', 'FontSize',font_text);
  xlabel(axwv2, 'Time (s)', 'FontSize',font_text);
  ylabel(axwv2, 'Displ. (norm)', 'FontSize',font_text);

  axwv3 = axes(fig3,'Units','Normalized','OuterPosition',[0. .0 1 .33]);
  plot(axwv3,(1:Nf)'*k, out_st2);
  title(axwv3,'Waveform string 2', 'FontSize',font_text);
  xlabel(axwv3,'Time (s)', 'FontSize',font_text);
  ylabel(axwv3,'Displ. (norm)', 'FontSize',font_text);

  % fig3.Children

end
% ---------------------------------------





% ---------------------------------------
% 6.4 Energy
if flag_energy

    H_tot  = H_string1 + H_string2  + H_plate + H_coupling + H_hammer;
    H_norm = 1 - H_tot/max(H_tot(1));
    % H_norm = abs(H_tot - H_tot(1))/H_tot(1);

    coll_start = find(sign(collision_history) == 1, 1, 'first');
    coll_end   = find(sign(collision_history) == 1, 1, 'last');
    coll_dur   = coll_end - coll_start +1;

    fig4 = figure('Name','Normalized energy','Units', 'Normalized', 'Name', 'Energy terms', 'OuterPosition', [0.0500    0.5133    0.4139    0.4611]);

    axenrg1 = axes(fig4,'Units','Normalized','OuterPosition',[0. .66 1 .33], 'FontSize', 14);
    plot(axenrg1,1:Nf,H_norm,'linewidth',1.2); title(axenrg1,'Normalised energy');
    ylabel(axenrg1,'1-H_{tot}/H_{0}');%xlabel('Time (samples)');
    axenrg1.FontSize = 14;

    axenrg2 = axes(fig4,'Units','Normalized','OuterPosition',[0. .33 1 .33], 'FontSize', 14);
    scatter(axenrg2,1:coll_dur,H_norm(coll_start:coll_end),2.);
    title(axenrg2,'During collision');
    xlabel(axenrg2,'Time (samples)');ylabel(axenrg2,'1 - H_{tot}/H_{0}');
    axenrg2.FontSize = 14;


    axenrg3 = axes(fig4,'Units','Normalized','OuterPosition',[0. .0 1 .33], 'FontSize', 14);
    plot(axenrg3,1:Nf, Q);title(axenrg3,'Q losses sig0_{st1} and sig1_{st1}');
    hold on
    plot(axenrg3,1:Nf, Qp);
    axenrg3.FontSize = 14;

    % fig4.Children

    %%%%%%%%%%%%%%%%%%%%%

    fig5 = figure('Units', 'Normalized', 'Name', 'Energy terms', 'OuterPosition', [0.5500    0.0344    0.4333    0.9400]);

    axE1 = axes(fig5,'Units','Normalized','OuterPosition',[0. .83 1 .16], 'FontSize', 14);
    plot(axE1,H_hammer(coll_start:coll_end));
    ylabel(axE1,'H (J)');title(axE1, 'Hammer 1');

    axE2 = axes(fig5,'Units','Normalized','OuterPosition',[0. .66 1 .16], 'FontSize', 14);
    plot(axE2,tvec,H_string1);
    ylabel(axE2,'H (J)'); title(axE2, 'String 1');


    axE3 = axes(fig5,'Units','Normalized','OuterPosition',[0. .5 1 .16], 'FontSize', 14);
    plot(axE3,tvec,H_string2);
    ylabel(axE3,'H (J)');title(axE3, 'String 2');

    axE4 = axes(fig5,'Units','Normalized','OuterPosition',[0. .33 1 .16], 'FontSize', 14);
    plot(axE4,tvec,H_plate);
    ylabel(axE4,'H (J)');title(axE4, 'Plate');


    axE5 = axes(fig5,'Units','Normalized','OuterPosition',[0. .16 1 .16], 'FontSize', 14);
    plot(axE5,tvec, H_coupling);
    ylabel(axE5,'H (J)'); title(axE5, 'Coupling1');

    axE6 = axes(fig5,'Units','Normalized','OuterPosition',[0. .0 1 .16], 'FontSize', 14);
    plot(axE6,collision_history(coll_start:coll_end));
    title(axE6,'Coll.Hist - Interpenetration');

    % fig5.Children


end
% figure

% plot(H_tot - H_tot(1) -AcumQ);
% hold on
% plot(H_tot - H_tot(1)); title('Space Integrated H_{tot} - H_0 + k*Q = 0');
% hold on
% plot(AcumQ*4);

% figure();
% subplot(2,1,1)
% plot(Q); ylabel('String (J)');title('Losses');
% subplot(2,1,2)
% plot(Qp); ylabel('Plate (J)');
