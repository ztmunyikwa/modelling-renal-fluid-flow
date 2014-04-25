function kidney_code_v6
%THIS VERSION IS FOR SUSTAINED OR STEP CHANGES IN PRESSURE
clear all
close all
tic
format long e
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%INPUT PARAMETERS
%change as needed for simulations
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
tinc = 0.1;         % [s]   time increment for calling ODE solver; 
%                           Time vector is in increments of tinc
N = 161;            % [-]   number of discretization points (including
%                           boundaries) along the 1D length of the TAL
tfinal = 300;        % [s]  steady state simulation time for each P value 
C0 = 275e-3;        % [mmol/cm^3] C(0,t); boundary concentration of Cl-
%                           at loop bend
Vmax = 14.5e-6;     % [mmol/cm^2/s] 
%                           maximum transport rate of Cl- from TAL    
Km = 70e-3;         % [mmol/cm^3] 
%                           active solute transport Michaelis constant
p = 1.5e-5;         % [cm/s] trans-epithelial diffusion backleak 
%                           chloride permeability
tau0 = 4;           % [s]   time delay at juxtaglomerular apparatus (JGA)
%                           code sets the minimum tau0 to one time step
td = 1;             % [s]   time constant for diameter ODE
ta = 10;            % [s]   time constant for activation ODE
Pa = 1333.*[100 100 110 120 130 140 150 160 170 180];% [dynes/cm^2] 
%                           incoming pressure vector [constrolstate repeatconstrolstate otherpressures]
constantalpha = 1;  % [-]   1 for constant alpha = 0.2
%                           0 for pressure natriuresis
alpha = 0.2;        % [-]   fraction of single-nephron golumerular filtration rate reaching TAL
%                           this value is overwritten if constantalpha = 0
CLc = 32.32606517560781;
%Turn Myo or TGF on or off. Use 1 for on and 0 for off.
MyoOn = 1;
TGFOn = 1;
%name the .mat output file; saved at the end of the file
outputfilename = ['output_Pvector_tau' num2str(tau0) '_tfinal' num2str(tfinal) 's_td' num2str(td) '_ta' num2str(ta) '.mat'];


%Parameters for arterioleODE
Cpass_specified = 220;
Cprimepass_specified = 11.467;
Cactive_specified = 274.193;
Cprimeactive_specified = 0.75;
Cdpactive_specified = 0.384;
Cmyo_specified = .159;%0.159;
C2_specified = 1.8;%.421;
CTGF_specified = 3e5;%4e5;
Ctone_specified = 12.59182581486616;
D0_specified = 2*Cpass_specified/(100*1333);

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% MAIN FILE for modeling autoregulation in the kidney
%
% WhAM! A Research Collaboration Workshop for Women in Applied Mathematics, 
%   Dynamical Systems with Applications to Biology and Medicine, IMA, Sept. 2013
% Team: Anita T. Layton, Julia Arciero, Laura Ellwein, Ashlee N. Ford Versypt, 
%   and Elizabeth Makrides
%
% This code has two parts: 
%   (1) arterioleODE: ODE model for afferent arteriole diameter D and smooth muscle
%       activation A as a result of the myogenic response and the
%       tubuloglomerular feedback mechanism. The SNGFR is a key output of
%       this part of the model and is the inlet flow rate for the second
%       part of the model.
%   (2) LoopHenlePDE: PDE model for the mass transport of chloride ions through the
%       TAL of the Loop of Henle. The cloride concentration profile at the
%       macular densa is the key output of this part of the model and is
%       fedback to the first part of the model with a time delay.
%
% Modified from code_from_WhAM_Sept_2013/LoopHenlePDE_MyoTGFmechFB_td.m
% Updated Jan 2014 at Duke during group reunion
%   Key changes from Dec 2013 kidney_code_v2.m
%   (1) Uses v2 ode45 solver over each time increment for the dynamic response
%       of the PDE. This numerical method does NOT use Strang splitting.
%   (2) The steady state initial values for the chloride concentration
%       distribution are defined via f_ss rather than LoopHenleconstant.m.
%       This also affects the initial & boundary conditions and the 
%       steady-state CMD values should be equal to Cop.
%   (3) Changed beta from 0.12 to 30/355.3977998342052
%   (4) Specify initial Q rather than Qc
%   (5) Four time intervals are used for (i) open-loop steady state,
%       (ii) closed-loop steady state with feedback, (iii) pulse with 
%       constant Q, and (iv) closed-loop dynamics
%   (6) Eliminated jdelay. CMD starts at index 1 of the NThistory vector,
%       which corresponds to index i=1 at tstart in the i=1:NT vector.
%   (7) Reintialized the .mat files for the control state with the revised 
%       chloride concentration as a reference/operating/control value.
%       This Cop was the steady state solution for f_ss at the operating
%       flow rate of SNGFR Q = 30nL/min, alpha = 0.2, F = 6nL/min.
%   Key changes from kidney_code_v4.m
%   (1) Pressure natriuresis including for updating alpha(Pa)
%   Key changes from kidney_code_v5.m
%   (1) Removed time periods 1-3 and pulse. Explored dynamics of changes in
%       pressure Pa.
% NOTE keep Pa vector of the form [constrolstate repeatconstrolstate otherpressures]
%   so that the plots for the control state start with a steady-state
%   chloride distribution
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Part 1: arteriole ODE
% Modified from code_from_WhAM_Sept_2013/FiveCompMRSRCR2.m written by
%   Julia & Laura
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Control State from FiveCompMRSRCR2.m in Sigmoid_code with Jan. 2014 CLc
% updated in FiveCompMRSRCR2.m and GetCStateMRSRCR2.m
%-------------------------------------------------------------------------
% function  [SSDTtPAQla, SSDTtPAQsa, Pa, Condsa, Condla, Actsa, Actla, Tensionsa, Tensionla,myosa,myola,tautermsa,tautermla,...
%         Shearsa, Shearla, Stonesa, Stonela, Diamsa, Diamla, Flowla, Flowsa, CRsignalSA, CRsignalLA, StonevalueSA, StonevalueLA,...
%         Cglobalsa, Prmsa, Prmla, Presssa, Pressla,consja,flow_tissue,extract,satstarta, satstartla,satstartsa,satstartc,satstartsv] ...
%         = FiveCompMRSRCR2(Myo_prox,Myo_dist,TGF_prox,TGF_dist,multiplier,capdensity)

% flag = 1;
%command line: FiveCompMRSRCR2(1,1,1,1,1,100)
FiveCompMRSRCR2options = [1,1,1,1,1,100];
%Fop = 30nL/min*.2 = 6nL/min;
%[xi,Si] = ode45(@f_ss,x,C0,[],Vmax,Km,p,Fop,r,C0,A1,A2,A3,L);
%Si(N) was hardcoded into FiveCompConsMRSRCR2 and GetCstate for determining
%the control state
% Myo_prox = FiveCompMRSRCR2options(1); Myo_dist = FiveCompMRSRCR2options(2);
% TGF_prox = FiveCompMRSRCR2options(3); TGF_dist = FiveCompMRSRCR2options(4);
multiplier = FiveCompMRSRCR2options(5); capdensity = FiveCompMRSRCR2options(6);
% %  Read in all of the input variables, creating MAT file FiveCompConsMRSRCR.mat and loading it into the workspace
% %  Also reads in the the levels of oxygen demand (M0)
% ReadInConsMRSRCR2(flag,Myo_prox,Myo_dist,TGF_prox,TGF_dist)
% starts with steady-state concentration gradient
load FiveCompConsMRSRCR2.mat

% %  Calculation of the complete control state, creating MAT file CStateValsMRSRCR.mat and loading it into the workspace
% GetCStateMRSRCR2(flag,multiplier,capdensity)
load CStateValsMRSRCR2.mat
% Qdistc
% Ddistc

% %  Read in pressure distribution
% ReadInPaMRSRCR2(Pac) with 1 P value: Psteps = 1; SDMS = 1;
% load Pavals2.mat

% %  Places control state in appropriate location of SS matrix
% InitPaWMRSRCR2(SMDS, Pa, Pac, NDS1, NDS2)
load InitPaWValsMRSRCR2.mat

% %  Name Start/End values for start of loop (so that correct values are referenced in loop if Pa goes up and down through entries)
% Start1 = Start;
% End1 = End;
% Again1 = Again;
% Inc1 = Inc;
% 
% %  initialize (or re-initialize) start/end/inc values
% Start = Start1;
% End = End1;
% Again = Again1; 
% Inc = Inc1;

%Loop over pressure vector
% while Again >= 0             %  Going through once if Again is 0.5, twice if 1.5  
%     
%     for pindex=Start:Inc:End
pindex = Start;
% index of pressure vector

% Call the function where the ODE solvers are used to find D and A
%-------------------------------------------------------------------------
% pre_MyoTGFmechFB(SSDTtPAQdist, pindex, Start, multiplier,capdensity);
%   Initializing iteration values and constants

if pindex > Start
    UoD = 1;
else
    UoD = -1;
end

%   Starting with the the appropriate initial values to start the iteration
%   depending on whether this is the first step from the control state or
%   a subsequent metabolic conductance step

if pindex == Start
    
    %     Dprox = Dproxc;
    Ddist = Ddistc;
    %     Tprox = Tproxc;
    Tdist = Tdistc;
    %     tauprox = tauproxc;
    %     taudist = taudistc;
    %     Pprox = 0.5 * (Pac + Pmid);
    Pdist = 0.5*(Pac + Pend);
    %     Aprox = Aproxc;
    Adist = Adistc;
    %     Qprox = Qproxc;
    Qdist = Qdistc;
    
else
    
    %     Dprox = SSDTtPAQla(i-UoD,1);
    Ddist = SSDTtPAQdist(pindex-UoD,1);
    %     Tprox = SSDTtPAQla(pindex-UoD,2);
    Tdist = SSDTtPAQdist(pindex-UoD,2);
    %     tauprox = SSDTtPAQla(pindex-UoD,3);
    %     taudist = SSDTtPAQsa(pindex-UoD,3);
    %     Pprox = SSDTtPAQla(pindex-UoD,4);
    Pdist = SSDTtPAQdist(pindex-UoD,3);
    %     Aprox = SSDTtPAQla(pindex-UoD,5);
    Adist = SSDTtPAQdist(pindex-UoD,4);
    %     Qprox = SSDTtPAQla(pindex-UoD,6);
    Qdist = SSDTtPAQdist(pindex-UoD,5);
    
end


%initial diameter conditions
% Dproxold = Dprox;
Ddistold = Ddist;

Qnum =  1/Cdistc;
Pdrop = Pa(pindex) - Pend;

u = [Ddistold Adist];
avgupast = mean(u,1);
toldsa = 10^-7;
tolasa = 10^-4;
toldla = toldsa;
tolala = tolasa;

%dynamic method of finding activation and diameter; loop through small time interval first (for non-oscillatory solutions)
% diameter_distal = [];
% activation_distal = [];
% % diameterla = [];
% % activationla = [];
% totaltime = [];

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Part 2: Loop of Henle PDE
% Modified from LoopHenle_PDE_simpleFBdelay.m written by Ashlee
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Dependent variable: tubular fluid Cl- concentration in thick ascending limb (TAL) 
%   of the loop of Henle in a single nephron
% C(x,t)        % [mmol/cm^3]     

%Parameters
%-------------------------------------------------------------------------
% Ref. Layton, Pitman, and Moore, 1991 eq. 12-13, Table 1 with backleak
L = 0.5;        % [cm] length of TAL
%Spatial discretization parameters
%-------------------------------------------------------------------------
dx = L/(N-1);   % [-]; spatial discretization step size
r = 10e-4;      % [cm] luminal radius of TAL
A = pi*r^2;     % [cm^2] luminal cross-sectional area 
% alpha_nominal = 0.20;   % [-] fraction of single-nephron golumerular filtration rate reaching TAL
Qop = 30;       % [nL/min] closed-loop steady-state single-nephron golumerular filtration rate
% Cop = 31.92e-3;    % [mmol/cm^3] steady-state Cl- concentration at the MD
% DQ = 18;        % [nL/min] TGF-mediated range of single-nephron golumerular filtration rate
% k = 0.24e3;     % [cm^3/mmol] scaling coefficient for TGF response
% Qc = 5.9233e-6; % [cm^3/s] afferent artiole flow rate before glomerular filtration
beta = 30/355.3977998342052;  % [-] fraction of artiole flow filtered in single nephron glomerulus
% NOTE: beta denominator for control state with CLc = 32.32606517560781 mM
% Q = Qc*beta;    % [nL/min] steady-state single-nephron gloumerular filtration rate
% Q=Qop;            % [nL/min] open-loop steady-state single-nephron gloumerular filtration rate
% Qc=Q/beta/60/10^6;% [cm^3/s] afferent artiole flow rate before glomerular filtration
%F(t)           % [cm^3/s] tubular fluid flow rate, can be a function of time
% F = alpha_nominal*beta*Qc;  % [cm^3/s]  tubular fluid flow rate

%Independent variables
%-------------------------------------------------------------------------
%t              % [s] time
x = 0:dx:L;     % [cm] axial position along TAL

%Specialized functions
%-------------------------------------------------------------------------
%Ce(x)          % [mmol/cm^3] interstitial Cl- concentration
% Ref. Layton, Pitman, and Moore, 1991 eq. 12-13
A3 = 2;
CeN = 150e-3;  
A1 = ( 1 - CeN / C0 ) ./ ( 1 - exp( -A3 ) );
A2 = 1 - A1;
Ce(1:N) = C0 .* ( A1 .* exp( -A3 .* x(1:N) ./ L ) + A2 ) ;

%Initial condition
Activation(1) = Adistc;
diameter_finali(1) = Ddistc;
% [xi,Si] = ode45(@f_ss,x,C0,[],Vmax,Km,p,F,r,C0,A1,A2,A3,L);

% C_0=Si;
% % C_0 = LoopHenlePDEconstant(C0,Q,N);% C0ss(2:N);  % [mmol/cm^3]  
% %Boundary condition
% C_0(1) = C0;
%Boundary condition
C_0(1) = C0;

%Initial condition (update the 
%No initial concentration gradient
C_0(2:N) = C0; % [mmol/cm^3, mM] 


%Compute stochastic effects of tau on index for read appropriating CMD at
%stochastically varying lag time
% Note: this is deterministic for now
%-------------------------------------------------------------------------
%jdelay = ones(NT,1);
% jvector = 1:NThistory;
% jdelay(1:NT) = jvector(1+int_tau0:NT+int_tau0) - int_tau0;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Combine parts 1 and 2 of model with feedback between the parts
% Simultaneously solve arterioleODE and LoopHenlePDE systems of ODEs
%
% Modified from FiveCompMRSRCR2.m and LoopHenle_PDE_simpleFBdelay.m
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

tic
%Concatenate vectors dCdt, dDdt, dAdt
%-------------------------------------------------------------------------
% load pre_MyoTGF.mat

%check if specified parameters for arterioleODE match those for the control
%state from ReadInConsMRSRCR2(flag,Myo_prox,Myo_dist,TGF_prox,TGF_dist) 
%& FiveCompConsMRSRCR2.mat
if Params_dist(1) ~= Cpass_specified
    Params_dist(1) = Cpass_specified
end
if Params_dist(2) ~= Cprimepass_specified
    Params_dist(2) = Cprimepass_specified
end
if Params_dist(3) ~= Cactive_specified
    Params_dist(3) = Cactive_specified
end
if Params_dist(4) ~= Cprimeactive_specified
    Params_dist(4) = Cprimeactive_specified
end
if Params_dist(5) ~= Cdpactive_specified
    Params_dist(5) = Cdpactive_specified
end
if Params_dist(6) ~= Cmyo_specified
    Params_dist(6) = Cmyo_specified

end
if Params_dist(7) ~= C2_specified
    Params_dist(7) = C2_specified
end
if Params_dist(8) ~= CTGF_specified
    Params_dist(8) = CTGF_specified
end
if Params_dist(9) ~= Ctone_specified
    Params_dist(9) = Ctone_specified
end
if Params_dist(10) ~= D0_specified
    Params_dist(10) = D0_specified
end
if MyoOn_dist~=MyoOn
    MyoOn_dist=MyoOn
end
if TGFOn_dist~=TGFOn
    TGFOn_dist=TGFOn
end


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Initialization at control state
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
pindex = 1;
P_alpha=Pa(pindex);
P0 = Pac;
if constantalpha == 0
    alpha=1-0.8*(1-0.45*(max(P_alpha,80)-P0)/P0);  % [-] fraction of single-nephron golumerular filtration rate reaching TAL
end
F = alpha*beta*Qdist;% [cm^3/s]  tubular fluid flow rate
Qc=Qdist;% [cm^3/s] afferent artiole flow rate before glomerular filtration
Q=Qc*beta*60*10^6;% [nL/min] open-loop steady-state single-nephron gloumerular filtration rate
%-------------------------------------------------------------------------
%Initial steady state concentration gradient in Loop of Henle at control
%pressure
%-------------------------------------------------------------------------
tstart = 0;     % [s] starting simulation time
t0 = tstart;    % [s] starting time for first call to ODE solver
tf = t0 + tinc; % [s] time after one ODE solver call increment
time = tstart:tinc:tfinal;
NT = length(time);% [-] number of time points
%Initial values coming into ODE solver call loop
Csoln = ones(NT,N);
Csoln(1,:) = C_0(:);
%preallocate vectors
QC= zeros(NT,length(Pa));
Fi = ones(NT,length(Pa));
Activationvector = ones(NT,length(Pa));
diametervector = ones(NT,length(Pa));
CMDplotvector = ones(NT,length(Pa));
tspan = [t0 tf];
%CMD history vector, chloride concentration at the macula densa
%-------------------------------------------------------------------------
int_tau0 = round(tau0/tinc);
if int_tau0 == 0
    int_tau0 = 1;
end
NThistory = int_tau0 + NT;
CMD = zeros(NThistory,1);
CMD(1:int_tau0+1) = C_0(N);
% for tt = 1:NT
%     Csoln(tt,:) = C_0(:);
% end
for i = 2:NT
    [tis,Cis] = ode45(@LoopHenlePDE,tspan,C_0,[],N,dx,A,Vmax,Km,p,Ce,r,F);
    %Concatenate the solution vectors  
    a=size(Cis,1);
    Csoln(i,:) = Cis(a,:);
    %update times & initial condition for next call to ODE solver
    t0= tf;
    tf= t0 + tinc; 
    tspan = [t0 tf];
    C_0 = Cis(a,:);
    CMD(int_tau0+i) = Cis(a,N);
end



CMDplotvector(1:NT,1)=Csoln(:,N);
%-------------------------------------------------------------------------
%Steady state control state at a constant chloride concentration
%-------------------------------------------------------------------------
tstart = 0;     % [s] starting simulation time
t0 = tstart;    % [s] starting time for first call to ODE solver
tf = t0 + tinc; % [s] time after one ODE solver call increment
tspan = [t0 tf];
uis = [Ddistc Adistc];
for i = 2:NT
    [tis,uis] = ode15s(@arterioleODE,tspan,uis(end,:),[],Pa(1),Qnum,Pdrop,...
        multiplier,capdensity,Cis(end)/1000,td,Pac,Pend,...
        MyoOn_dist,TGFOn_dist,ndistc,mu_dist,Leff_dist,Qdistc,...
        Pdistc,Params_dist,CLc,ta);
    %Concatenate the solution vectors  
    a=size(uis,1);
    usteady(i,:) = uis(a,:);
    %update times & initial condition for next call to ODE solver
    t0 = tf;
    tf = t0 + tinc; 
    tspan = [t0 tf];
end
    
    diameter_distalsteady = usteady(:,1);
    activation_distalsteady = usteady(:,2);
    
    Ddistnewisteady = diameter_distalsteady;

    Cdististeady = (pi.*Ddistnewisteady.^4*ndistc)./(128*mu_dist*Leff_dist);

%     Qdeni = 1/Cdisti;
%     Qdisti = Qdistc*((Pa(pindex)-Pend)/(Pac-Pend))*(Qnum/Qdeni);
    Qdististeady = (Pa(pindex)-Pend)*Cdististeady;

%     deltaPdisti = Pa(pindex) - Pend;
    Pdististeady = 0.5*(Pa(pindex)+Pend);

    Tdististeady = Pdististeady*Ddistnewisteady/2;

    % % %control vals
if MyoOn_dist == 1
    Tsigdististeady  = Tdististeady ;
else
    %   Tsigla = interp1(DiameterLA,TensionLA,Dlanew,'linear');
    Tsigdististeady  = Tdistc; %control state
%     Tsigdist = Tdist;  %zeroed out
end

    Adististeady = activation_distalsteady;

    Activationsteady = Adististeady;
    perfusionisteady = Qdististeady;
    diameter_finalisteady = Ddistnewisteady;

    QCsteady = perfusionisteady;
    Fisteady = alpha*beta*QCsteady;
    QC(:,1) = QCsteady;
    Fi(:,1) = Fisteady;
    Activationvector(:,1) = Activationsteady;
    diametervector(:,1) = diameter_finalisteady;
%initial conditions for pressure loop
    u(end,:) = usteady(end,:);
    Y0 = C_0;
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Pressure loop
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
for pindex=2:length(Pa)
    %-------------------------------------------------------------------------
    %Time parameters
    %-------------------------------------------------------------------------
    tstart = 0;     % [s] starting simulation time
    t0 = tstart;    % [s] starting time for first call to ODE solver
    tf = t0 + tinc; % [s] time after one ODE solver call increment
    %NT = 1001;      % [-] number of time points
    tspan = [t0 tf];

    P_alpha=Pa(pindex);
    P0 = Pac;
    if constantalpha == 0
        alpha=1-0.8*(1-0.45*(max(P_alpha,P0)-P0)/P0);  % [-] fraction of single-nephron golumerular filtration rate reaching TAL
    end

    %ODE solver loop to accomodate time delays in signal CMD. The overall 
    %Delta t value is specified as tinc for the time between each index tt
    %the ODE solvers can used finer adaptive time points within the desired
    %tspan. Update A & D ODEs each Delta t.

%     CMD = zeros(NThistory,1);
%     CMD(1:int_tau0+1) = C_0(N);
    CMD(1:int_tau0+1) = CMD(NT:NThistory);

    %   Initial condition for pindex loop is specified as the values from the previous final state     
    tt = 1;

        diameter_distal(tt) = u(end,1);
        activation_distal(tt) = u(end,2);

        Ddistnewi = diameter_distal(tt);

        Cdisti = (pi*Ddistnewi^4*ndistc)/(128*mu_dist*Leff_dist);

    %     Qdeni = 1/Cdisti;
    %     Qdisti = Qdistc*((Pa(pindex)-Pend)/(Pac-Pend))*(Qnum/Qdeni);
        Qdisti = (Pa(pindex)-Pend)*Cdisti;

    %     deltaPdisti = Pa(pindex) - Pend;
        Pdisti = 0.5*(Pa(pindex)+Pend);

        Tdisti = Pdisti*Ddistnewi/2;

        % % %control vals
        if MyoOn_dist == 1
            Tsigdisti = Tdisti;
        else
            %   Tsigla = interp1(DiameterLA,TensionLA,Dlanew,'linear');
            Tsigdisti = Tdistc; %control state
        %     Tsigdist = Tdist;  %zeroed out
        end

        Adisti = activation_distal(tt);

        SSDTtPAQdisti(pindex,1) = Ddistnewi;
        SSDTtPAQdisti(pindex,2) = Tsigdisti;
        SSDTtPAQdisti(pindex,3) = Pdisti;
        SSDTtPAQdisti(pindex,4) = Adisti;
        SSDTtPAQdisti(pindex,5) = Qdisti;

        Activation(tt) = Adisti;
        perfusioni(tt) = SSDTtPAQdisti(pindex,5);
        diameter_finali(tt) = SSDTtPAQdisti(pindex,1);
        QC(tt,pindex) = perfusioni(tt);
        Fi(tt,pindex) = alpha*beta*QC(tt,pindex);
        Activationvector(tt,pindex) = Activation(tt);
        diametervector(tt,pindex) = diameter_finali(tt);;
        CMD(int_tau0+tt) = Y0(N);
        CMDplotvector(tt,pindex)=Y0(N);
    %Concatenate the solution vectors    
        Csoln(tt,:) = Y0;
        %update times & initial condition for next call to ODE solver
        t0 = tf;
        tf = t0 + tinc; 
        tspan = [t0 tf];
        
    for tt = 2:NT %time index
        %---------------------------------------------------------------------
        % arterioleODE for Delta t
        %---------------------------------------------------------------------
        [t,u] = ode15s(@arterioleODE,tspan,u(end,:),[],Pa(pindex),Qnum,Pdrop,...
            multiplier,capdensity,CMD(tt)/1000,td,Pac,Pend,...
            MyoOn_dist,TGFOn_dist,ndistc,mu_dist,Leff_dist,Qdistc,...
            Pdistc,Params_dist,CLc,ta);

        diameter_distal(tt) = u(end,1);
        activation_distal(tt) = u(end,2);

        Ddistnewi = diameter_distal(tt);

        Cdisti = (pi*Ddistnewi^4*ndistc)/(128*mu_dist*Leff_dist);

    %     Qdeni = 1/Cdisti;
    %     Qdisti = Qdistc*((Pa(pindex)-Pend)/(Pac-Pend))*(Qnum/Qdeni);
        Qdisti = (Pa(pindex)-Pend)*Cdisti;

    %     deltaPdisti = Pa(pindex) - Pend;
        Pdisti = 0.5*(Pa(pindex)+Pend);

        Tdisti = Pdisti*Ddistnewi/2;

        % % %control vals
        if MyoOn_dist == 1
            Tsigdisti = Tdisti;
        else
            %   Tsigla = interp1(DiameterLA,TensionLA,Dlanew,'linear');
            Tsigdisti = Tdistc; %control state
        %     Tsigdist = Tdist;  %zeroed out
        end

        Adisti = activation_distal(tt);

        SSDTtPAQdisti(pindex,1) = Ddistnewi;
        SSDTtPAQdisti(pindex,2) = Tsigdisti;
        SSDTtPAQdisti(pindex,3) = Pdisti;
        SSDTtPAQdisti(pindex,4) = Adisti;
        SSDTtPAQdisti(pindex,5) = Qdisti;

        Activation(tt) = Adisti;
        perfusioni(tt) = SSDTtPAQdisti(pindex,5);
        diameter_finali(tt) = SSDTtPAQdisti(pindex,1);
        QC(tt,pindex) = perfusioni(tt);
        Fi(tt,pindex) = alpha*beta*QC(tt,pindex);
        Activationvector(tt,pindex) = Activation(tt);
        diametervector(tt,pindex) = diameter_finali(tt);
        %---------------------------------------------------------------------
        % LoopHenlePDE for Delta t
        %---------------------------------------------------------------------
        [ti,Yi] = ode45(@LoopHenlePDE,tspan,Y0,[],N,dx,A,Vmax,Km,p,Ce,r,Fi(tt,pindex));

        %---------------------------------------------------------------------
        % finalize this time step
        %---------------------------------------------------------------------
        a=size(Yi,1);

    %Cl- concentration feedback to afferent artery through macula densa (MD) and
    %juxtaglomerular apparatus (JGA)    
        CMD(int_tau0+tt) = Yi(a,N);
        CMDplotvector(tt,pindex)=Yi(a,N);
    %Concatenate the solution vectors    
        Csoln(tt,:) = Yi(a,:);
        %update times & initial condition for next call to ODE solver
        t0 = tf;
        tf = t0 + tinc; 
        tspan = [t0 tf];
        Y0 = Yi(a,:);
    %     u0 = u(a,:);
    end
end
% figure(1)
% plot(x/L,Csoln(1,:).*1000,'g');
% hold on
% plot(x/L,Csoln(NT,:).*1000,'r');
% hold off
% legend('t = 0s',['t = ' num2str(tfinal) 's'],'Location','SouthWest')
% xlabel('x/L')
% ylabel('[ Cl^- ] (mM)')
% 
% figure(2)
% plot(time,Csoln(:,1).*1000,'g')
% hold on
% plot(time,Csoln(:,(N-1)/2+1).*1000,'b')
% plot(time,Csoln(:,N).*1000,'r')
% hold off
% xlabel('t')
% ylabel('[ Cl^- ] (mM)')
% legend('x/L = 0',['x/L = ' num2str(x((N-1)/2+1)/L)],['x/L = ' num2str(x(N)/L)])
% 
% figure(3)
% plot(x/L,Ce*1000)
% xlabel('x/L')
% ylabel('C_e (mM)')
% % check: Ce(N) == CeN
% figure(20)
% QC(1)=Qdistc;
% plot(time,QC*60*10^6)
% xlabel('t (s)')
% ylabel('Q_{AA} (nL/min)')
% 


previnterval = 100; %[s] maximum width between sustained peaks 
QCtol = 1e-3;
QCavg=zeros(length(Pa),1);
Aavg=zeros(length(Pa),1);
Davg=zeros(length(Pa),1);
QCpeakdiff=zeros(length(Pa),1);
timepeakdiff=zeros(length(Pa),1);
for pindex = 1:length(Pa)
max_xvals=findmaxima(QC(:,pindex)./QC(end,1));
%remove the last time point from the list of max vals in case it's not
%really a peak but simply a final time with a positive slope.
if max_xvals(end) == NT
    max_xvals = max_xvals(1:end-1);
end
% time(max_xvals(end))
% time(max_xvals(end-1))
    if length(max_xvals)>1
        QCavg(pindex) = average(QC(max_xvals(end-1)+1:max_xvals(end),pindex)./QC(end,1));
        Davg(pindex) = average(diametervector(max_xvals(end-1)+1:max_xvals(end),pindex));
        Aavg(pindex) = average(Activationvector(max_xvals(end-1)+1:max_xvals(end),pindex));
        QCpeakdiff(pindex) = abs((QC(max_xvals(end-1),pindex)-QC(max_xvals(end),pindex))./QC(end,1));
        if QCpeakdiff(pindex) > QCtol
            ['Not at steady-state for pindex = ' num2str(pindex) '. Extend tfinal and replot.']
            beep
        end
        timepeakdiff(pindex) = abs(time(max_xvals(end))-time(max_xvals(end-1)));
        if timepeakdiff(pindex) >previnterval
            ['The period between the final two maxima is longer than ' num2str(previnterval) 's for pindex = ' num2str(pindex) '.']
            QCavgdiff = abs(QC(end,pindex)./QC(end,1)-QC(end-1,pindex)./QC(end,1));
            if QCavgdiff < QCtol
                QCavg(pindex) = QC(end,pindex)./QC(end,1);
            end
            beep
        end
    else
        ['Solution only has one local maxima at the initial condition for pindex = ' num2str(pindex) '.']
    end
end

figure(31)
plot(time,CMDplotvector(:,2)*1000,'k')


temp = horzcat(time', CMDplotvector(:,2:end)*1000);

csvwrite('figure31.dat', temp, 1, 0)
hold on
plot(time,CMDplotvector(:,3)*1000,'b')
plot(time,CMDplotvector(:,4)*1000,'r')
plot(time,CMDplotvector(:,5)*1000,'g')
plot(time,CMDplotvector(:,6)*1000,'k--')
plot(time,CMDplotvector(:,7)*1000,'b--')
%% 
plot(time,CMDplotvector(:,8)*1000,'r--')
plot(time,CMDplotvector(:,9)*1000,'g--')
%% 
plot(time,CMDplotvector(:,10)*1000,'k:')
%% 
hold off
xlabel('t (s)')
ylabel('CMD (mM)')
legend(['P = ' num2str(Pa(2)/1333) 'mmHg'], ['P = ' num2str(Pa(3)/1333) 'mmHg'],['P = ' num2str(Pa(4)/1333) 'mmHg'],...
    ['P = ' num2str(Pa(5)/1333) 'mmHg'], ['P = ' num2str(Pa(6)/1333) 'mmHg'],['P = ' num2str(Pa(7)/1333) 'mmHg'],...
    ['P = ' num2str(Pa(8)/1333) 'mmHg'], ['P = ' num2str(Pa(9)/1333) 'mmHg'],['P = ' num2str(Pa(10)/1333) 'mmHg'], 'Location','Best')



figure(33)
newQC = QC * 60 *10^6*beta;
temp = horzcat(time', newQC(:,2:end));

csvwrite('figure33.dat', temp, 1, 0)
plot(time,QC(:,2)*60*10^6*beta,'k')
hold on
plot(time,QC(:,3)*60*10^6*beta,'b')
plot(time,QC(:,4)*60*10^6*beta,'r')
plot(time,QC(:,5)*60*10^6*beta,'g')
plot(time,QC(:,6)*60*10^6*beta,'k--')
plot(time,QC(:,7)*60*10^6*beta,'b--')
plot(time,QC(:,8)*60*10^6*beta,'r--')
plot(time,QC(:,9)*60*10^6*beta,'g--')
plot(time,QC(:,10)*60*10^6*beta,'k:')
hold off
xlabel('t (s)')
ylabel('Q (nL/min)')
legend(['P = ' num2str(Pa(2)/1333) 'mmHg'], ['P = ' num2str(Pa(3)/1333) 'mmHg'],['P = ' num2str(Pa(4)/1333) 'mmHg'],...
    ['P = ' num2str(Pa(5)/1333) 'mmHg'], ['P = ' num2str(Pa(6)/1333) 'mmHg'],['P = ' num2str(Pa(7)/1333) 'mmHg'],...
    ['P = ' num2str(Pa(8)/1333) 'mmHg'], ['P = ' num2str(Pa(9)/1333) 'mmHg'],['P = ' num2str(Pa(10)/1333) 'mmHg'], 'Location','Best')




figure(21)

temp = horzcat(time', Activationvector(:,2:end));

csvwrite('figure21.dat', temp, 1, 0)

plot(time,Activationvector(:,2),'k')
hold on
plot(time,Activationvector(:,3),'b')
plot(time,Activationvector(:,4),'r')
plot(time,Activationvector(:,5),'g')
plot(time,Activationvector(:,6),'k--')
plot(time,Activationvector(:,7),'b--')
plot(time,Activationvector(:,8),'r--')
plot(time,Activationvector(:,9),'g--')
plot(time,Activationvector(:,10),'k:')
hold off
xlabel('t (s)')
ylabel('Activation')
legend(['P = ' num2str(Pa(2)/1333) 'mmHg'], ['P = ' num2str(Pa(3)/1333) 'mmHg'],['P = ' num2str(Pa(4)/1333) 'mmHg'],...
    ['P = ' num2str(Pa(5)/1333) 'mmHg'], ['P = ' num2str(Pa(6)/1333) 'mmHg'],['P = ' num2str(Pa(7)/1333) 'mmHg'],...
    ['P = ' num2str(Pa(8)/1333) 'mmHg'], ['P = ' num2str(Pa(9)/1333) 'mmHg'],['P = ' num2str(Pa(10)/1333) 'mmHg'], 'Location','Best')


figure(22)
newdiametervector = diametervector * 1e4;
temp = horzcat(time', newdiametervector(:,2:end));

csvwrite('figure22.dat', temp, 1, 0)
plot(time,diametervector(:,2)*1e4,'k')
hold on
plot(time,diametervector(:,3)*1e4,'b')
plot(time,diametervector(:,4)*1e4,'r')
plot(time,diametervector(:,5)*1e4,'g')
plot(time,diametervector(:,6)*1e4,'k--')
plot(time,diametervector(:,7)*1e4,'b--')
plot(time,diametervector(:,8)*1e4,'r--')
plot(time,diametervector(:,9)*1e4,'g--')
plot(time,diametervector(:,10)*1e4,'k:')
hold off
xlabel('t (s)')
ylabel('diameter (\mu m)')
legend(['P = ' num2str(Pa(2)/1333) 'mmHg'], ['P = ' num2str(Pa(3)/1333) 'mmHg'],['P = ' num2str(Pa(4)/1333) 'mmHg'],...
    ['P = ' num2str(Pa(5)/1333) 'mmHg'], ['P = ' num2str(Pa(6)/1333) 'mmHg'],['P = ' num2str(Pa(7)/1333) 'mmHg'],...
    ['P = ' num2str(Pa(8)/1333) 'mmHg'], ['P = ' num2str(Pa(9)/1333) 'mmHg'],['P = ' num2str(Pa(10)/1333) 'mmHg'], 'Location','Best')



figure(23)
legendStrings = {['Pressure (mmHg)'], ['Steady State Normalized SNGFR']}



temp = horzcat(Pa(2:end)'./1333, QCavg(2: end));

csvwrite('figure23.dat', temp, 1, 0)
plot(Pa(2:end)./1333,QCavg(2:end))
xlabel('Pressure (mmHg)')
ylabel('Steady State Normalized SNGFR')



figure(24)

temp = horzcat(Pa(2:end)'./1333, Aavg(2: end));
csvwrite('figure24.dat', temp, 1, 0)
plot(Pa(2:end)./1333,Aavg(2:end))
xlabel('Pressure (mmHg)')
ylabel('Steady State Activation')

legendStrings = {['Pressure (mmHg)'], ['Steady State Activation']}



figure(25)
temp = horzcat(Pa(2:end)'./1333, Davg(2: end).*1e4);
csvwrite('figure25.dat', temp, 1, 0)
plot(Pa(2:end)./1333,Davg(2:end).*1e4)
xlabel('Pressure (mmHg)')
ylabel('Steady State Diameter')
legendStrings = {['Pressure (mmHg)'], ['Steady State Diameter']}


% figure(32)
% plot(time,Csoln(:,1).*1000,'g')
% hold on
% plot(time,Csoln(:,(N-1)/2+1).*1000,'b')
% plot(time,Csoln(:,N).*1000,'r')
% hold off
% xlabel('t (s)')
% ylabel('[ Cl^- ] (mM)')
% title('Steady state chloride concentration at control flow rate')
% legend('x/L = 0',['x/L = ' num2str(x((N-1)/2+1)/L)],['x/L = ' num2str(x(N)/L)])

save(outputfilename)
% 
% %         load DAsolver2.mat
% %         
%         
% %         post_MyoTGFmechFB(SSDTtPAQdist, pindex, Start, multiplier,capdensity);
% % % %control vals
% if MyoOn_dist == 1
%     Tsigdisti = Tdisti;
% else
%     %   Tsigla = interp1(DiameterLA,TensionLA,Dlanew,'linear');
%     Tsigdisti = Tdistc; %control state
% %     Tsigdist = Tdist;  %zeroed out
% end
% %         load DAvals2.mat 
% %         load DAallvals2.mat
%         
% %         flow_tissue(pindex) = perfusion;
% %         Ds(pindex) = diameter_final;   
toc

%---------------------------------
function dCdt = LoopHenlePDE(t,C,N,dx,A,Vmax,Km,p,Ce,r,F)
%DCDT

%Boundary condition 
dCdt(1) = 0;    %constant
%Method of lines
%Discretize the spatial derivative with upwind finite difference method
%Solve characteristic ODEs with Matlab ODE solver
dCdt(2:N) = ( -F .* ( C(2:N) - C(1:(N-1) ) ) ./ dx - 2*pi*r .* ( Vmax .* C(2:N) ./ ...
    ( Km + C(2:N) ) + p .* ( C(2:N)  - Ce(2:N)' ) ) ) ./ A;

dCdt = dCdt';

function avg = average(x)
avg = sum(x)/length(x);