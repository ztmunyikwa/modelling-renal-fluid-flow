function kidney_code_v4
%USE THIS VERSION FOR TRANSIENT PULSES STARTING FROM CONTROL STATE.
%DO NOT CHANGE PRESSURE

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
tfinal = 800;        % [s]   final simulation time; 60s simulation takes 
%                           45s to solve
C0 = 275e-3;        % [mmol/cm^3] C(0,t); boundary concentration of Cl-
%                           at loop bend
Vmax = 14.5e-6;     % [mmol/cm^2/s] 
%                           maximum transport rate of Cl- from TAL    
Km = 70e-3;         % [mmol/cm^3] 
%                           active solute transport Michaelis constant
p = 1.5e-5;         % [cm/s] trans-epithelial diffusion backleak 
%                           chloride permeability
pulseduration = 15.7;%[s]   duration of the pulse interval; set to 0 for 
%                           immediate step change in values
Qpulsepercent = 1;  % [-]   percent that Qpulse will be changed from Q value
tau0 = 4;           % [s]   time delay at juxtaglomerular apparatus (JGA)
%                           code sets the minimum tau0 to one time step
td = 1;             % [s]   time constant for diameter ODE
ta = 100;            % [s]   time constant for activation ODE
%Note: ta should always be greater than td
CLc = 32.32606517560781;

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
Pa = [133300 133300];        % [dynes/cm^2] incoming pressure

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
%-------------------------------------------------------------------------
%Time parameters
%-------------------------------------------------------------------------
tstart = 0;     % [s] starting simulation time
t0 = tstart;    % [s] starting time for first call to ODE solver
tf = t0 + tinc; % [s] time after one ODE solver call increment
time = tstart:tinc:tfinal;
NT = length(time);
%NT = 1001;      % [-] number of time points
tspan = [t0 tf];
% ATOL=1e-6;
% RTOL=1e-3;

%Parameters
%-------------------------------------------------------------------------
% Ref. Layton, Pitman, and Moore, 1991 eq. 12-13, Table 1 with backleak
L = 0.5;        % [cm] length of TAL
%Spatial discretization parameters
%-------------------------------------------------------------------------
dx = L/(N-1);   % [-]; spatial discretization step size
r = 10e-4;      % [cm] luminal radius of TAL
A = pi*r^2;     % [cm^2] luminal cross-sectional area 
alpha = 0.20;   % [-] fraction of single-nephron golumerular filtration rate reaching TAL
Qop = 30;       % [nL/min] closed-loop steady-state single-nephron golumerular filtration rate
% Cop = 31.92e-3;    % [mmol/cm^3] steady-state Cl- concentration at the MD
% DQ = 18;        % [nL/min] TGF-mediated range of single-nephron golumerular filtration rate
% k = 0.24e3;     % [cm^3/mmol] scaling coefficient for TGF response
% Qc = 5.9233e-6; % [cm^3/s] afferent artiole flow rate before glomerular filtration
beta = 30/355.3977998342052;  % [-] fraction of artiole flow filtered in single nephron glomerulus
% NOTE: beta denominator for control state with CLc = 32.32606517560781 mM
% Q = Qc*beta;    % [nL/min] steady-state single-nephron gloumerular filtration rate
Q=Qop;            % [nL/min] open-loop steady-state single-nephron gloumerular filtration rate
Qc=Q/beta/60/10^6;% [cm^3/s] afferent artiole flow rate before glomerular filtration
%F(t)           % [cm^3/s] tubular fluid flow rate, can be a function of time
F = alpha*beta*Qc;  % [cm^3/s]  tubular fluid flow rate

%Time periods
%first: open-loop steady state
openloopsstime = 5;
openloopssindex = 1+round(openloopsstime/tinc);
%second: closed-loop steady state
closedloopsstime = 5+openloopsstime;
closedloopssindex = 1+round(closedloopsstime/tinc);
%third: constant pulse, specified as Fpulse in constant open-loop
pulset0= pulseduration+closedloopsstime;%15.7
pulseindex = 1+round(pulset0/tinc)
% tpulse=time(pulseindex);
Qpulse = (1+Qpulsepercent/100)*Q
FnLpulse = alpha*Qpulse;  % [nL/min] tubular fluid flow rate
Fpulse = FnLpulse/60/10^6;
%fourth: closed-loop dynamics

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

%preallocate vectors
QC= zeros(NT,1);
Fi = ones(NT,1);

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

%ODE solver loop to accomodate time delays in signal CMD. The overall 
%Delta t value is specified as tinc for the time between each index tt
%the ODE solvers can used finer adaptive time points within the desired
%tspan. Update A & D ODEs each Delta t.
%-------------------------------------------------------------------------
%Initial steady state concentration gradient in Loop of Henle at control
%pressure
%-------------------------------------------------------------------------
tstartsteady = 0;     % [s] starting simulation time
t0steady = tstartsteady;    % [s] starting time for first call to ODE solver
tfsteady = t0steady + tinc; % [s] time after one ODE solver call increment
tfinalsteady=100;
timesteady = tstartsteady:tinc:tfinalsteady;
NTsteady = length(timesteady);
%NT = 1001;      % [-] number of time points
tspansteady = [t0steady tfsteady];
for i = 2:NTsteady
    [tis,Cis] = ode45(@LoopHenlePDE,tspansteady,C_0,[],N,dx,A,Vmax,Km,p,Ce,r,F);
    %Concatenate the solution vectors  
    a=size(Cis,1);
    Csteady(i,:) = Cis(a,:);
    %update times & initial condition for next call to ODE solver
    t0steady = tfsteady;
    tfsteady = t0steady + tinc; 
    tspansteady = [t0steady tfsteady];
    C_0 = Cis(a,:);
end

'Csteady (mM)'
Csteady(end,N)*1000
figure(31)
plot(timesteady,Csteady(:,N)*1000)
xlabel('t (s)')
ylabel('CMD (mM)')
title('Steady state chloride concentration at control flow rate')
figure(32)
plot(timesteady,Csteady(:,1).*1000,'g')
hold on
plot(timesteady,Csteady(:,(N-1)/2+1).*1000,'b')
plot(timesteady,Csteady(:,N).*1000,'r')
hold off
xlabel('t (s)')
ylabel('[ Cl^- ] (mM)')
title('Steady state chloride concentration at control flow rate')
legend('x/L = 0',['x/L = ' num2str(x((N-1)/2+1)/L)],['x/L = ' num2str(x(N)/L)])
%CMD history vector, chloride concentration at the macula densa
%-------------------------------------------------------------------------
int_tau0 = round(tau0/tinc);
if int_tau0 == 0
    int_tau0 = 1;
end
NThistory = int_tau0 + NT;
CMD = zeros(NThistory,1);
CMD(1:int_tau0+1) = C_0(N);
%Initial values coming into ODE solver call loop
Csoln = ones(NT,N);
Csoln(1,:) = C_0(:);
% for tt = 1:NT
%     Csoln(tt,:) = C_0(:);
% end
%-------------------------------------------------------------------------
%Steady state control state at a constant chloride concentration
%-------------------------------------------------------------------------
tstartsteady = 0;     % [s] starting simulation time
t0steady = tstartsteady;    % [s] starting time for first call to ODE solver
tfsteady = t0steady + tinc; % [s] time after one ODE solver call increment
tfinalsteady=100;
timesteady = tstartsteady:tinc:tfinalsteady;
NTsteady = length(timesteady);
%NT = 1001;      % [-] number of time points
tspansteady = [t0steady tfsteady];
uis = [Ddistc Adistc];
for i = 2:NTsteady
    [tis,uis] = ode15s(@arterioleODE,tspansteady,uis(end,:),[],Pa(pindex),Qnum,Pdrop,...
        multiplier,capdensity,Csteady(end,N)/1000,td,Pac,Pend,...
        MyoOn_dist,TGFOn_dist,ndistc,mu_dist,Leff_dist,Qdistc,...
        Pdistc,Params_dist,CLc,ta);
    %Concatenate the solution vectors  
    a=size(uis,1);
    usteady(i,:) = uis(a,:);
    %update times & initial condition for next call to ODE solver
    t0steady = tfsteady;
    tfsteady = t0steady + tinc; 
    tspansteady = [t0steady tfsteady];
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
    figure(33)
    plot(timesteady,Fisteady*60*10^6/alpha)
    xlabel('t (s)')
    ylabel('Q (nL/min)')
    title('SNGFR at control state with steady C(L,t)')
    
%---------------------------------------------------------------------
% First time period
% Open-loop steady state: Already computed above by taking the long time 
% solution of the steady-state chloride distribution using a uniform  
% initial condition C_0. The results are just copied here to visualize 
% the open loopp steady state.
% Does not call arteriole here because the feedback loop is not closed
%---------------------------------------------------------------------
for i = 2:openloopssindex
    C1i(1:N) = Csoln(i-1,:);
    %Cl- concentration feedback to afferent artery through macula densa (MD) and
    %juxtaglomerular apparatus (JGA)    
    CMD(int_tau0+i) = C1i(N);
    %Concatenate the solution vectors  
    Csoln(i,:) = C1i(1:N);
    %update times & initial condition for next call to ODE solver
    t0 = tf;
    tf = t0 + tinc; 
    tspan = [t0 tf];
    C_0 = Csoln(i,:);
    diameter_distal(i) = Ddistc;
    activation_distal(i) = Adistc;

%     diameter_distal = [diameter_distal; u(:,1)];
%     activation_distal = [activation_distal; u(:,2)];
%     
    Ddistnewi = diameter_distal(end);

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


    Adisti = activation_distal(end);

    SSDTtPAQdisti(pindex,1) = Ddistnewi;
    SSDTtPAQdisti(pindex,2) = Tsigdisti;
    SSDTtPAQdisti(pindex,3) = Pdisti;
    SSDTtPAQdisti(pindex,4) = Adisti;
    SSDTtPAQdisti(pindex,5) = Qdisti;

    Activation(i) = Adisti;
    perfusioni(i) = SSDTtPAQdisti(pindex,5);
    diameter_finali(i) = SSDTtPAQdisti(pindex,1);   
    QC(i) = perfusioni(i);
    Fi(i) = alpha*beta*QC(i);
end
%---------------------------------------------------------------------
% Second time period
% Closed-loop steady state: Starts with the steady-state control state
% solution computed above. Closes the feedback loop so that updated values
% of CMD are fed into the arterioleODE function dynamically. This period 
% should be constant if a steady-state existed in the initial conditions.
%---------------------------------------------------------------------
u(end,:) = uis(end,:);
for tt = (openloopssindex+1):closedloopssindex
    %---------------------------------------------------------------------
    % arterioleODE for Delta t
    %---------------------------------------------------------------------
    [t,u] = ode15s(@arterioleODE,tspan,u(end,:),[],Pa(pindex),Qnum,Pdrop,...
        multiplier,capdensity,CMD(tt)/1000,td,Pac,Pend,...
        MyoOn_dist,TGFOn_dist,ndistc,mu_dist,Leff_dist,Qdistc,...
        Pdistc,Params_dist,CLc,ta);
    
    diameter_distal(tt) = u(end,1);
    activation_distal(tt) = u(end,2);
    
    Ddistnewi = diameter_distal(end);

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

    Adisti = activation_distal(end);

    SSDTtPAQdisti(pindex,1) = Ddistnewi;
    SSDTtPAQdisti(pindex,2) = Tsigdisti;
    SSDTtPAQdisti(pindex,3) = Pdisti;
    SSDTtPAQdisti(pindex,4) = Adisti;
    SSDTtPAQdisti(pindex,5) = Qdisti;

    Activation(tt) = Adisti;
    perfusioni(tt) = SSDTtPAQdisti(pindex,5);
    diameter_finali(tt) = SSDTtPAQdisti(pindex,1);

    QC(tt) = perfusioni(tt);
    Fi(tt) = alpha*beta*QC(tt);

    [ti,C2i] = ode45(@LoopHenlePDE,tspan,C_0,[],N,dx,A,Vmax,Km,p,Ce,r,Fi(tt));
    a=size(C2i,1);
%Cl- concentration feedback to afferent artery through macula densa (MD) and
%juxtaglomerular apparatus (JGA)    
    CMD(int_tau0+tt) = C2i(a,N);
%Concatenate the solution vectors    
    Csoln(tt,:) = C2i(a,:);
    %update times & initial condition for next call to ODE solver
    t0 = tf;
    tf = t0 + tinc; 
    tspan = [t0 tf];
    C_0 = C2i(a,:);
end
%---------------------------------------------------------------------
% Third time period
% Constant pulse: Specified Fpulse as a constant independent of arterioleODE
% Flow in this period should be constant, while the CMD value should
% change.
%---------------------------------------------------------------------
for i = (closedloopssindex+1):pulseindex
    [ti,C3i] = ode45(@LoopHenlePDE,tspan,C_0,[],N,dx,A,Vmax,Km,p,Ce,r,Fpulse);
    a=size(C3i,1);
    %Cl- concentration feedback to afferent artery through macula densa (MD) and
    %juxtaglomerular apparatus (JGA)    
    CMD(int_tau0+i) = C3i(a,N);
    %Concatenate the solution vectors  
    Csoln(i,:) = C3i(a,:);
    %update times & initial condition for next call to ODE solver
    t0 = tf;
    tf = t0 + tinc; 
    tspan = [t0 tf];
    C_0 = Csoln(i,:);

    diameter_distal(i) = diameter_distal(i-1);
    activation_distal(i) = activation_distal(i-1);

%     diameter_distal = [diameter_distal; u(:,1)];
%     activation_distal = [activation_distal; u(:,2)];
%     
    Ddistnewi = diameter_distal(end);

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

    Adisti = activation_distal(end);

    SSDTtPAQdisti(pindex,1) = Ddistnewi;
    SSDTtPAQdisti(pindex,2) = Tsigdisti;
    SSDTtPAQdisti(pindex,3) = Pdisti;
    SSDTtPAQdisti(pindex,4) = Adisti;
    SSDTtPAQdisti(pindex,5) = Qdisti;

    Activation(i) = Adisti;
    perfusioni(i) = SSDTtPAQdisti(pindex,5);
    diameter_finali(i) = SSDTtPAQdisti(pindex,1);   
    QC(i) = perfusioni(i);
    Fi(i) = alpha*beta*QC(i);
end
%initial condition
Y0 = C_0;
%---------------------------------------------------------------------
% Fourth time period
% Closed-loop dynamics
%---------------------------------------------------------------------
for tt = (pulseindex+1):NT %time index
    %---------------------------------------------------------------------
    % arterioleODE for Delta t
    %---------------------------------------------------------------------
    [t,u] = ode15s(@arterioleODE,tspan,u(end,:),[],Pa(pindex),Qnum,Pdrop,...
        multiplier,capdensity,CMD(tt)/1000,td,Pac,Pend,...
        MyoOn_dist,TGFOn_dist,ndistc,mu_dist,Leff_dist,Qdistc,...
        Pdistc,Params_dist,CLc,ta);
    diameter_distal(tt) = u(end,1);
    activation_distal(tt) = u(end,2);
    
    Ddistnewi = diameter_distal(end);

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

    Adisti = activation_distal(end);

    SSDTtPAQdisti(pindex,1) = Ddistnewi;
    SSDTtPAQdisti(pindex,2) = Tsigdisti;
    SSDTtPAQdisti(pindex,3) = Pdisti;
    SSDTtPAQdisti(pindex,4) = Adisti;
    SSDTtPAQdisti(pindex,5) = Qdisti;

    Activation(tt) = Adisti;
    perfusioni(tt) = SSDTtPAQdisti(pindex,5);
    diameter_finali(tt) = SSDTtPAQdisti(pindex,1);

    QC(tt) = perfusioni(tt);
    Fi(tt) = alpha*beta*QC(tt);
    %---------------------------------------------------------------------
    % LoopHenlePDE for Delta t
    %---------------------------------------------------------------------
    [ti,Yi] = ode45(@LoopHenlePDE,tspan,Y0,[],N,dx,A,Vmax,Km,p,Ce,r,Fi(tt));
    
    %---------------------------------------------------------------------
    % finalize this time step
    %---------------------------------------------------------------------
    a=size(Yi,1);

%Cl- concentration feedback to afferent artery through macula densa (MD) and
%juxtaglomerular apparatus (JGA)    
    CMD(int_tau0+tt) = Yi(a,N);
%Concatenate the solution vectors    
    Csoln(tt,:) = Yi(a,:);
    %update times & initial condition for next call to ODE solver
    t0 = tf;
    tf = t0 + tinc; 
    tspan = [t0 tf];
    Y0 = Yi(a,:);
%     u0 = u(a,:);
end
%  save DAsolver2.mat

figure(1)
plot(x/L,Csoln(1,:).*1000,'g');
hold on
plot(x/L,Csoln(NT,:).*1000,'r');
hold off
legend('t = 0s',['t = ' num2str(tfinal) 's'],'Location','SouthWest')
xlabel('x/L')
ylabel('[ Cl^- ] (mM)')

figure(2)
plot(time,Csoln(:,1).*1000,'g')
hold on
plot(time,Csoln(:,(N-1)/2+1).*1000,'b')
plot(time,Csoln(:,N).*1000,'r')
hold off
xlabel('t')
ylabel('[ Cl^- ] (mM)')
legend('x/L = 0',['x/L = ' num2str(x((N-1)/2+1)/L)],['x/L = ' num2str(x(N)/L)])

figure(3)
plot(x/L,Ce*1000)
xlabel('x/L')
ylabel('C_e (mM)')
% check: Ce(N) == CeN

% %Cl- concentration at final time at x = L
% ['C(L,tfinal) = ' num2str(C(NT,N)*1000) ' mM']
%Cl- concentration at final time at the macula densa with time delay in JGA
['C_{MD} = C(L,tfinal-tau) = ' num2str(CMD(NThistory).*1000) ' mM']
['tfinal-tau = ' num2str(tfinal-tau0) 's']

Qvector = ones(1,NT);
for ii=1:openloopssindex
    Qvector(ii) = Q;
end
for ii = (openloopssindex+1):closedloopssindex
    Qvector(ii)=QC(ii)*beta*60*10^6;
end
for ii = (closedloopssindex+1):pulseindex
    Qvector(ii) = Qpulse;
end
for ii = (pulseindex+1):NT
    Qvector(ii) = QC(ii)*beta*60*10^6;
end

figure(5)
plot(time,Qvector)
xlabel('t (s)')
ylabel('Q (nL/min)')
% axis([0 tfinal 20 40])
figure(6)
plot(time,CMD((int_tau0+1):(int_tau0+NT)).*1000)
xlabel('t (s)')
ylabel('CMD (mM)')
% axis([0 tfinal 20 40])

figure(20)
QC(1)=Qdistc;
plot(time,QC*60*10^6)
xlabel('t (s)')
ylabel('Q_{AA} (nL/min)')

figure(21)
plot(time,Activation)
xlabel('t (s)')
ylabel('Activation')

figure(22)
plot(time,diameter_finali*1e4)
xlabel('t (s)')
ylabel('diameter (\mu m)')

save(['output_tau' num2str(tau0) '_tfinal' num2str(tfinal) 's_td' num2str(td) '_ta' num2str(ta) '_Qpulse' num2str(Qpulse) '.mat'])
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

