function [V1,V2] = TcellDoublePump(Iinj, dt, sim20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two-compartment leech T cell model with transient Na, 
% high-voltage-activated K and slow M-type K current and a Na+/K+ pump.
% Written by Kevin Sandbote (kevin.sandbote@uni-oldenburg.de)
% Version 0.1
% March 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +++ input
%  Iinj : vector for injected current input into somatic compartment [pA]
%  dt : time step [ms]
%  pA: vector holding a value for: [EK,EL1,EL2,GL1,GL2,GM,Gc,ENa,gNa,gK]
%         in that order!
%  sim20: flag for the simulation of the 20th trial, inh. current is
%  applied if set to true
% +++ output
%  V1 : membrane potential of the soma [mV]
%  V2 : membrane potential of the spike-initiation zone [mV]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note: somatic compartment has only leak and no Na, K, or KM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% assigning membrane parameters
% Compartment sizes
% dSoma = 40; %[um] diameter of soma
% S1 = pi*dSoma^2; % [um^2] surface area of the somatic membrane
% rSIZ = 10/2; %[um] radius SIZ
% hSIZ = 10; %[um] hight SIZ
% S2 = 2*rSIZ*pi*(rSIZ+hSIZ); % [um^2] surface area of the SIZ
S1 = 5000; %[um^2] surface area of the somatic membrane
S2 = 500; % [um^2] surface area of the SIZ

% capacitances
Cm = 1.0;  % [uF/cm^2] membrane capacitance density 
c2 = Cm * S2 * 1e-8; % [uF] membrane capacitance of spike-initiation zone
c1 = Cm * S1 * 1e-8; % [uF] membrane capacitance of somatic compartment

% conductance densities
GN = 660; % [mS/cm^2] transient Na conductance density
GK = 24; % [mS/cm^2] high-voltage activated K conductance density
GM = 16; % [mS/cm^2] M-type K conductance density
GL1 = 0.35; % [mS/cm^2] leak conductance density (spike-initiation zone)
GL2 = 0.25;  % [mS/cm^2] leak conductance density (soma)

% reversdal potentials
EN = 40; % [mV] Na reversal potential 
EK = -70; % [mV] K reversal potential 
EL1 = -15; % [mV] leak reversal potential (spike-initiation zone)
EL2 = 0; % [mV] leak reversal potential (soma)

% Na/K-pump current parameters
Imax = 800e-6; % [uA] max pump current (=800[pA])
kIN = 10 * 0.06e-6; % conversion factor from INa into [Na] [mM/pA.ms]
kIP =  2 * 0.06e-6; % conversion factor from Ipump into [Na] [mM/pA.ms]

% conductances 
gN = GN * S2 * 1e-8; % [mS] Na conductance of spike-initiation zone
gK = GK * S2 * 1e-8; % [mS] K conductance of spike-initiation zone
gM = GM * S2 * 1e-8; % [mS] M-type K conductance of spike-initiation zone
gL2= GL2* S2 * 1e-8; % [mS] leak conductance of spike-initiation zone
gL1= GL1 * S1 * 1e-8; % [mS] leak conductance of somatic compartment
gC = 90 * 1e-6; % [mS] conductance between the two compartments

% activation and inactivation of KM
v0 = 37;
vk = 4;
TM = 350;

% activation and inactivation of Na
v0Na = 20;
vNa = 8;
vTNa = 16;
TNa = 0.75;

v0H = 36;
vH = 5;
vTH = 10;
TH = 7.5;

%% data vectors and initial values
Ntotal = length(Iinj); % number of data points
V1 = zeros(1,Ntotal); % [mV] membrane potential of somatic compartment
V2 = zeros(1,Ntotal); % [mV] membrane potential of spike-initiation zone


% estimate initial potential value 
cinit = 0.0; % [mM] initial Na concentration (measured from rest) 
vtest = -80:0.01:-20;
IN0 = gN * infM(vtest,v0Na,vNa).^4 .* infH(vtest,v0H,vH) .* (EN-vtest); % [uA]
IK0 = gK * infN(vtest).^2                .* (EK-vtest); % [uA]
IL0 = gL1                                .* (EL1-vtest); % [uA]
IM0 = gM * infZ(vtest,v0,vk).^2                .* (EK-vtest); % [uA]
IP0 = -Imax * pumpNa(cinit); % [uA]
Itotal = IN0 + IK0 + IL0 + IM0 + IP0 + Iinj(1)*1e-6; % [uA]
[~,idx] = min(abs(Itotal)); % get the point where all currents balance
Vinit = vtest(idx); 


% if the 20th trial is to be simulated, the current needed to put the
% membrane potential to the right value is calculated here
if sim20 ==1
    
    vtestInh = -48.8;
    
    IN0 = gN * infM(vtestInh,v0Na,vNa).^4 .* infH(vtestInh,v0H,vH) .* (EN-vtestInh); % [uA]
    IK0 = gK * infN(vtestInh).^2                .* (EK-vtestInh); % [uA]
    IL0 = gL1                                .* (EL1-vtestInh); % [uA]
    IM0 = gM * infZ(vtestInh,v0,vk).^2                .* (EK-vtestInh); % [uA]
    IP0 = -Imax * pumpNa(cinit); % [uA]
    ItotalInh = IN0 + IK0 + IL0 + IM0 + IP0 + Iinj(1)*1e-6; % [uA]
    CurrDiff = ItotalInh - Itotal(idx);
    Iperm = -CurrDiff;
else 
    Iperm = 0;
end


% initial values
V1(1) = Vinit; 
V2(1) = Vinit;
m = infM(Vinit,v0Na,vNa);
h = infH(Vinit,v0H,vH);
n = infN(Vinit);
z = infZ(Vinit,v0,vk);
c = cinit;


%% calculate membrane response step-by-step 
for j=1:Ntotal-1
    
    % ionic currents (soma): g[mS] * V[mV] = I[uA]
    IL1 = gL1 *(EL1-V1(j)); % leak current
    IC  = gC *(V2(j)-V1(j)); % current from compartment 2 to compartment 1

    % ionic currents (spike-initiation zone): g[mS] * V[mV] = I[uA]
    IN = gN * m^4 * h * (EN-V2(j)); % Na current
    IK = gK * n^2     * (EK-V2(j)); % K current
    IM = gM * z^2     * (EK-V2(j)); % KM current
    IL2= gL2          *(EL2-V2(j)); % leak current
    IP = -Imax * pumpNa(c);         % [?A]

    
    % derivatives: I[uA] / C[uF] * dt[ms] = dv[mV]
    dv1_dt = ( Iinj(j)*1e-6 + IL1 + IC + Iperm ) / c1; 
    dv2_dt = ( IN + IK + IM + IL2 + IP - IC  ) / c2; 
    dm_dt = ( infM(V2(j), v0Na, vNa) - m ) / tauM(V2(j), v0Na, vTNa, TNa);
    dh_dt = ( infH(V2(j), v0H, vH) - h ) / tauH(V2(j), v0H, vTH, TH);
    dn_dt = ( infN(V2(j)) - n ) / tauN(V2(j));
    dz_dt = ( infZ(V2(j), v0, vk) - z ) / tauZ(V2(j), v0, vk, TM);
    dc_dt = ( kIN*IN + kIP*3*IP )* 1e6;

    % calculate next step 
    V1(j+1) = V1(j) + dv1_dt * dt; 
    V2(j+1) = V2(j) + dv2_dt * dt; 
    m = m + dm_dt * dt; 
    h = h + dh_dt * dt; 
    n = n + dn_dt * dt; 
    z = z + dz_dt * dt; 
    c = c + dc_dt * dt;

end

%% activation/inactivation functions
function x = infM(v,v0Na,vNa) % steady-state function for Na activation 
  x = 1 ./ ( 1 + exp(-(v+v0Na)/vNa) ); %20=v0 & vNA = 8 
function x = tauM(v,v0Na,vTNa,TNa) % time scale [ms] for Na activation 
  x = TNa*( 0.1 + 2 ./ ( exp(-(v+v0Na)/vTNa) + exp((v+v0Na)/vTNa) ) ); % v0Na = 20 & vTNa = 16& TNa = 0.75
function x = infH(v, v0H, vH) % steady-state function for Na inactivation 
  x = 1 ./ ( 1 + exp( (v+v0H)/vH) ); % v0H = 36 vTH = 5
function x = tauH(v,v0H,vTH,TH) % time scale [ms] for Na inactivation 
  x = TH*( 0.1 + 2 ./ ( exp(-(v+v0H)/vTH) + exp((v+v0H)/vTH) ) ); %v0H = 36 & vTH = 10 & TH = 7.5 
function x = infN(v) % steady-state function for K activation 
  x = 1 ./ ( 1 + exp(-(v+20)/8) ); 
function x = tauN(v) % time scale [ms] for K activation 
  x = 4.0*( 0.1 + 2 ./ ( exp(-(v+20)/16) + exp((v+20)/16) ) ); 
function x = infZ(v,v0,vk) % steady-state function for M-type K activation 
  x = 1 ./ ( 1 + exp(-(v+v0)/vk) ); 
function x = tauZ(v,v0,vk,TM) % time scale [ms] for M-type K activation 
  x = TM * ( 1 + ( 2 ./ ( exp(-(v+v0)/2*vk) + exp((v+v0)/2*vk) ) ) ); 
function x = pumpNa(ci) % Na/K-pump activation
  x = 1 ./ ( 1 + exp(-(ci-18.0)/18.0) ).^3; % ci[mM] 