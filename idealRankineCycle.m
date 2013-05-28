% Sree Prasanna Rajagopal, 
% [Mechanical Engineering Department, IIT Guwahati] February 2013
% XSteam.m by Magnus Holmgren, www.x-eng.com 
% used for all property calculation purposes

% IDEAL RANKINE CYCLE

function [ x2,eff ] = idealRankineCycle(P1,P2,T1)
% ideal rankine cycle analysis

% subscripts and associated states
% 1 - turbine inlet
% 2 - turbine outlet
% 3 - condensor outlet
% 4 - boiler inlet

% P1 - condenser pressure
% P2 - boiler pressure

% turbine inlet
h1 = XSteam('h_pT',P2,T1);
s1 = XSteam('s_ph',P2,h1);

% entropy remains the same through turbine ( ideal )
s2=s1;
h2 = XSteam('h_ps',P1,s2);
x2 = XSteam('x_ps',P1,s2);

if (x2 < 0.88 )
    disp('Steam quality less than 88%.');
    fprintf('\n\n\n');
end

% condenser outlet
T3 = XSteam('Tsat_P',P1);
h3 = XSteam('hL_P',P1);
v3 = XSteam('vL_p',P1);

% calculating h4 from pump work
h4 = h3 + v3*(P2-P1)*1e2;

% calculating pump work, turbine work and heat input
Wt = h1 - h2;
Wp = h4 - h3;
Qin = h1 - h4;

% calculating efficiency
eff = (Wt - Wp)/Qin;


    