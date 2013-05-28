% Sree Prasanna Rajagopal, 
% [Mechanical Engineering Department, IIT Guwahati] February 2013
% XSteam.m by Magnus Holmgren, www.x-eng.com 
% used for all property calculation purposes

% IDEAL REHEAT REGENERATIVE WITH ONE OPEN FWH

function [ x4,eff ] = idealReheatRegenerative(T1,P1,P3,P2)

% the fraction of steam extracted is not a variable since it is an ideal
% reheat regenrative cycle being considered and the fraction extracted is
% such that the mixture is a saturated liquid.

% subscripts and associated states
% 1 - turbine inlet
% 2 - turbine outlet
% 3 - reheat out
% 4 - condenser in, second turbine out
% 5 - condenser out, pump in
% 6 - pump out, P3
% 7 - mixture final state 
% 8 - pump exit, P2

% P1 - condenser pressure
% P2 - boiler pressure
% P3 - reheat pressure

% turbine inlet         
h1 = XSteam('h_pT',P2,T1);
s1 = XSteam('s_ph',P2,h1);

% entropy remains the same through turbine ( ideal )
% reheat initial state
s2 = s1;
h2 = XSteam('h_ps',P3,s2);

% reheat final state
% known temperature and pressure
h3 = XSteam('h_pT',P3,T1);
s3 = XSteam('s_pT',P3,T1);

% condenser inlet
s4 = s3;
h4 = XSteam('h_ps',P1,s4);
x4 = XSteam('x_ps',P1,s4);

% Checking steam quality
if (x4 < 0.88 )
    disp('Steam quality less than 88%.');
end

% condenser outlet
T5 = XSteam('Tsat_P',P1);
h5 = XSteam('hL_P',P1);
v5 = XSteam('vL_p',P1);

% pump work to h6
h6 = h5 + v5*(P3-P1)*1e2;

% saturated state at P3
h7 = XSteam('hL_p',P3);
v7 = XSteam('vL_p',P3);

% second pump out, P2
h8 = h7 + v7*(P2-P3)*1e2;

% mixture and fraction calculation according to assumption that the final
% mixture is saturated at P3.
% expression: y(h2 - h7) = (1-y)(h7 - h6)
y = (h7 - h6)/(h2 - h6);

% calculating pump work, turbine work and Qin
Wp = (1-y)*(h6 -h5) + h8 - h7;
Wt = h1- h2 + (1-y)*(h3 - h4);
Qin = h1 - h8 + (1-y)*(h3 - h2);

%calculating efficiency
eff = (Wt - Wp)/Qin;

