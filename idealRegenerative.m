% Sree Prasanna Rajagopal, 
% [Mechanical Engineering Department, IIT Guwahati] February 2013
% XSteam.m by Magnus Holmgren, www.x-eng.com 
% used for all property calculation purposes

% IDEAL REGENERATIVE RANKINE CYCLE WITH ONE OPEN FWH

function [ x3,eff] = idealRegenerative(T1,P1,P3,P2)

% ideal reheat rankine cycle analysis

% subscripts and associated states
% 1 - turbine inlet
% 2 - turbine outlet
% 3 - condenser in
% 4 - condenser out, pump in
% 5 - pump out, mix start
% 6 - mixture saturated state, pump in
% 7 - final pump out, boiler in

% P1 - condenser pressure
% P2 - boiler pressure
% P3 - extraction pressure

% turbine inlet           
h1 = XSteam('h_pT',P2,T1);
s1 = XSteam('s_ph',P2,h1);


% entropy remains the same through turbine ( ideal )
% extraction pressure P1
s2 = s1;
h2 = XSteam('h_ps',P3,s2);

%condenser inlet
% known temperature and pressure
s3 = s1;
h3 = XSteam('h_ps',P1,s3);
s3 = XSteam('s_ph',P1,h3);
x3 = XSteam('x_ps',P1,s3);

% Checking steam quality
if (x3 < 0.88 )
    disp('Steam quality less than 88%.');
end

% condenser outlet
h4 = XSteam('hL_p',P1);
v4 = XSteam('vL_p',P1);

% calculating h6 from pump work
h5 = h4 + v4*(P3-P1)*1e2;

% state at h6 - saturated
h6 = XSteam('hL_p',P3);
v6 = XSteam('vL_p',P3);

% h7 after pump work
h7 = h6 + v6*(P2-P3)*1e2;
% mixture enthalpy expression: (1-y)h5 + (y)h2 = h6
% simplifying y = (h5 - h6)/(h5 - h2)
y = (h5 - h6)/(h5 - h2);

% calculating pump work, turbine work and Qin
Wp = (1-y)*(h5 - h4) + h7 - h6;
Wt = h1 - h2 + (1-y)*(h2 - h3);
Qin = h1 - h7;

%calculating efficiency
eff = (Wt - Wp)/Qin;