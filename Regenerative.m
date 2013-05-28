% Sree Prasanna Rajagopal, 
% [Mechanical Engineering Department, IIT Guwahati] February 2013
% XSteam.m by Magnus Holmgren, www.x-eng.com 
% used for all property calculation purposes

% IDEAL REGENERATIVE RANKINE WITH ONE OPEN AND ONE CLOSED FWH

function [ x4, eff] = Regenerative(T1,P1,P4,P3,P2)

% the fraction of steam extracted is not a variable since it is an ideal
% reheat regenrative cycle being considered and the fraction extracted is
% such that the mixture is a saturated liquid and calculated from the above
% condition

% subscripts and associated states
% 1 - turbine inlet
% 2 - closed fwh extraction
% 3 - open fwh extraction
% 4 - condenser in
% 5 - condenser out
% 6 - pump out, P3
% 7 - mixture, open fwh out
% 8 - pump to boiler pressure, P2
% 9 - closed fwh exit, saturated liquid
% 11 - closed feed heater liquid pumped to boiler pressure 
% 10 - feed water after heating by closed fwh. Ideally h10 = h9 
% 12 - after mixing at boiler pressure

% P1 - condenser pressure
% P2 - boiler pressure
% P3 - extraction pressure for closed fwh
% P4 - extraction pressure for open fwh

% turbine inlet         
h1 = XSteam('h_pT',P2,T1);
s1 = XSteam('s_ph',P2,h1);

% entropy remains the same through turbine ( ideal )
% extraction point for closed fwh. 'z' fraction extracted,
s2 = s1;
h2 = XSteam('h_ps',P3,s2);

% extraction point for open fwh. 'y' fraction extracted and '1-y-z' expands
% to condenser pressure
s3 = s1;
h3 = XSteam('h_ps',P4,s3);
x3 = XSteam('x_ps',P4,s3);

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
h6 = h5 + v5*(P4-P1)*1e2;

% 7 - saturated liquid state at P4
h7 = XSteam('hL_p',P4);
v7 = XSteam('vL_p',P4);

% enthalpy mixing in open fwh
% expression: (1-y-z)h6 + yh3 = (1-z)h7
% cannot be calculated now. closed fwh fraction 'z' needed
exp1 = '(1-y-z)*h6 + y*h3 = (1-z)*h7)';

% second pump to boiler pressure
h8 = h7 + v7*(P2 - P4)*1e2;

% heating by closed fwh
% ideal assumption that final feed water enthalpy is equal to final
% extracted water enthalpy
% expression: (h8)*(1-z) + (h2)*z = h10*(1-z) + h9*(z)
h9 = XSteam('hL_p',P3);
v9 = XSteam('vL_p',P3);
h10 = h9;
z = (h10 - h8)/(h2 - h9 + h10 - h8);

% calculating y
y = (z*(h6 - h7) + h7 - h6)/(h3 - h6);

% extracted water pumped to boiler pressure
h11 = h9 + v9*(P2 - P3)*1e2;

% mixing at boiler pressure
% expression: (1-z)h12 + zh10 = 1(h13)
h12 = h10*(1 - z) + z*h11;


% calculating pump work, turbine work and Qin
Wp = (1-y-z)*(h6 - h5) + (1-z)*(h8 - h7) + z*(h11 - h9);
Wt = h1 - h2 + (1-z)*(h2 - h3) + (1-y-z)*(h3 - h4);
Qin = h1 - h12;

%calculating efficiency
eff = (Wt - Wp)/Qin;

