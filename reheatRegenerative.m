% Sree Prasanna Rajagopal, 
% [Mechanical Engineering Department, IIT Guwahati] February 2013
% XSteam.m by Magnus Holmgren, www.x-eng.com 
% used for all property calculation purposes

% IDEAL REHEAT REGENERATIVE RANKINE WITH ONE OPEN AND ONE CLOSED FWH

function [ x5, eff] = reheatRegenerative(T1,P1,P4,P3,P2)

% the fraction of steam extracted is not a variable since it is an ideal
% reheat regenrative cycle being considered and the fraction extracted is
% such that the mixture is a saturated liquid and calculated from the above
% condition

% subscripts and associated states
% 1 - turbine inlet
% 2 - turbine outlet
% 3 - reheat out
% 4 - open fwh extraction point
% 5 - condenser in
% 6 - condenser out, pump in
% 7 - pump out, P3
% 8 - mixture final state after open fwh
% 9 - pump exit, P2
% 11 - closed fwh exit
% 10 - pump from state 11 to P2
% 12 - final state after heating by closed fwh; ideally h12 = h11
% 13 - boiler inlet after mixing h12 and h10

% P1 - condenser pressure
% P2 - boiler pressure
% P3 - reheat pressure, extraction pressure for closed fwh
% P4 - extraction pressure for open fwh

% turbine inlet         
h1 = XSteam('h_pT',P2,T1);
s1 = XSteam('s_ph',P2,h1);

% entropy remains the same through turbine ( ideal )
% reheat initial state.
% extraction point for closed fwh. 'z' fraction extracted, '1-z' gets
% reheated
s2 = s1;
h2 = XSteam('h_ps',P3,s2);

% reheat final state
% known temperature and pressure
h3 = XSteam('h_pT',P3,T1);
s3 = XSteam('s_pT',P3,T1);

% extraction point for open fwh. 'y' fraction extracted and '1-y-z' expands
% to condenser pressure
s4 = s3;
h4 = XSteam('h_ps',P4,s4);
x4 = XSteam('x_ps',P4,s4);

% condenser inlet
s5 = s3;
h5 = XSteam('h_ps',P1,s5);
x5 = XSteam('x_ps',P1,s5);

% Checking steam quality
if (x5 < 0.88 )
    disp('Steam quality less than 88%.');
end

% condenser outlet
T6 = XSteam('Tsat_P',P1);
h6 = XSteam('hL_P',P1);
v6 = XSteam('vL_p',P1);

% pump work to h6
h7 = h6 + v6*(P4-P1)*1e2;

% saturated state at P4
h8 = XSteam('hL_p',P4);
v8 = XSteam('vL_p',P4);

% enthalpy mixing in open fwh
% expression: (1-y-z)h7 + yh4 = (1-z)h8
% cannot be calculated now. closed fwh fraction 'z' needed
exp1 = '(1-y-z)*h7 + y*h4 = (1-z)*h8)';

% second pump
h9 = h8 + v8*(P2 - P4)*1e2;

% heating by closed fwh
% ideal assumption that final feed water enthalpy is equal to final
% extracted water enthalpy
% expression: (h9)*(1-z) + (h2)*z = h12*(1-z) + h11*(z)
h11 = XSteam('hL_p',P3);
v11 = XSteam('vL_p',P3);
h12 = h11;
z = (h12 - h9)/(h2 - h11 + h12 - h9);

% calculating y
y = (z*(h7 - h8) + h8 - h7)/(h4 - h7);

% extracted water pumped to boiler pressure
h10 = h11 + v11*(P2 - P3)*1e2;

% mixing at boiler pressure
% expression: (1-z)h12 + zh10 = 1(h13)
h13 = h12*(1 - z) + z*h10;


% calculating pump work, turbine work and Qin
Wp = (1-y-z)*(h7 - h6) + (1-z)*(h9 - h8) + z*(h10 - h11);
Wt = h1 - h2 + (1-z)*(h3 - h4) + (1-y-z)*(h4 - h5);
Qin = h1 - h13 + (1-z)*(h3 - h2);

%calculating efficiency
eff = (Wt - Wp)/Qin;

