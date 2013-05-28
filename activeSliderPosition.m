function [ outputPosition ] = activeSliderPosition(V) 

% define material properties
syms t;
syms youngs inertia linkL k

youngs = 10e9;
inertia = 10e-3;
linkL = 0.5;
k = 10e12;

% find transfer function from relation: G(s) = k^2*l/(6E^2*I^2)

syms transferFunction;
transferFunction = k^2*linkL/(6*youngs^2*inertia^2);

outputPosition = transferFunction*V^2;
% where V is the voltage applied

ezplot(outputPosition);

