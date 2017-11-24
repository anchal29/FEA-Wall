% SDOF system tests
clc;
clear variables;
P = load('Elcentro.txt');
g = 9.81;
m = 2000;%Let
k = 2*10^6;%Let
% Using the wall data
m = 6375;
% m = 2.5*10^9;
k = 892275.971; 
k = 2*10^3;
zeta = 0;
wn = sqrt(k/m);
c = 2*m*wn*zeta;
%c = 0.05;
for i = 1:length(P)
    P(i,2) = (P(i,2)*m*g);
end

choice = input('Select the case:\n1.Average Acceleration Method, or\n2.Linear Acceleration Method\nYour choice(1/2):        ');
if(choice == 1)
    gamma = 0.5;
    beta = 0.25;
else
    gamma = 0.5;
    beta = 1/6;
end
timeStep = 0.02;

kcap = k + (gamma/(beta*timeStep))*c + m/(beta*timeStep^2);
a = m/(beta*timeStep)  + (gamma*c / beta);
b = m/(beta*2)+timeStep*((gamma/(2*beta))-1);
U(1) = 0;
V(1) = 0;
A(1) = P(1,2)/m;
for i = 1:1:length(P)-1
   
    deltaP(i)  = P(i+1,2) - P(i,2); 
    deltaPcap(i) = deltaP(i) + (a*V(i)) + (b*A(i));
    deltaU(i) = deltaPcap(i)  / kcap;
    deltaV(i) = (gamma*deltaU(i)/(beta*timeStep))-(gamma*V(i)/beta)+(timeStep*A(i)*(1-(gamma/(2*beta))));
    deltaA(i) = (deltaU(i)/(beta*(timeStep^2)))-(V(i)/(beta*timeStep))-(A(i)/(2*beta));
    U(i+1) = U(i) + deltaU(i);
    V(i+1) = V(i) + deltaV(i);
    A(i+1) = A(i) + deltaA(i);
end
t= 0:0.02:0.02*length(P)-0.02;
% figure;
plot(t,U, 'linewidth', 1.3);
hold on;
%1/(pi*sqrt(2)*(sqrt(1-(2*beta))));

%% CDM
disp('Central difference method');

kcap = (m/(timeStep)^2) + (c/(2*timeStep));
a = m/(timeStep^2)  - (c /(2* timeStep));
b = k -(2*m/timeStep^2);



Uminus1 = P(1,2)*timeStep^2/m;
pcap(1) = P(1,2) - a*Uminus1;
U(1) = 0;
U(2) = pcap(1)/kcap;

for i=2:1:length(P)
    pcap(i) = P(i,2) - a*U(i-1) - b*U(i);
    U(i+1) = pcap(i)/kcap;
end

%%
plot(0:timeStep:timeStep*length(P),U, 'linewidth', 1.3);
