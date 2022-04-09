%% Part 1
%
clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')

meshsizex = 6;
meshsizey = 6;



vinlow = 0.1;
vinhigh = 10;

% Components

Cap = 0.25;
R1 = 1;
R2 = 2;
L = 0.2;



R3 = 10;
alpha = 100;
R4 = 0.1;
Ro = 1000;
omg = 10;

% C Matrix
%C=zeros(6,6);
C = zeros(meshsizex,meshsizey);
C(2,1) = -Cap;
C(2,2) = Cap;
C(6,6) = L;

% G Matrix
G = zeros (meshsizex, meshsizey);
G(1,1) = 1;
G(2,1) = -1/R1;
G(2,2) = (1/R1) + (1/R2);
G(2,6) = -1;
G(3,3) = 1/R3;
G(3,6) = 1;
G(4,3) = -alpha/R3;
G(4,4) = 1;
G(5,4) = -R4;
G(5,5) = R4 - (1/Ro);
G(6,2) = 1;
G(6,3) = -1;
%%
% The  C and G matrix
C;
G;
%%
F=zeros(6,1);
Vv1 = zeros(10,1);
Vv2 = zeros(100,1);
gain = zeros(10,1);
vin = linspace(-10,10,100);
for i=1:100
   %F(6) = i;
   F(1) = vin(i);
   V=G\F;
   
   Vv1(i) = V(5);
   Vv2(i) = V(3);
end
figure (1)
plot(vin,Vv1)
title('DC SWEEP');
xlabel('voltage sweep');
ylabel(' output Voltage');
hold on
plot(vin,Vv2)

%THE AC 
Vv1 = zeros(10,1);
Vv2 = zeros(100,1);

step=100;
omg = 2*pi*linspace(0,20,step);

for i=1:100
omg =2*pi*linspace(0,20,step);
F(1) = vin(i);

s=i*omg(i);
    x = G + (s.*C);
    V = x\F;
 
     Vv1(i) = V(5);
    Vv2(i) = V(3);
     gain(i) = 20 * log(abs(Vv1(i))/abs(V(1)));
end
figure (2)
plot (omg,abs(Vv1))
title(' 0 to 100Hz');
hold on;
figure (2)
plot(omg,abs(Vv2))
legend('Vo', 'V3');
xlabel('w');
ylabel('V');
figure(3)
plot(omg, gain);
title('Gain vs omega');
xlabel('omega');
ylabel('gain');

s = 1i*pi;

%C case plot the gain as function of random perturbations on C using a normal distribution
for i = 1:1000
    randomcap = Cap + 0.05*randn(); % random perturbations on C
    C(2, 1) = randomcap; 
    C(2, 2) = -randomcap;
    C(6, 6) = L;

    H = G +(s.*C);
    V = H\F;
    Vv1(i) = abs(V(5));
    gain(i) = 20*log10(abs(Vv1(i))/abs(V(1)));
end

figure(4)
histogram(gain)
title('Gain distribution');
grid on


% C dV/dt + GV = F
% C*(V_i - V_(i-1))/delta t + G*V_i = F
% (C/delta t + G) * V_i + C/delta t*V_(i-1) = F
% C/delta t * V_(i-1) - F = (C/delta t + G)*V_i
% (C/delta t + G)  \ (C/delta t * V_(i-1) - F) = V_i  <-


Simtime = 1;
numSteps = 1000;
deltaT = Simtime/numSteps;
% H vector

% C Matrix
%C=zeros(6,6);
C = zeros(meshsizex,meshsizey);
C(2,1) = -Cap;
C(2,2) = Cap;

C(6,6) = L;

H = (C./deltaT) + G;

 V = zeros(6,1);
 Vp = V;


% F vector
F = zeros(1,6);

Vv1 = zeros(10,1);
Vv2 = zeros(100,1);
% Time simulation - step function

newtime = linspace(1,numSteps,numSteps);

for i = 1:numSteps
    
    % F vector
    
    if (i == 20)
        F(1) = 1;
    end
    newtime(i) = i*deltaT;
    V = H\(((C * Vp)./deltaT) + F');

    Vv1(i) = V(1);
    Vv2(i) = V(5);
    
    Vp = V;
end

 figure(5)
   plot(newtime,Vv1,'-r');
    hold on
   plot(newtime,Vv2,'g');
pause(0.01);
legend('Vin', 'Vout');
title('Time Simulation vs Step function');
xlabel('Time');
ylabel('Voltage');
figure(6)
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv2)))),'-g');
hold on 
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv1)))),'-b');
legend('Vout', 'Vin');
xlabel('Frequency');
ylabel('Voltage (dB)');
title(' Frequency vs voltage');


% Time simulation - sinusoidal function

Simtime = 1;
numSteps = 1000;
deltaT = Simtime/numSteps;
% H vector
H = (C./deltaT) + G;

Vv1 = zeros(10,1);
Vv2 = zeros(100,1);
Vp = zeros(6,1);

newtime = linspace(1,numSteps,numSteps);
for i = 1:numSteps
    
    % F vector
    
    F(1) = sin(2 * pi * (1/0.03) * newtime(i) * deltaT);
    
    V = H\(((C * Vp)./deltaT) + F');
  
    
    
   Vv1(i) = V(1);
    Vv2(i) = V(5);
    Vp = V;
end

 figure(7)
    plot(newtime,Vv1,'r');
    hold on
   plot(newtime,Vv2,'g');
pause(0.01);
legend('Vin', 'Vout');
title('Time Simulation vs Step function');
xlabel('Time');
ylabel('Voltage');
legend('Vin', 'Vout');
title('The Time Simulation vs The Sinusoidal function');
xlabel('Time');
ylabel('Voltage');

figure(8)
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv2)))),'g');
hold on 
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv1)))),'b');
legend('Vout', 'Vin');
xlabel('Frequency');
ylabel('Voltage (dB)');
title(' Frequency vs voltage');


% plot(ti,g)
% g = exp(-(ti-0.3).^2/0.1^2)
%A guassian pulse with a magnitude of 1, std dev. of 0.03s and a delay of 0.06s.

Simtime = 1;
numSteps = 1000;
deltaT = Simtime/numSteps;
% H vector
H = (C./deltaT) + G;

Vv1 = zeros(10,1);
Vv2 = zeros(100,1);
Vp = zeros(6,1);

gp=sin(2*pi*0.03)*exp(-(newtime-0.06-0.5).^2/(2*0.03^2));
newtime = linspace(1,numSteps,numSteps);
for i = 1:numSteps
    
    % F vector
    if (i == 40)
        F(1) = gp(i);
    end
    
   
    
    V = H\(((C * Vp)./deltaT) + F');
  
    
    
   Vv1(i) = V(1);
    Vv2(i) = V(5);
    Vp = V;
end

 figure(9)
    plot(newtime,Vv1,'r');
    hold on
   plot(newtime,Vv2,'g');
pause(0.01);
legend('Vin', 'Vout');
title('Time Simulation vs guassian pulse');
xlabel('Time');
ylabel('Voltage');
legend('Vin', 'Vout');
title('The Time Simulation vs The guassian pulse');
xlabel('frequency');
ylabel('Voltage');

figure(10)
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv2)))),'g');
hold on 
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv1)))),'b');
legend('Vout', 'Vin');
xlabel('Frequency');
ylabel('Voltage (dB)');
title(' Frequency vs voltage');





