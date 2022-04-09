

clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')

meshsizex = 6;
meshsizey = 6;

% Voltage Range

Vmin = 0.1;
Vmax = 10;

% Components

Cap = 0.25;
R1 = 1;
R2 = 2;
L = 0.2;
% R3 = R3finder(Vmin,Vmax,20);
R3 = 10;
alpha = 100;
R4 = 0.1;
Ro = 1000;
omg = 10;

% Noise components
In = 0.001;
Cn =  0.00001;

% C Matrix
C = zeros(meshsizex,meshsizey);
C(2,1) = -Cap;
C(2,2) = Cap;
C(3,3) = Cn;
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
G(6,3) = -1;%%
% (a) Updated C and G matrices
C;
G;

Simtime = 1;
numSteps = 1000;
deltaT = Simtime/numSteps;
% H vector
H = (C./deltaT) + G;

 V = zeros(6,1);
 Vp = V;

% F vector
F = zeros(1,6);
%F(3)=In;
 F(3) = In*randn();

Vv1 = zeros(10,1);
Vv2 = zeros(100,1);
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

 figure(1)
  plot(newtime,Vv1,'-r');
    hold on
   plot(newtime,Vv2,'-b');
pause(0.01);
legend('Vin', 'Vout');
title('Time Simulation vs  The Vout noise source');
figure(2)
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv2)))),'-g');
hold on 
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv1)))),'-b');
legend('Vout', 'Vin');
xlabel('Frequency');
ylabel('Voltage (dB)');
title(' Frequency vs voltage');


%vary Cn to see how the bandwidth changes C(6,6) = 1e-10
C = zeros(meshsizex,meshsizey);
C(2,1) = -Cap;
C(2,2) = Cap;
C(6,6) = 1e-20;
Simtime = 1;
numSteps = 1000;
deltaT = Simtime/numSteps;
% H vector
H = (C./deltaT) + G;

 V = zeros(6,1);
 Vp = V;

% F vector
F = zeros(1,6);
F(3)=In;
 %F(3) = In*randn();

Vv1 = zeros(10,1);
Vv2 = zeros(100,1);
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

 figure(3)
  plot(newtime,Vv1,'-r');
    hold on
   plot(newtime,Vv2,'-b');
pause(0.01);
legend('Vin', 'Vout');
xlabel('Frequency');
ylabel('Voltage (dB)');
title('Time Simulation vs  The Vout noise source');
figure(4)
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv2)))),'-g');
hold on 
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv1)))),'-b');
legend('Vout', 'Vin');
xlabel('Frequency');
ylabel('Voltage (dB)');
title(' Frequency vs voltage');

%vary Cn to see how the bandwidth changes C(6,6) = 1e-7

C(6,6) = 1e-7;
Simtime = 1;
numSteps = 1000;
deltaT = Simtime/numSteps;
% H vector
H = (C./deltaT) + G;

 V = zeros(6,1);
 Vp = V;

% F vector
F = zeros(1,6);
F(3)=In;
 %F(3) = In*randn();

Vv1 = zeros(10,1);
Vv2 = zeros(100,1);
newtime = linspace(1,numSteps,numSteps);

for i = 1:numSteps
    
    % F vector
    
    if (i == 40)
        F(1) = 1;
    end
    newtime(i) = i*deltaT;
    V = H\(((C * Vp)./deltaT) + F');

    Vv1(i) = V(1);
    Vv2(i) = V(5);
    
    Vp = V;
    
end

 figure(5)
  plot(newtime,Vv1,'r');
    hold on
   plot(newtime,Vv2,'b');
pause(0.01);
legend('Vin', 'Vout');
xlabel('Frequency');
ylabel('Voltage (dB)');
title('Time Simulation vs  The Vout noise source');
figure(6)
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv2)))),'g');
hold on 
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv1)))),'b');
legend('Vout', 'Vin');
xlabel('Frequency');
ylabel('Voltage (dB)');
title(' Frequency vs voltage');

%vary Cn to see how the bandwidth changes C(6,6) = 1e-2

C(6,6) = 1e-2;
Simtime = 1;
numSteps = 1000;
deltaT = Simtime/numSteps;
% H vector
H = (C./deltaT) + G;

 V = zeros(6,1);
 Vp = V;

% F vector
F = zeros(1,6);
F(3)=In;
 %F(3) = In*randn();

Vv1 = zeros(10,1);
Vv2 = zeros(100,1);
newtime = linspace(1,numSteps,numSteps);

for i = 1:numSteps
    
    % F vector
    
    if (i == 40)
        F(1) = 1;
    end
    newtime(i) = i*deltaT;
    V = H\(((C * Vp)./deltaT) + F');

    Vv1(i) = V(1);
    Vv2(i) = V(5);
    
    Vp = V;
    
end

 figure(8)
  plot(newtime,Vv1,'r');
    hold on
   plot(newtime,Vv2,'b');
pause(0.01);
legend('Vin', 'Vout');
xlabel('Frequency');
ylabel('Voltage');
title('Time Simulation vs  The Vout noise source');
figure(9)
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv2)))),'g');
hold on 
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv1)))),'b');
legend('Vout', 'Vin');
xlabel('Frequency');
ylabel('Voltage (dB)');
title(' Frequency vs voltage');

%vary the time step 

C(6,6) = 1e-4;
% Simtime = 1;
% numSteps = 1000;
% deltaT = Simtime/numSteps;
deltaT=1e-5;
% H vector
H = (C./deltaT) + G;

 V = zeros(6,1);
 Vp = V;

% F vector
F = zeros(1,6);
F(3)=In;
 %F(3) = In*randn();

Vv1 = zeros(10,1);
Vv2 = zeros(100,1);
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

 figure(10)
  plot(newtime,Vv1,'r');
    hold on
   plot(newtime,Vv2,'b');
pause(0.01);
legend('Vin', 'Vout');
xlabel('Frequency');
ylabel('Voltage (dB)');
title('Time Simulation vs  The Vout noise source');
figure(11)
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv2)))),'g');
hold on 
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv1)))),'b');
legend('Vout', 'Vin');
xlabel('Frequency');
ylabel('Voltage (dB)');
title(' Frequency vs voltage');

%vary time
C(6,6) = 1e-4;
% Simtime = 10000;
% numSteps = 1000;
% deltaT = Simtime/numSteps;
deltaT=1e-8;
% H vector
H = (C./deltaT) + G;

 V = zeros(6,1);
 Vp = V;

% F vector
F = zeros(1,6);
F(3)=In;
 %F(3) = In*randn();

Vv1 = zeros(10,1);
Vv2 = zeros(100,1);
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

 figure(12)
  plot(newtime,Vv1,'r');
    hold on
   plot(newtime,Vv2,'b');
pause(0.01);
legend('Vin', 'Vout');
xlabel('Frequency');
ylabel('Voltage (dB)');
title('Time Simulation(1e-4) vs  The Vout noise source');
figure(13)
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv2)))),'g');
hold on 
plot(linspace(-500,500,1000),fftshift(20*log(abs(fft(Vv1)))),'b');
legend('Vout', 'Vin');
xlabel('Frequency');
ylabel('Voltage (dB)');
title(' Frequency vs voltage');

