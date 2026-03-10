%%
% Nume si prenume:  Deac Lorena-Daniela
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 3; 
n = 13; 

%% Process data (fixed, do not modify)
c1 = (1000+n*300)/10000;
c2 = (1.15+2*(m+n/10)/20);
a1 = 2*c2*c1;
a2 = c1;
b0 = (1.2+m+n)/5.5;

rng(m+10*n)
x0_slx = [2*(m/2+rand(1)*m/5); m*(n/20+rand(1)*n/100)];

%% Experiment setup (fixed, do not modify)
Ts = 10*c2/c1/1e4*1.5; % fundamental step size
Tfin = 30*c2/c1*10; % simulation duration

gain = 10;
umin = 0; umax = gain; % input saturation
ymin = 0; ymax = b0*gain/1.5; % output saturation

whtn_pow_in = 1e-6*5*(((m-1)*8+n/2)/5)/2*6/8; % input white noise power and sampling time
whtn_Ts_in = Ts*3;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(10); % input quantizer (DAC)

whtn_pow_out = 1e-5*5*(((m-1)*25+n/2)/5)*6/80*(0.5+0.3*(m-2)); % output white noise power and sampling time
whtn_Ts_out = Ts*5;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(9); % output quantizer (ADC)

u_op_region = (m/2+n/5)/2; % operating point

%% Input setup (can be changed/replaced/deleted)
wf= 1/6.54;% de la proiectul trecut 1/T1 ,cu T1 usor de citit aproximativ din raspunsul la treapta (y63 de exemplu)
fmin=wf/2/pi/10;
fmax=wf/2/pi*10;
Ain=1.78; %arbitrar peste niovelul de zgomot 

%% Data acquisition (use t, u, y to perform system identification)
out = sim("DEAC_Lorena_circuit_hidraulic_chirp_R2025a.slx");

t = out.tout;
u = out.u;
y = out.y;

plot(t,u,t,y)
shg

%% System identification
%%yst=(7.24811+5.49716)/2
%%ust=(3.82812+0.234375)/2

%%yst=(6.84091+5.9858)/2
yst=(8.14394+4.68277)/2
ust=u_op_region
k=yst/ust

w1=pi/(310.489-304.056)%585.288-583.117  447.659-445.197  818.593-816.726 576.195-572.819
Delta_T1= 307.186-304.056 %818.226-816.726  443.238-440.627 581.782-579.654 575.093-572.819  575.093-572.819

phi1=rad2deg(-w1*Delta_T1)

w2=pi/(322.736-316.593)     %814.864-812.978  436.186-431.63 576.326-572.834 814.864-812.978 583.059-579.721
Delta_T2=319.708-316.593    %814.351-812.978  438.716-434.237 578.42-574.991

phi2=rad2deg(-w2*Delta_T2)%ar trebui sa ajung pe undeva pe la -88,-90,-91 grade ca sa fie bine 
%%
Ay1=(8.14394-4.68277)/2   %6.80019-6.02652  7.61458-5.1714
Au1=(3.84766-0.234375)/2
Ay2=(8.0625-4.7642)/2  %3.80859-0.3125  3.81836-0.263672
Au2=(3.79883-0.253906)/2
M=Ay1/Au1
Im=-M;%%M
%%
wn=w1 %ala care imi da defazaj bun 
zeta=-k/2/Im  %zeta ar trebui sa dea peste 1 din cauza polilor reali

H=tf(k*wn^2,[1,2*zeta*wn,wn^2])
zpk(H)

%%
%spatiul starilor +constante de timp

T1=1/0.167 % iau valoarea mai mica din functia de transfer pentru ca T1 e cel mare 
T2=1/1.428 %iau valoarea mai mica din functia de transfer pentru ca T2 e cel mai mic 

A=[0,1;
-1/T1/T2,-(1/T1+1/T2)];
B=[0;k/T1/T2];
C=[1,0];
D=0;

sys=ss(A,B,C,D);
ysim2=lsim(sys,u,t,[y(1);6]);


figure
plot(t,u,t,y,t,ysim2)

%Calculul erorilor

J=1/sqrt(length(t))*norm(y-ysim2);
eMPN=norm(y-ysim2)/norm(y-mean(y))*100;
%%
%functie de transfer
H=tf(sys);
zpk(H);

ysim1=lsim(H,u,t);
figure;
plot(t,y,t,ysim1)
%%
%spatiul starilor cu tita si wn

A=[0,1;
-wn^2,-2*zeta*wn];
B=[0;k*wn^2];
C=[1,0];
D=0;

sys=ss(A,B,C,D);
ysim=lsim(sys,u,t,[y(1);6]);


figure
plot(t,u,t,y,t,ysim)

%Calculul erorilor

J1=1/sqrt(length(t))*norm(y-ysim);
eMPN1=norm(y-ysim)/norm(y-mean(y))*100;

