%%GRUPPO S PROGETTO 1A
%%Stefano Favale, Gabriele Nunziati, Giorgia Quarta, Jacopo Jop
clear all, clc

%% Parameters
Id=0.5;     %[kg*m^2] Momento di Inerzia del drone rispetto all'asse di rotazione
beta=0.3;   %[] Coefficiente di attrito dinamico dovuto all'aria
aa=0.3;     %[m] Distanza dai punti di applicazione delle forze al baricentro
Fv=-8.7;    %[N] Forza del vento

%Fv=-7.7;    %[N] Forza del vento MAX
%Fv=-9.7     %[N] Forza del vento MIN


%%System dynamics
% dot{x_1}=x_2
% dot{x_2}mm=(-beta*x_2+aa/2*sin(x_1)*Fv+aa*u_1)/Id
% y=x_1

%Matrici del sistema linearizzato
A=[0,1; (aa*sqrt(2)*Fv)/(4*Id), -beta/Id];
B=[0;aa/Id];
C=[1,0];
D=0;

%Equilibrio

u_e= (-sqrt(2)*Fv)/4;
y_e= pi/4;       %=x_1

%% Equivalent transfer function
s=tf('s');
GG=(aa/Id)/(s^2+s*beta/Id-aa*sqrt(2)*Fv/(4*Id));
%GG=(4*aa)/(4*s^2*Id+4*s*beta-aa*sqrt(2)*Fv);

figure(1)
bode(GG)
figure(2)
step(GG)
figure(3)
margin(GG)

%% Requirements/specifications
% Gradino di ingresso--> w=W 1(t)
W=-pi/4;    % W dell'ingresso
% -> Sovraelongazione percentuale massima S%<=1%=0.01

%Sovraelongazione percentuale:
S_100= 0.01;
xi= sqrt(log(S_100)^2/(pi^2+log(S_100)^2)); %xi DESIDERATA (derivante dalla S-percententuale)

%per stare + tranquilli aumentiamo qui per un Mf aumentato 

Mf=xi*100; %(in good approximation) <--MARGINE DI FASE (DEFINIZIONE)
Mf=83;     %Per essere conservativo prendo l'intero successivo
display(Mf)

% -> Settling time (Tempo di assestamento)
%Ta,1<=0.3 [s]
T_a1Max=0.3;

%T_a1=4.6/(xi*omega_c) < 0.3
%omega_c > 4.6/(0.3*xi)
omega_cMin=4.6/(T_a1Max*xi);% <-- RICAVIAMO omega_c min dal tempo di assestamento massimo

% Measure Noise requirements 
% (Parametri dal testo)
omega_n=12000;
A_n=-29.542;    %<-- ABBATTIMENTO di 30 VOLTE

%% Graphical constraints of the step response
figure(2); 
T_simulation=20;                %secondi massimi rappresentati (da 0 a T_simulation)
%step(GG,T_simulation,'b')
step(GG*W,T_simulation,'b')     %<-- PLOTTO la risposta della G(s) col NOSTRO gradino, non quello unitario
hold on

% add overshoot constraint
patch([0,T_simulation,T_simulation,0],[W*(1+S_100),W*(1+S_100),W-1,W-1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
hold on; 
%ylim([0,W+1]);
ylim([W-1,0]);

%W - 1: Guardiamo il comportamento da 0 a -pi/4 -1
% add Settling time constraint
patch([T_a1Max,T_simulation,T_simulation,T_a1Max],[W*(1-0.01),W*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a1Max,T_simulation,T_simulation,T_a1Max],[W*(1+0.01),W*(1+0.01),W-1,W-1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

hold off
Legend=["G(s) Step Response";"Overshoot Constraint"; "Settling time Cons"];
legend(Legend);

%% Graphical constraints on Bode's plot
omega_plot_min=10^(-2);
omega_plot_max=10^6;

omega_cMax=omega_n;

[Mag,phase,omega]=bode(GG,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Nominal plant TF Bode diagram")

% % %% Adjusting the gain we cannot acheive the right performance
% % G_temp=GG/0.5*W;
% % 
% % [Mag,phase,omega]=bode(G_temp,{omega_plot_min,omega_plot_max},'k');
% % figure()
% % patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
% % patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
% % Legend_dB=["A_n";
% %     "\omega_{cMin}"
% %     "G(s)" ];
% % legend(Legend_dB);
% % hold on;
% % margin(Mag,phase,omega);grid on;
% % patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on
% % Legend_arg=["G(s)";"M_f"] ;
% % legend(Legend_arg)
% % hold off
% % step(G_temp)

%% Inserting the internal model

%Viene inserito un polo a w=0 per assicurarci un errore a regime nullo
%La nostra R_statica(s) sarà quindi 1/s
%La G_estesa è definita come il prodotto tra la R_statica(s) e la G(s)
GG_e=GG/s;
[Mag,phase,omega]=bode(GG_e,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G_e(j\omega)" ];

legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G_e(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Extended plant TF Bode diagram")

%% Prima parte della rete dinamica
%recupero di fase maggiore di 90°, devo avere due zeri o due reti anticipatrici
%Aggiungo uno zero prima di omega_c_min

%we can add for free a zero
tau_z=1/(omega_cMin*0.1); %= 1/omega_z

display(omega_cMin);
display(1/tau_z);

R_1=1+tau_z*s; %GG_e*R_1=G(s)* R_1/s; thus respecting the causality property

[Mag,phase,omega]=bode(GG_e*R_1,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G_e R_1(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G_e R_1(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Plant with PI TF Bode diagram")

%% Rete anticipatrice
omega_cStar=(150); 

%R_ant=(1+tau*s)/(1+atau*s), a<1
%tau= (M_star - cos(phi_star))/(omega_cStar * sin(phi_star))

[Mag_omega_cStar,phase_omega_cStar,omega_cStar]=bode(GG_e*R_1,omega_cStar);

phi_star= (Mf+5) - 180 -(phase_omega_cStar);
display(phi_star);

%phi star tra 0 e 90 gradi
% so to fulfill cos(phi_star)>1/M_star
M_star = 50;
display(1/M_star)
display(cosd(phi_star))

tau= (M_star - cosd(phi_star))/(omega_cStar * sind(phi_star));
atau=(cosd(phi_star) - 1/M_star)/(omega_cStar * sind(phi_star));

R_ant=(1+tau*s)/(1+atau*s);
figure()
bode(R_ant)
title("Rete anticipatrice")

%% 
[Mag,phase,omega]=bode(GG_e*R_1*R_ant,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G_e R_1 R_{ant}(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G_e R_1 R_{ant}(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Plant with PI+ R_{ant} TF Bode diagram")
zoom on

figure()
LL=GG*R_1/s*R_ant;
step(LL/(1+LL))

%% Define the gain to adjust the cross frequency
gain=1/abs(evalfr(GG_e*R_1*R_ant,1i*omega_cStar));

CC=2.6749;
LL=GG_e*R_1*R_ant*gain*CC;

% 2.6749 c dal control system design


[Mag,phase,omega]=bode(LL,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "L(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["L(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Open loop TF Bode diagram")
zoom on

figure()
step(W*LL/(1+LL),5)

% add Settling time constraint
patch([0,T_simulation,T_simulation,0],[W*(1+S_100),W*(1+S_100),W-1,W-1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
hold on; 

patch([T_a1Max,T_simulation,T_simulation,T_a1Max],[W*(1-0.01),W*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a1Max,T_simulation,T_simulation,T_a1Max],[W*(1+0.01),W*(1+0.01),W-1,W-1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

hold off
Legend=["G(s) Step Response";"Overshoot Constraint"; "Settling time Cons"];
legend(Legend);
grid on
zoom on % ylim([0.94 1.07]);xlim([0 2])

%% Plot graphic with GG, L 
%Confrontiamo il grafico iniziale della G(s) con quello finale 
%comprendente il nostro regolatore
%Attenzione: la legenda non è corretta

gain=1/abs(evalfr(GG_e*R_1*R_ant,1i*omega_cStar));

L = tf(LL,{omega_plot_min,omega_plot_max},'k');
G = tf(GG,{omega_plot_min,omega_plot_max},'b');



figure();
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["";"M_f"] ;
legend(Legend_arg)


bode(L,G);
hold off
title("Open loop TF Bode diagram")
zoom on

%% Simulink addictional part

% Define the controller coefficients to run the simulink file
%R=gain*gain_CSD*R_ant*R_1/s; %R cioè Rs * Rd * gain * CC

R=(1/s)*R_1*R_ant*gain*CC;
display(R);

[n_r,d_r]=tfdata(R);
num_r=n_r{1};
den_r=d_r{1};

%[n_g,d_g] = tfdata(GG);
%num_g=n_g{1};
%den_g=d_g{1};

%initial condition

x0=[pi/4;0];

%Caso LINEARE:
%x0_lin=[0;0];

%Caso NON LINEARE:
%x0_non_lin=[pi/4;0];

%Scegliere x0 in base al caso da simulare


open("simulink_1a_S_definitivo");

