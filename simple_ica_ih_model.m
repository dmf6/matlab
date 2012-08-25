close all;
clear all;

% capacitance in nF
c= 1;

% conductances in uS
gl = 0.05;
gca = 0.05;
gh = 0.03;

% reversal potentials in mV
el = -68.0;
eca = 120;
eh = -20;

vhalf_minf = -55;
k_minf = 4;

vhalf_mtau = -55;
k_mtau = 4;
tau1m = 140;
tau2m = -120;

vhalf_hinf = -60;
k_hinf = 4;

k_htau = 5;
vhalf_htau = -60;
tau1h = 300;
tau2h = -50;

vhalf_r = -55;
k_r = 4;
tau1r = 980;
tau2r = 150;
iapp = 0;

rinf = @(v) 1/(1+exp((v-vhalf_r)/k_r));
rtau = @(v) tau1r + tau2r/(1+exp(-(v-vhalf_r)/k_r));
minf = @(v) 1./(1+exp(-(v-vhalf_minf)/k_minf));
hinf = @(v) 1./(1+exp((v-vhalf_hinf)/k_hinf));
htau = @(v) (tau1h+tau2h/(1+exp(-(v-vhalf_htau)/k_htau)));

Tmax = 2000;
dt = 0.1;

%initial vectors 
t_vec = zeros(floor(Tmax/dt), 1);
v_vec = zeros(floor(Tmax/dt), 1);
r_vec = zeros(floor(Tmax/dt), 1);
h_vec = zeros(floor(Tmax/dt), 1);

v = -60;
r = 0;
m = minf(v); % set m = minf
h = 0.2;
v_vec(1) = v;
r_vec(1) = 0;
h_vec(1) = h;


%figure(1)
for i = 1:length(t_vec) - 1
    t_vec(i+1) = i * dt;
    
    k1_v = -gl*(v-el) - gh*r*(v-eh) - gca*minf(v)*h*(v-eca) + iapp;

    k1_v = k1_v / c;
    k1_r = (rinf(v)-r) / rtau(v);
    k1_h = (hinf(v)-h) / htau(v);
    a1_v = v + k1_v * dt; % y+k1
    a1_r = r + k1_r * dt;
    a1_h = h + k1_h * dt;
    
    k2_v = -gl*(a1_v-el) - gh*a1_r*(a1_v-eh) - gca*minf(a1_v)*a1_h*(a1_v-eca) + iapp;
    k2_v = k2_v / c;
    k2_r = (rinf(a1_v)-a1_r) / rtau(a1_v);
    k2_h = (hinf(a1_v)-a1_h) / htau(a1_v);
    
    v = v + (k1_v + k2_v) * dt / 2;
    r = r + (k1_r + k2_r) * dt / 2;
    h = h + (k1_h + k2_h) * dt / 2;
    
    v_vec(i+1) = v;
    r_vec(i+1) = r;
    h_vec(i+1) = h;
    
%     vv = -100:0.1:100;
% 	bnlc = hinf(vv);
%     subplot(1, 2, 1)
%     hold on
% 	plot(vv, bnlc, '-or', 'Linewidth', 1);
%     hold on
%     plot(v_vec(1:i), h_vec(1:i), '-ok', 'Linewidth', 3);
% 	
%     hold on
%     vnlc = (iapp-gl.*(vv-el)-gh.*r*(vv-eh))./(gca.*minf(vv).*(vv-eca));
% 	plot(vv, vnlc, '--ob', 'Linewidth', 2);
%     hold on
%     %hold off
% 	axis([-80 -20 -0.1 1.1])
%     xlabel('V (mV)');
%     ylabel('h');
%     
%     subplot(1, 2, 2)
%     plot(t_vec(1:i), v_vec(1:i), '-ok', 'Linewidth', 3);
%     axis([0 Tmax -70 -20]) ;
%     xlabel('Time (ms)');
%     ylabel('V (mV)');
% 	pause(0.001); drawnow;
end
% 
% figure(1)
% 
% subplot(1, 2, 1);
% plot(t_vec,v_vec,'b-', 'Linewidth', 4);
% axis([0 Tmax -70 -20]);
% 
% subplot(1, 2, 2);
% xlabel('V (mV)');
% ylabel('h');
% r = linspace(0, 1, 10);
% for i = 1:length(r) 
%     vv = -100:0.1:100;
%     hold on
%     bnlc = hinf(vv);
%     plot(vv, bnlc, '-r', 'Linewidth', 1);
%     vnlc = (iapp-gl.*(vv-el)-gh.*r(i).*(vv-eh))./(gca.*minf(vv).*(vv-eca));
%     plot(vv, vnlc, '--b', 'Linewidth', 1);
%     
% end
% %plot(v_vec, h_vec, 'ko', 'Linewidth', 1);
% subplot(1, 2, 2);
% h1 = plot(v_vec(1), h_vec(1),'k-', v_vec(2), h_vec(2), 'ko', 'Linewidth', 3);
% legend('Trajectory', 'V-nullcline', 'h-nullcline', 'Location', 'NE');
% pause(10e-15);
% axis([-80 -20 -0.1 1.1]) ; % specify a strange axis, that is not changed
% 
% for i=2:numel(t_vec) -1
%     % update plots
%     %subplot(1, 2, 2);
%     axis([-80 -20 -0.1 1.1]) ; % specify a strange axis, that is not changed
%     set(h1(1),'xdata',v_vec(1:i),'ydata',h_vec(1:i)) ;
%     set(h1(2),'xdata',v_vec(i),'ydata',h_vec(i)) ;
%     pause(10e-15); drawnow ; % visibility
% end
plot(t_vec, v_vec, 'ko', 'Linewidth', 1);

%linearize the system about the fixed point with these parameters 
f = @(v) -gl*(v-el) - gh*rinf(v)*(v-eh) - gca*(minf(v))*hinf(v)*(v-eca);                      
vinf = fzero(f, -65);

%assign state variables from y vector components
v = vinf;

% Find elements of Jacobian matrix of partial derivatives at fixed point
rprimeinf = (-(1/k_r)*exp((v-vhalf_r)/k_r))/((1+exp((v-vhalf_r)/k_r))^2);
mprimeinf = ((1/k_minf)*exp(-(v-vhalf_minf)/k_minf))/((1+exp(-(v-vhalf_minf)/k_minf))^2);
hprimeinf = (-(1/k_hinf)*exp((v-vhalf_hinf)/k_hinf))/((1+exp((v-vhalf_hinf)/k_hinf))^2);


%% Jacobian elements
Fv = -gl - gh*rprimeinf*(v) - gca*(mprimeinf*hinf(v) + hprimeinf*minf(v))*(v-eca) + minf(v)*hinf(v);
Fr = -gh*(v-eh);
Fh = -gca*minf(v)*(v-eca);
Gv = rprimeinf/rtau(v);
Gr = -1/rtau(v);
Gh = 0;
Hv = hprimeinf/htau(v);
Hr = 0;
Hh = -1/htau(v);
% Jacobian
jac = [Fv Fr Fh; Gv Gr Gh; Hv Hr Hh];

% use eigenvalues of jac to determine stability of linear system
eig(jac)

%% rescalings
% parameter 1
gL = Fv;

% parameter 2
g1 = gh*rprimeinf*(v-eh);

% parameter 3
g2 = gca*minf(v)*(v-eca)*hprimeinf;

amp = 5;
f = 10; 

tau1 = rtau(v);
     tau2 = htau(v);
     gammaL = (gL*tau1)/c;
     gamma1 = (g1*tau1)/c;  
     gamma2 = (g2*tau1)/c;
     sigma = tau1/tau2;
     T = 1000/tau1;
     Ain = (amp*tau1)/c;
v= 0;  
w1 =  0;
w2 = 0;
v_vec = zeros(floor(Tmax/dt), 1);
w1_vec = zeros(floor(Tmax/dt), 1);
w2_vec = zeros(floor(Tmax/dt), 1);
v_vec(1) = v;
w1_vec(1) = w1;
w2_vec(1) = w2;

for i = 1:length(t_vec) - 1
    t_vec(i+1) = (i * dt)/tau1;
    %%%%%%%%%%%%%%%%%%%%%%%%   RK2      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% k1 f(t, y) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    v = v_vec(i);
    t = t_vec(i);
   
   
    Iin = Ain*sin(2*pi*f*t/T); 

    % define ODE system
    k1_v = (Iin-gammaL*v - gamma1*w1 - gamma2*w2);
    k1_w1 = sigma*(v - w1);
    k1_w2 = (v - w2);
   
    a1_v = v + k1_v * dt; %y+k1
    a1_w1 = w1 + k1_w1 * dt;
    a1_w2 = w2 + k1_w2 * dt;
    
    %%%%%%%%%%%%%%%%%%%%%% k2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rescaling
    k2_v = (Iin-gammaL*a1_v - gamma1*a1_w1 - gamma2*a1_w2);
    k2_w1 = sigma*(k1_v - k1_w1);
    k2_w2 = (k1_v - k1_w2);
   
    v = v + (k1_v + k2_v) * dt / 2;
    w1 = w1 + (k1_w1 + k2_w1) * dt / 2;
    w2 = w2 + (k1_w2 + k2_w2) * dt / 2;
    
    v_vec(i+1) = v;
    w1_vec(i+1) = w1;
    w2_vec(i+1) = w2;
end

figure(1)
plot(t_vec, v_vec,'-or','LineWidth', 2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','b',...
            'MarkerSize', 5);


