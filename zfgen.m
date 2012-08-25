clear all;
close all;
out = importdata('2008_10_20_20_38_48.asc');
%load 'out';
fmin=0.1; fmax=4;
%iscale=.001; %(nA/pA)
iscale=1.;
t=out(:,1);
ii=iscale*out(:,2); % To convert to nA.
vv=out(:,3); % already in mV
v=vv - mean(vv);
i=ii - mean(ii);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = t(2)-t(1);
L = length(i);
T = L * dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,1);
plot(t,ii);
xlim([0,0.99*T]);
axis([0,T,-2, 2]);
xlabel('Time (s)');
ylabel('I (nA)');
grid on 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,2);
plot(t,v);
%axis([0,0.99*T,-10, 10]);
xlabel('Time (s)');
ylabel('V (mV)');
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% difficult to identify frequency components by looking at the original
% signal. We thus convert to the frequency domain
Nfft = 2^nextpow2(L); % returns the smallest power of two that is greater than or equal to the absolute value of L 
                      % Next power of 2 from length of v. For instance if
                      % we are calculating the FFT over 512 then Nfft will
                      % be 2^9 points. 9 is the lowest power that is
                      % greater than 512. Nfft is thus the number of
                      % samples over the sample period
Fs=L/(T/1000); % number of timesteps divided by the duration (in ms) gives 
               %us sampling rate e.g. 3 sampling intervals in a 30ms sweep would give 0.1Hz. Sample rate 
               % determines the highest frequency that can be represented
               % in the frequency domain. Number of points and rate
               % determine the range and precision of the DFT of the signal
f = Fs/2*linspace(0,1,Nfft/2); % one half of the function represents all the information. The frequency axis range covers half the spectrum = (sample rate)/2
V = fft(v,Nfft)/L; %do a Nfft point fft...we divide by the number of points to normalize the frequencies
Vs = abs(V(1:Nfft/2)); %because of symmetry around mid-point, it is common practice to show only the values corresponding to positive frequencies
%Vssm=smooth(Vs,0.04,'rloess');
I = fft(i,Nfft)/L;
Is = abs(I(1:Nfft/2));
subplot(3,2,3); 
%figure(2)


plot(f,Is,'.');
xlabel('Freq. [Hz]');
ylabel('|I(f)|');
grid on 
%axis([0,fmax,0,+Inf]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = V./I; %impedance calculation 
Zs=abs(Z(1:Nfft/2));
subplot(3,2,4)
%figure(2)
fcut=f(1:1:1049); %adjust to fit all frequencies that capture the bell curve
Zscut=Zs(1:1:1049);
%axis([fmin,fmax,0,20]);
%xlim([fmin,fmax]);
hold on

grid on
Zscutsm=smooth(Zscut,0.3,'rloess'); %A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
plot(fcut, Zscut, 'ob');
plot(fcut, Zscutsm, '-r','LineWidth', 2)
xlabel('Freq. [Hz]');
ylabel('|Z(f)|');
axis([fmin, fmax, -inf, +inf]);
%legend('Original Data','Smoothed Data Using ''rloess''',...
%'Location','NE')
phase = angle(Z(1:1:1049));
phasesm=smooth(phase,0.04,'rloess');
subplot(3,2,5:6);
plot(fcut,phase, 'or','LineWidth', 2,...
           'MarkerEdgeColor','k',...
           'MarkerFaceColor','b',...
           'MarkerSize', 5);
%plot(fcut,abs(Phasesm))
xlabel('Frequency (Hz)')
ylabel('Phase (Radians)')
axis([0.1, fmax, -inf, +inf]);
grid on