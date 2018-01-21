%% Green's function for 1-D chain
clear all; clc;
%% Number of atoms in each section
NumD = 3;
%% Constants - mes - Unit:atomic m unit
ml = 20.0;  md = 20.0;  mr = 20.0; % in AMU(atomic m unit)
% Lattice Constant
a = 5e-10;  % lattice constant: unit [m]
% Unit Conversion of the m: times 1.66e-27kg = 1.66e-24g
amu = 1.66e-27;
ml = ml*amu; mr = mr*amu; md = md*amu;
%% Constants - k Constants
kl = 50.0;  kd = 40.0; kr = 60.0; % in Unit: N/m
% simple mixing rule
kld = sqrt(kl*kd);  krd = sqrt(kr*kd);
%% Define the calculation range
OMax = 8e13;   % Max frequency
OMin = 0e13;      % Min frequency
ONum = 100;    % Interval Number
OInterval = (OMax - OMin)/ONum; % Frequency interval
freq = zeros(ONum,1);
T = zeros(ONum,1);
DOS = zeros(ONum,1);
%% Undecided Vectors/Matrices
% The interaction:t1;t2; dimensions are DxL and DxR;
%% Delta and frequency.
for j = 1:ONum
freq(j) = OMin + j*OInterval;
Omega = freq(j);
delta = Omega^2*1e-20;

%% Construct Harmonic Matrix and Interaction Matrix
% Assume the atoms is from Left to right arrange.
% Conjugate Transpose: conj(A');

% Device Harmonic Matrix: Checked.
HD = zeros(NumD, NumD);
for i = 1:NumD
    if i-1>=1
        HD(i,i-1) = -kd/md;
    end
    if i+1<=NumD
        HD(i,i+1) = -kd/md;
    end
    HD(i,i) = 2*kd/md;
end

HD(1,1) = (kd+kld)/md;
HD(NumD,NumD) = (kd+krd)/md;
%% END OF HARMONIC MATRICES

%% Analytical Calculation for the bulk surface Green's function: 1x1 matrix: seem to work
%  Left
t = kl./ml;
v = Omega.^2 + 1i*delta;
tmp1 = (v-2*t-sqrt((2*t-v).^2-4*t*t))/(2*t*t);
tmp2 = (v-2*t+sqrt((2*t-v).^2-4*t*t))/(2*t*t);
if imag(tmp1)>0
    glcbs = tmp1;
else
    glcbs = tmp2;
end
%  Right
t = kr./mr;
tmp1 = (v-2*t-sqrt((2*t-v).^2-4*t*t))/(2*t*t);
tmp2 = (v-2*t+sqrt((2*t-v).^2-4*t*t))/(2*t*t);
if imag(tmp1)>0
    grcbs = tmp1;
else
    grcbs = tmp2;
end
%% End of surface green's function

%% The surface green's function
g1s = 1/(v-(kld+kl)/ml-(kl/ml).^2*glcbs);
g2s = 1/(v-(krd+kr)/mr-(kr/mr).^2*grcbs);
%% Submatrices of the interaction Matrix
t1 = zeros(NumD,1);                % N x 1 Matrix
t1(1,1) = -kld/sqrt(ml*md);        % the only nonzero term
t2 = zeros(NumD,1);                % N x 1 Matrix
t2(NumD,1) = -krd/sqrt(mr*md);     % the only nonzero term
% t1 = -kld/sqrt(ml*md);
% t2 = -krd/sqrt(mr*md);
%% Calculate the self-energies and Gamma
Sigma1 = t1*g1s*t1';                    % (N x 1) x (1 x 1) x (1 x N) = N x N Matrix
Sigma2 = t2*g2s*t2';                    % (N x 1) x (1 x 1) x (1 x N) = N x N Matrix
% Sigma1 = zeros(NumD);
% Sigma2 = zeros(NumD);
% Sigma1(1,1) = t1*g1s*t1;
% Sigma2(NumD,NumD) = t2*g2s*t2;
Gamma1 = 1i*(Sigma1 - Sigma1');         % N x N Matrix
Gamma2 = 1i*(Sigma2 - Sigma2');         % N x N Matrix
%% End

%% Green's function in the device
Gd = inv(Omega.^2*eye(NumD)-HD-Sigma1-Sigma2);
%% Transmission
T(j) = real(trace(Gamma1*Gd*Gamma2*Gd'));
A = 1i*(Gd-Gd');
DOS(j) = abs(trace(A)*Omega/pi/a);
end
% Plot

fp = fopen('./AGF_analytical.dat','a+');
fprintf(fp, 'Variables = "frequency","Transmission"\n');
fprintf(fp,'Zone T = "My code with Analytic gs",I = %d, DataPacking = Point\n',ONum);
for i = 1:ONum
    fprintf(fp,'%e   %e\n',freq(i),T(i));
end
fclose(fp);
figure(1)
plot(freq,T,'b-');
title('Transmission Function for 1D atom chain');
xlabel('frequency,rad/s');
ylabel('Transmission function');
figure(2)
plot(freq,DOS,'-');
title('DOS for the device');
xlabel('frequency,rad/s');
ylabel('DOS');