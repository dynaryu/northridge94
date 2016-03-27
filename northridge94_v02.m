%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: H. RYU, J. KIM & J. BAKER (2009)                                %
%        "A PROBABILISTIC METHOD FOR THE MAGNITUDE ESTIMATION OF A        %
%         HISTORICAL DAMAGING EARTHQUAKE USING STRUCTURAL FRAGILITY       %
%         FUNCTIONS"                                                      %
%         Validation example: 1994 Northridge earthquake, pp.523-525      %
%                                                                         %
% by: Eduardo Charters Morais, PhD student, BME, Budapest, Hungary        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; tic

% Validation Example (M=6.7, Rep=20km)
Tn=0;    % Singhal & Kiremidjan (1998)-T=0.30s, though for PGA=Sa(T=0s)
Rep=20;  % average distance to Epicenter (in km)
F=0;     % Strike-Slip fault! (1-reverse; 0.5-reverse/oblique; 0-other)
HW=1;    % Hanging wall effect true-1, false-0
S=1;     % (Deep soil-S=1; rock-S=0;)
Fv=1.0;  % Soil Amplification Factor, HAZUS MH MR-4, p.5-58

% COEFFICIENT FOR STRUCTURAL FAILURE DATA
% Number of damages - Singhal & Kiremidjan (1988)
n=[174 4 6]; nt=sum(n);
% Combinatorial operation
syms n0 n1 n2 ntot
F1=subs(n0*factorial(n0-1),n(1));
F2=subs(n1*factorial(n1-1),n(2));
F3=subs(n2*factorial(n2-1),n(3));
Ft=subs(ntot*factorial(ntot-1),nt);
Coeff=Ft/(F1*F2*F3);

% PROBABILITY OF AN EVENT GIVEN A MAGNITUDE - PE|m
nit1=201;
nit2=301;
M=zeros(1,nit1);
PEm=zeros(1,nit1);
PGAar=zeros(1,nit1);
fx=zeros(1,nit1);
for i=1:nit1
    M(i)=5+(i-1)*3/nit1;
    [PGAmu,PGAs,PGAs0]=Abrahamson(Rep,M(i),Tn,F,HW,S);
    PGAar(i)=PGAmu;
    
    %from HAZUS technical manual, HAZUS MH MR-4;
    fact=interp1([1 5 6 7 8 10],[3.3 3.3 1.8 1.2 0.85 0.85],M(i))*1.5/Fv;   
    % Equivalent-PGA HAZUS MH MR-4, pages 5-64 to 5-67
    PGAds=[0.16 0.23]*fact;  % Pre-Code:[0.10 0.12] Low-Code:[0.12 0.15]
                             % Mod.-Code:[0.16 0.23] High-Code:[0.21 0.35]
    PGA=0:4/nit2:4;
    dpga=4/nit2;
    xPGA=PGA+dpga/2;
    A=normcdf(log(xPGA/PGAds(1))/0.64);
    B=normcdf(log(xPGA/PGAds(2))/0.64);
    C=normpdf(log(xPGA/PGAmu)/(PGAs));
    fx(i)=sum(((((1-A).^n(1)).*(A-B).^n(2)).*B.^n(3)).*C);
    PEm(i)=Coeff*fx(i)*(Rep*1000)*dpga;

end

% PROBABILITY OF A MAGNITUDE GIVEN AN EVENT - Pm|E
% Uniform distribution of magnitude
ml=5; mu=8; fm=1/(mu-ml);

% Limits of Distance, in km
Ri=5;
Rf=30;

% PLOTS
% PE|m
figure(1);
plot(M,PEm); grid on; hold on;
% Confirm Attenuation Relationship results
% figure(2);
% loglog(M,PGAas); grid on; hold on;
% Pm|E
% figure(3);
% plot(M,fm*ones(length(M))); hold on;
% plot(M,PmE); grid on;

t1=toc