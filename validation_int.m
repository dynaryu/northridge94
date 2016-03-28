%function validation(GMpdf)

% Validation of the proposed method 
Mv = 5.0:0.05:8.0; % and use G-R attenuation as a prior fM dist.
Rv = 20.0; % scalar value for Rv
Epsilon = -3.5:0.01:3.5;

%if (nargin < 1) % default fsamp
%    disp('Uniform dist of Magnitude');
%end
GMpdf(1) = 1; GMpdf(2) = Mv(1); GMpdf(3) = Mv(end);
%GMpdf(1) = 3; GMpdf(2) = 0.8; 

% Given event 
% ns = 1; nd = 2 (0 none, 1 slight, 2 moderate)
% n0 = 184-17+7 = 174; ; in the none damage state 
% n1 = 4 ; in the slight damage state
% n2 = 6 ; in the moderate damage state
% nt = 184 ; the total number of structures 
nBldg = [ 174 4 6 0 0 ]; % none(1) slight(2) moderate(3) extensive(4) complete(5) 
nt = sum(nBldg);

% construct fragility function
% median PGA|ds for C1L High seismic code level
PGA_Rds = [0.21 0.35 0.70 1.37]; % slight(1) moderate(2) extensive(3) complete(4) 
beta_SPGA = 0.64; 

% Spectrum shape ratio R_{PGA/SA1} - WUS ROCK (Site Class B) 
% From HAZUS Technical manual
rpga = [ 5 6 7 8 ; 3.3 1.8 1.2 0.85 ]'; % r = 20.0 km 
rpgaM = repmat( [ 5 6 7 8 ], [4,1]); % Magnitude
rpgaR = repmat( [ 10 20 40 80 ]', [1,4]); % Distance
rpgaMR = [ 3.8 2.1 1.5 0.85 ; 3.3 1.8 1.2 0.85 ; 2.9 1.6 1.05 0.80 ; 3.2 1.7 1.0 0.75];

% PGA(ds) = PGA(R,ds) * R(PGA/SA1)*(1.5/Fv)
% R(PGA/SA1) = fn_spec_shape_ratio(M,R)
jR=1;
Fv = 1.0; % Assume that all building are located in the site class B(Rock site)
PGA_ds = zeros(length(PGA_Rds),length(Mv),length(Rv));
% slight(1) moderate(2) extensive(3) complete(4) 
for iM=1:length(Mv);

%    rs_ratio(iM,jR) = fn_spec_shape_ratio(Mv(iM),Rv(jR)); % compute PGA/SA1 using A&S attenuation relationship
%     rs_ratio(iM,jR) = INTERP1(rpga(:,1),rpga(:,2),Mv(iM),'linear'); % table lookup for R=20km
      rs_ratio(iM,jR) = interp2(rpgaM,rpgaR,rpgaMR,Mv(iM),Rv(jR),'linear'); % table lookup for M,R

    for ids=1:4;
        PGA_ds(ids,iM,jR) = PGA_Rds(ids).*rs_ratio(iM,jR).*1.5./Fv;
    end
end

% PDF of epsilon
Fepsilon = normpdf(Epsilon); % [1 x nEpsilon]

% initialize
PE_MR = zeros(length(Mv),length(Rv)); % 

% P(E|M,R) = int[P(E|EPS,M,R)*f(EPS|M,R)]dEPS]
for iM=1:length(Mv); %(iM
    [pga sigma_atten] = abrahamson_atten(Mv(iM),Rv(1),0.01); % median PGA, std of lnPGA given M,R   
    lnPGA = log(pga) + Epsilon.*sigma_atten; 
    
% compute P(DS|PGA,M,R)
% none(1)  
    p_ds(:,1) = 1-normcdf((lnPGA-log(PGA_ds(1,iM,jR)))./beta_SPGA);

% slight(2) moderate(3) extensive(4) 
    for ids=2:4 
        p_ds(:,ids) = normcdf((lnPGA-log(PGA_ds(ids-1,iM,jR)))./beta_SPGA) ... 
                     -normcdf((lnPGA-log(PGA_ds(ids,iM,jR)))./beta_SPGA);
    end
%complete(5)    
    p_ds(:,5) = normcdf((lnPGA-log(PGA_ds(4,iM,jR)))./beta_SPGA) ;

    lninteps = log(p_ds)*nBldg'; % [nEps*nds]x[nds*1] 
    inteps = prod(175:184)./(prod(4:1)*prod(6:1)).*exp(lninteps).*Fepsilon';
%    log_inteps = n0.*log(p_none) + n1.*log(p_slight) + n2.*log(p_moderate);
%    inteps = exp(log_inteps);

    PE_MR(iM,jR) = trapz(Epsilon,inteps);  % [nMv]
end %iM)

%%%%%%%%
% Mv
%%%%%%%%
opt_M = GMpdf(1);
if opt_M == 1 ; %  option 1: unifrom dist. ML, MU 
   Mpdf = unifpdf(Mv,Mv(1),Mv(end));

elseif opt_M == 2; % option 2: pdf of M based on previous studies
   mean_M = GMpdf(2); %5.8;
   sig_M = GMpdf(3); %0.5;
   [ mulnM, siglnM ] = norm2lognorm(mean_M,sig_M);
   Mpdf = lognpdf(Mv,mulnM,siglnM);

elseif opt_M == 3; % option 2: pdf of M based on G-R law
%   bvalue = 0.6; %, 0.8, 0.9
   bvalue = GMpdf(2);
   beta_value = bvalue.*log(10);
   Mpdf = beta_value*exp(-beta_value.*Mv)./(exp(-beta_value.*Mv(1))-exp(-beta_value.*Mv(end)));
end
pE_M = PE_MR';

% P(E) = int P(E|M)*f(M)
pE = trapz(Mv,pE_M.*Mpdf); 
if pE == 0 
   fM_E = 0;
else
   fM_E = pE_M.*Mpdf./pE;
end
% calc expv, stdv
expv = trapz(Mv,Mv.*fM_E); stdv = sqrt(trapz(Mv,Mv.*Mv.*fM_E) - expv.*expv);   

%subplot(3,1,1); plot(Mv,pE_M,'linewidth',1.5); xlabel('Magnitude'); ylabel('P(E|m)'); grid on;  
%subplot(3,1,2); plot(Mv,Mpdf,'linewidth',1.5); xlabel('Magnitude'); ylabel('f(m)'); grid on;
%subplot(3,1,3); plot(Mv,fM_E,'linewidth',1.5); xlabel('Magnitude'); ylabel('f(m|E)'); grid on; title(strcat('E(m|E)=',num2str(expv,2),', \sigma(m|E)=',num2str(stdv,2))); 

% P(E|m)
figure;
plot(Mv,pE_M,'k-',Mv(1)*ones(1,5),linspace(0.0,pE_M(1),5),'k--',Mv(end).*ones(1,5),linspace(pE_M(end),0.0,5),'k--','linewidth',2.0); 
xlabel('Magnitude','FontName','Times New Roman','FontSize',14); 
ylabel('{\it P(E|m)}','FontName','Times New Roman','FontSize',14); 
set(gca,'FontName','Times New Roman','FontSize',14);
xlim([4.9 8.1]); ylim([0.0 0.6]);
set(gca,'XTickLabel',{'5.0' '5.5' '6.0' '6.5' '7.0' '7.5' '8.0'})
set(gca,'YTick',[0.0:0.1:0.6],'YTickLabel',{'0' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6'})
grid on; 


% f(M|E)
%Mpdf1v = [ linspace(0.0,max(Mpdf),5) Mpdf linspace(max(Mpdf),0.0,5) ];
%Mv1 =    [ Mv(1)*ones(1,5)          Mv     Mv(end).*ones(1,5) ];
figure;
plot(Mv,Mpdf,'k-+',Mv,fM_E,'k-',Mv(1)*ones(1,5),linspace(0.0,Mpdf(1),5),'k--',Mv(end).*ones(1,5),linspace(Mpdf(end),0.0,5),'k--',Mv(1)*ones(1,5),linspace(0.0,fM_E(1),5),'k--',Mv(end).*ones(1,5),linspace(fM_E(end),0.0,5),'k--','linewidth',2.0); 
%plot(Mv1,Mpdf1v,'b-',,'linewidth',2.0); 
xlabel('Magnitude','FontName','Times New Roman','FontSize',14); 
ylabel('Probability density','FontName','Times New Roman','FontSize',14); 
set(gca,'FontName','Times New Roman','FontSize',14); 
legend('Prior distribution: {\it f_{M}(m)}','Posterior distribution: {\it f_{M|E}(m|E)}',2)
xlim([4.9 8.1]); ylim([0.0 1.0]); grid on;
set(gca,'XTickLabel',{'5.0' '5.5' '6.0' '6.5' '7.0' '7.5' '8.0'})
set(gca,'YTick',[0.0:0.2:1.0],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'})

