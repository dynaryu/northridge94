function [Smu,Ss,PGAs0,f1] = Abrahamson(Rrup,M,Tn,F,HW,S)

% Interpolation of Abrahamson Coefficients
[c4,a1,a2,a3,a4,a5,a6,a9,a10,a11,a12,a13,c1,c5,n1,b5,b6]=interpolation1(Tn);
R=sqrt(Rrup^2+c4^2);

% Calculation of f1
if M<=c1
    f1=a1+a2*(M-c1)+a12*(8.5-M)^n1+(a3+a13*(M-c1))*log(R);
elseif M>c1
    f1=a1+a4*(M-c1)+a12*(8.5-M)^n1+(a3+a13*(M-c1))*log(R);
end

% Calculation of f3 (Fault-type factor)
if M<=5.80
    f3=a5;
elseif (M>5.8)&&(M<c1)
    f3=a5+(a6-a5)/(c1-5.8);
elseif M>=c1
    f3=a6;
end

% Calculation of fHWm
if M<=5.50
    fHWm=0;
elseif (M>5.5)&&(M<6.5)
    fHWm=M-5.5;
elseif M>=6.5
    fHWm=1;
end

% Calculation of fHWr
if Rrup<=4
    fHWr=0;
elseif (Rrup>4)&&(Rrup<=8)
    fHWr=a9*(Rrup-4)/4;
elseif (Rrup>8)&&(Rrup<=18)
    fHWr=a9;
elseif (Rrup>18)&&(Rrup<=25)
    fHWr=a9*(1-(Rrup-18)/7);
elseif Rrup>=25
    fHWr=0;
end

% Calculation of Standard Deviation
if M<=5.0
    Ss=b5;
elseif (M>5.0)&&(M<7.0)
    Ss=b5-b6*(M-5);
elseif M>=7
    Ss=b5-2*b6;
end

% Calculation of f4 (Hanging wall effect)
f4=fHWm*fHWr;

% Calculation of PGA
[PGAs0,PGAstd]=PGAcalc(Rrup,R,M,F,HW);

% Calculation of f5
f5=a10+a11*log(PGAs0+c5);   % Site response

Smu=exp(f1+F*f3+HW*f4+S*f5);

end

