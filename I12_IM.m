function [I1, I2] = I12_IM(j, i, DZ, DT, FCm, FSm, FCp, FSp, Cm, Cp, Sm, Sp)
% function [I1, I2] = integrate(j, i, DZ, DT)

j = j+1;

SQRT2PIT = sqrt(2*pi*DT);
DT2 = 2*DT;

Z = i*DZ;
i = i + 1;

% Jmj = (j*DZ-i*DZ)/SQRT2PIT;
% Jmjm1 = ((j-1)*DZ-i*DZ)/SQRT2PIT;
% Jpj = (j*DZ+i*DZ)/SQRT2PIT;
% Jpjm1 = ((j-1)*DZ+i*DZ)/SQRT2PIT;
% 
% Kmj = (j*DZ-i*DZ).^2/(4*DT);
% Kmjm1 = ((j-1)*DZ-i*DZ).^2/(4*DT);
% Kpj = (j*DZ+i*DZ).^2/(4*DT);
% Kpjm1 = ((j-1)*DZ+i*DZ).^2/(4*DT);

% I1 = SQRT2PIT*(fresnelc(Jmj) - fresnelc(Jmjm1) - fresnelc(Jpj) + fresnelc(Jpjm1) + ...
%      1i*(-fresnels(Jmj) + fresnels(Jmjm1) + fresnels(Jpj) - fresnels(Jpjm1)));
%  
% I2 = SQRT2PIT*Z*fresnelc(Jmj) - SQRT2PIT*Z*fresnelc(Jmjm1) + SQRT2PIT*Z*fresnelc(Jpj) - ...
%      SQRT2PIT*Z*fresnelc(Jpjm1) + 1i*(DT2*(cos(Kmj) - cos(Kmjm1)) + SQRT2PIT*Z*(-fresnels(Jmj) + fresnels(Jmjm1))) + ...
%      1i*(DT2*(-cos(Kpj) + cos(Kpjm1)) + SQRT2PIT*Z*(-fresnels(Jpj) + fresnels(Jpjm1))) + DT2*(sin(Kmj) - sin(Kmjm1)) + ...
%      DT2*(-sin(Kpj) + sin(Kpjm1));

I1 = SQRT2PIT*(FCm(j,i) - FCm(j-1,i) - FCp(j,i) + FCp(j-1,i) + ...
     1i*(-FSm(j,i) + FSm(j-1,i) + FSp(j,i) - FSp(j-1,i)));
 
I2 = SQRT2PIT*Z*FCm(j,i) - SQRT2PIT*Z*FCm(j-1,i) + SQRT2PIT*Z*FCp(j,i) - ...
     SQRT2PIT*Z*FCp(j-1,i) + 1i*(DT2*(Cm(j,i) - Cm(j-1,i)) + SQRT2PIT*Z*(-FSm(j,i) + FSm(j-1,i))) + ...
     1i*(DT2*(-Cp(j,i) + Cp(j-1,i)) + SQRT2PIT*Z*(-FSp(j,i) + FSp(j-1,i))) + DT2*(Sm(j,i) - Sm(j-1,i)) + ...
     DT2*(-Sp(j,i) + Sp(j-1,i));
end

