function [I1, I2] = I12_IM(zj, Z, DZ, DT)

ZjmZ = (zj - Z)/sqrt(DT*2*pi);
ZjmDZmZ = (zj - DZ - Z)/sqrt(DT*2*pi);
ZjpZ = (zj + Z)/sqrt(DT*2*pi);
ZjmDZpZ = (zj - DZ + Z)/sqrt(DT*2*pi); 
SqrZjmZd4DT = (zj - Z).*(zj - Z)/(4*DT);
SqrZjmDZmZd4DT  = (zj - DZ - Z).*(zj - DZ - Z)/(4*DT);
SqrZjpZd4DT = (zj + Z).*(zj + Z)/(4*DT);
SqrZjmDZpZd4DT = (zj - DZ + Z).*(zj - DZ + Z)/(4*DT);

Coeff = 1/2*sqrt(1i/(pi*DT));

I1 = DT.^(1/2).*(2.*pi).^(1/2).*((-1).*fresnelc(ZjmDZmZ)+fresnelc( ...
  ZjmDZpZ)+fresnelc(ZjmZ)+(-1).*fresnelc(ZjpZ)+(sqrt(-1)*(-1)).*(( ...
  -1).*fresnels(ZjmDZmZ)+fresnels(ZjmDZpZ)+fresnels(ZjmZ)+(-1).* ...
  fresnels(ZjpZ)));

I2 = (-1).*DT.^(1/2).*(2.*pi).^(1/2).*Z.*fresnelc(ZjmDZmZ)+(-1).*DT.^( ...
  1/2).*(2.*pi).^(1/2).*Z.*fresnelc(ZjmDZpZ)+DT.^(1/2).*(2.*pi).^( ...
  1/2).*Z.*fresnelc(ZjmZ)+DT.^(1/2).*(2.*pi).^(1/2).*Z.*fresnelc( ...
  ZjpZ)+sqrt(-1).*(2.*DT.*((-1).*cos(SqrZjmDZmZd4DT)+cos( ...
  SqrZjmZd4DT))+DT.^(1/2).*(2.*pi).^(1/2).*Z.*(fresnels(ZjmDZmZ)+( ...
  -1).*fresnels(ZjmZ)))+sqrt(-1).*(2.*DT.*(cos(SqrZjmDZpZd4DT)+(-1) ...
  .*cos(SqrZjpZd4DT))+DT.^(1/2).*(2.*pi).^(1/2).*Z.*(fresnels( ...
  ZjmDZpZ)+(-1).*fresnels(ZjpZ)))+2.*DT.*((-1).*sin(SqrZjmDZmZd4DT)+ ...
  sin(SqrZjmZd4DT))+2.*DT.*(sin(SqrZjmDZpZd4DT)+(-1).*sin( ...
  SqrZjpZd4DT));

I1 = I1*Coeff;
I2 = I2*Coeff;

% j = j+1;
% 
% SQRT2PIT = sqrt(2*pi*DT);
% DT2 = 2*DT;
% 
% Z = i*DZ;
% i = i + 1;

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

% I1 = SQRT2PIT*(FCm(j,i) - FCm(j-1,i) - FCp(j,i) + FCp(j-1,i) + ...
%      1i*(-FSm(j,i) + FSm(j-1,i) + FSp(j,i) - FSp(j-1,i)));
%  
% I2 = SQRT2PIT*Z*FCm(j,i) - SQRT2PIT*Z*FCm(j-1,i) + SQRT2PIT*Z*FCp(j,i) - ...
%      SQRT2PIT*Z*FCp(j-1,i) + 1i*(DT2*(Cm(j,i) - Cm(j-1,i)) + SQRT2PIT*Z*(-FSm(j,i) + FSm(j-1,i))) + ...
%      1i*(DT2*(-Cp(j,i) + Cp(j-1,i)) + SQRT2PIT*Z*(-FSp(j,i) + FSp(j-1,i))) + DT2*(Sm(j,i) - Sm(j-1,i)) + ...
%      DT2*(-Sp(j,i) + Sp(j-1,i));
end

