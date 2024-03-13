clear all;

Nz = 500;
DT = 0.1;
Zout = 47;
DZ = Zout/Nz
ZAxis = 0:Zout/Nz:Zout;
I = 250;
Z = I*DZ;
% Z=23.500000000000000;
% f = @(z) exp(-3*((z - Zout/2)./(Zout/2)).^2) + 1i*exp(-3*((z - Zout/2)./(Zout/2)).^2);
f = @(z) (cos(2*pi*z/Zout).^2 + 1i*sin(2*pi*z/Zout).^2);
fg = @(z) f(z).*G(z,Z,DT);

fs = f(ZAxis);

for i=1:Nz+1
    yfs(i) = fs(i)*G(ZAxis(i),Z,DT);
    yfg(i) = fg(ZAxis(i));
end


Ii = integral(fg,0,Zout)
It = trapz(ZAxis,fg(ZAxis))

[FCm, FSm, FCp, FSp, Cm, Cp, Sm, Sp] = coord_arg_fresnel_for_IM(DZ, DT, Zout);

I1_old=zeros(Nz,1);
I2_old=zeros(Nz,1);
I1=zeros(Nz,1);
I2=zeros(Nz,1);
f0=zeros(Nz,1);
f1=zeros(Nz,1);
for j=1:Nz       
    [I1(j), I2(j)] = I12_IM(j, I, DZ, DT, FCm, FSm, FCp, FSp, Cm, Cp, Sm, Sp);
    [I1_old(j), I2_old(j)] = I12_IM_old(j, I, DZ, DT);
    f0(j) = j.*fs(j)-(j-1).*fs(j+1);
    f1(j) = (fs(j+1)-fs(j))/DZ;
end

Is = 0;
Is_old = 0;
for i=1:Nz
    Is = Is + f0(i)*I1(i) + f1(i)*I2(i);
    Is_old = Is_old + f0(i)*I1_old(i) + f1(i)*I2_old(i);
end
Is = Is*1/2*sqrt(1i/(pi*DT))
Is_old = Is_old*1/2*sqrt(1i/(pi*DT))

plot(ZAxis, yfs, ZAxis, yfg)

