function [FCm, FSm, FCp, FSp, Cm, Cp, Sm, Sp] = coord_arg_fresnel_for_IM(DZ, DT, Zout)

NZ = Zout/DZ;

SQRT2PIT = sqrt(2*pi*DT);

FCm = zeros(NZ,NZ);
FSm = zeros(NZ,NZ);
FCp = zeros(NZ,NZ);
FSp = zeros(NZ,NZ);
Cm = zeros(NZ,NZ);
Cp = zeros(NZ,NZ);
Sm = zeros(NZ,NZ);
Sp = zeros(NZ,NZ);
Jm = zeros(NZ,NZ);
Jp = zeros(NZ,NZ);
Km = zeros(NZ,NZ);
Kp = zeros(NZ,NZ);

i = 251;
% for i=1:NZ+1
    for j=1:NZ+1
        Jm(j,i) = ((j-i)*DZ)/SQRT2PIT;
        Km(j,i) = ((j-i)*DZ).^2/(4*DT);
        Jp(j,i) = ((j-2+i)*DZ)/SQRT2PIT;
        Kp(j,i) = ((j-2+i)*DZ).^2/(4*DT);
        FCm(j,i) = fresnelc(Jm(j,i));
        FSm(j,i) = fresnels(Jm(j,i));
        FCp(j,i) = fresnelc(Jp(j,i));
        FSp(j,i) = fresnels(Jp(j,i));
        Cm(j,i) = cos(Km(j,i));
        Cp(j,i) = cos(Kp(j,i));
        Sm(j,i) = sin(Km(j,i));
        Sp(j,i) = sin(Kp(j,i));        
    end  
% end
end

