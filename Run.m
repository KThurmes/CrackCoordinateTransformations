clear
clc
z = -2.5-2.5i;
anglealpha = 0.7854;
rho = 5;
nu = 0.3; %typical for shale.
E = 35; %Typical value for shale: 1-70 GPa (Young's Modulus)
G = E/(2*(1+nu)); %Shear modulus. In GPa.
gravity = 9.81;
thisworld = World(G, rho, nu);
gravy = Gravity(rho, gravity, nu);
lstart = -4-4i;
lend = -1-1i;
cracken = Crack([-4-4i,-1+-1i],100);
cadet = HalfSpace(1,10,5);
thisworld = addelement(thisworld, gravy);
thisworld = addelement(thisworld, cracken);
thisworld = addelement(thisworld, cadet);
thisworld = solveworld(thisworld);


%Moment of truth

% T = traction(thisworld,z,anglealpha);
 normaltocrack = ContourMe_I_nint(-10,10,100,-10,10,100,@(z)traction(thisworld,z,anglealpha),10);
 line([real(lstart),real(lend)],[imag(lstart),imag(lend)]);
 tangenttocrack = ContourMe_R_nint(-10,10,100,-10,10,100,@(z)traction(thisworld,z,anglealpha),10);
 line([real(lstart),real(lend)],[imag(lstart),imag(lend)]);

