classdef Gravity < Element
    %GRAVITY In a world with a half-space boundary
    
    properties
        constants = [];
        locations = [];
        rho
        gravity
        kappa
        nu
        beta = 0;
        Tau12
        Tau11
    end
    
    methods
        
        function obj = Gravity(rho,g,nu)
            obj.rho = rho;
            obj.gravity = g;
            obj.nu = nu;
            obj.kappa = 3-4*obj.nu;
        end
        function obj = solveforconstants(obj, tau11all, tau12all)%update as needed for other elements
        end
        %change inputs when I know what info i need from other elements
        function Psi = calcPsi(obj, z)
            zbar = conj(z);
            Psi = z*calcPhiPrime(obj,z)+1i*obj.rho*obj.gravity*(1-2*obj.nu)/(obj.kappa+1)*zbar^2;
        end
        function PhiPrime = calcPhiPrime(obj, z)
            PhiPrime = -2i * obj.rho * obj.gravity * z * (1-2*obj.nu)/(obj.kappa + 1);
        end
        function PhiBar = calcPhiBar(obj, z)
            zbar = conj(z);
            PhiBar = 1i*obj.rho*obj.gravity*zbar^2*(1-2*obj.nu)/(obj.kappa+1);
        end
        function Tau11 = calcTau11(obj,z,anglealpha)
            zbar = conj(z);
            Tau11 = 2i*(1-2*obj.nu)/(obj.kappa+1)*(z-zbar);
            obj.Tau11 = Tau11*exp(-1i*anglealpha);
        end
        function Tau12 = calcTau12(obj, z,anglealpha)
            zbar = conj(z);
            Tau12 = 2i*obj.rho*obj.gravity*1/(obj.kappa+1)*(z-zbar);
            obj.Tau12 = Tau12;
        end
        function B1 = calcB1(obj,z)
            B1 = obj.rho*obj.gravity*z;
        end
    end
    
end

