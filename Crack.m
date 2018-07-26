classdef Crack < Element
    %CRACK is the class for a crack in the system. It's a subclass of
    %Element.
    
    properties
        constants = [1,1];
        ends %a tuple
        p %the pressure in the crack
        locations %center of crack in z plane. Depends on ends.
        beta %the angle of the crack. Depends on ends.
        L %the length of the crack. Depends on the ends.
        dZdz %also dependent.
        Tau11
        Tau12
    end
    
    methods
        
        %Constructor
        function obj = Crack(ends,p)
            obj.ends = ends;
            obj.p = p;
            obj.locations = (ends(1)+ends(2))/2;
            obj.beta = angle(ends(2)-ends(1));
            obj.L = abs(ends(1)-ends(2));
            obj.dZdz = 2/obj.L*exp(-1i*obj.beta);
        end
        
        %Common methods
        function obj = solveforconstants(obj,tau11all, tau12all)%update this as needed for other elements
            %changes the value of obj.constants to reflect the updated
            %world.
            z = obj.locations;

            ToprowofA = [exp(-1i*obj.beta)*littlephiprime(obj,z), exp(-1i*obj.beta)*conj(littlephiprime(obj,z))];
            BottomrowofA = [0,exp(-1i*obj.beta)*1/(2*pi*1i)*F0prime(obj,z)];
            
            A = vertcat(ToprowofA, BottomrowofA);
            B = [-tau12all;obj.p-tau11all];
            X = A\B;
            
            obj.constants(1) = X(1);
            obj.constants(2) = X(2);

        end
        function Psi = calcPsi(obj, z)
            %calculates the Phi contribution from THIS ELEMENT ONLY
            c = obj.constants(1);
            Psi = c * bigPsi(obj, z);
        end
        function PhiPrime = calcPhiPrime(obj, z)
            %calculates the Phi contribution from THIS ELEMENT ONLY
            c = obj.constants(1);
            PhiPrime = c*bigPhiprime(obj,z);
        end
        function PhiBar = calcPhiBar(obj,z)
            %calculates the PhiBar contribution from THIS ELEMENT ONLY
            zbar = conj(z);
            PhiBar = obj.constants(1) * bigPhibar(obj,zbar);
        end    
        function Tau11 = calcTau11(obj, z,anglealpha)
            %calculates the Tau contribution from THIS ELEMENT ONLY
%             zbar = conj(z);
%             c = obj.constants(1) + 1i * obj.constants(2);
%             a = bigPhiprimeprime(obj,z);
%             b = bigPsiprime(obj,z);
%             Tau11 = -zbar*c*a+conj(c)*b;
            Z = ztoBigZ(obj,z);
            Zb = conj(Z);
            L = obj.L;
            c = obj.constants(1) + 1i * obj.constants(2);
            Tau11 = L/2*(Z-Zb)*c*littlephiprimeprime(obj, z)+exp(-1i*anglealpha)*(c*littlephiprime(obj,z)+conj(c)*littlepsiprime(obj,z));
            obj.Tau11 = Tau11;
        end       
        function Tau12 = calcTau12(obj, z,anglealpha)
            %calculates the Tau contribution from THIS ELEMENT ONLY
            %zbar = conj(z);
            c = obj.constants(1) + 1i * obj.constants(2);
            B1 = calcB1(obj,z);
            B1bar = conj(B1);
            Tau12 = c*exp(1i*anglealpha)*littlephiprime(obj,z)+exp(-1i*anglealpha)*conj(c)*conj(littlephiprime(obj,z));
            obj.Tau12 = Tau12;
        end
        function B1 = calcB1(obj,z)
            B1 = 0;
        end
        
        
        %Element-specific functions
        function Z = ztoBigZ(obj, z)
            z0 = obj.locations(1);
            Z = 2* exp(-1i*obj.beta)*(z-z0)/obj.L;
        end
        
        function F0z = F0(obj, z)
            Z = ztoBigZ(obj,z);
            F0z = log((Z-1)/(Z+1));
        end
        function F0zprime = F0prime(obj, z)
            Z = ztoBigZ(obj,z);
            F0zprime = ((1/(Z-1))-(1/(Z+1)))*obj.dZdz;
        end
        function F0zprimeprime = F0primeprime(obj, z)
            Z = ztoBigZ(obj,z);
            F0zprimeprime = -((Z-1)^-2-(Z+1)^-2)*(obj.dZdz)^2;
        end
        
        function lpsi = littlepsi(obj, z)
            %this one must be multiplied by the constant!
            lpsi = (1/(2i*pi))*F0(obj,z);
        end
        function lpsiprime = littlepsiprime(obj,z)
            %this one must be multiplied by the constant!
            lpsiprime = (1/2i*pi)*F0prime(obj,z);
        end
        function lpsiprime = littlepsiprimebar(obj,zbar)
            %this one must be multiplied by the constant!
            lpsiprime = (-1/2i*pi)*F0prime(obj,zbar);
        end
        
        function lphibar = littlephibar(obj,zbar)
            %must be multiplied by the constant!
            lphibar = (-1/(2i*pi))*F0(obj,zbar);
        end
        function lphiprime = littlephiprime(obj, z)
            %must be multiplied by the constant!
            lphiprime = (1/(2i*pi))*F0prime(obj, z);
        end
        function lphiprimebar = littlephiprimebar(obj, zbar)
            %must be multiplied by the constant!
            lphiprimebar = (1/(-2i*pi))*F0prime(obj, zbar);
        end
        function lphiprimeprime = littlephiprimeprime(obj, z)
            lphiprimeprime = (1/2i*pi)*F0primeprime(obj, z);
        end
        
        function bPhibar = bigPhibar(obj, zbar)
            %needs to be multiplied by the constant!
            bPhibar = exp(-1i*obj.beta)*littlephibar(obj, zbar);
        end
        function bPhiprime = bigPhiprime(obj, z)
            %needs to be multiplied by a constant!
            bPhiprime = exp(1i*obj.beta)*littlephiprime(obj, z);
        end
        function bPhiprimebar = bigPhiprimebar(obj, zbar)
            %needs to be multiplied by the constant!
            bPhiprimebar = exp(-1i*obj.beta)*littlephiprimebar(obj, zbar);
        end
        function bPhiprimeprime = bigPhiprimeprime(obj, z)
            %needs to be multiplied by the constant!
            bPhiprimeprime = exp(1i*obj.beta)*littlephiprimeprime(obj, z);
            
        end
        
        function bPsi = bigPsi(obj, z)
            %needs to be multiplied by the constant!
            Z = ztoBigZ(obj,z);
            zbar0 = conj(obj.locations);
            bPsi = obj.L/2*exp(-1i*obj.beta)*Z*bigPhiprime(obj, z)+zbar0*bigPhiprime(obj, z)+exp(-1i*obj.beta)*littlepsi(obj,z);
        end
        function bPsiprime = bigPsiprime(obj, z)
            %needs to be multiplied by the constant!
            Z = ztoBigZ(obj,z);
            zbar0 = conj(obj.locations);
            bPsiprime = obj.L/2*exp(-1i*obj.beta)*Z*obj.dZdz*bigPhiprimeprime(obj,z)+zbar0*bigPhiprimeprime(obj, z)+exp(-1i*obj.beta)*littlepsiprime(obj, z);
        end

    end
end

