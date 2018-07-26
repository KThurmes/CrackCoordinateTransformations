classdef HalfSpace < Element
    %HALFSPACE is the half-space element
    
    properties
        constants %this is going to be a matrix with (degree) columns and 2 rows.
        locations %locations around the circle where we will compute the taus
        chilocations %the chi version of locations
        d %the value used to comput chi
        degree %the degree of the taylor series
        N %The number of locations
        values %I don't really know how to describe this. It's an m by n 
        %matrix that holds values for exp(imtheta) values. It seemed 
        %silly to recompute them every time I went through this.
        beta = 0;
    end
    
    methods
        function obj = HalfSpace(d, N, degree)
            obj.d = d;
            obj.N = N;
            obj.degree = degree;
            obj.values = zeros(N,degree+1);
            obj = populateLocations(obj);
            obj.constants = zeros(2,obj.degree+1);
        end %done
        function obj = populateLocations(obj)
            N = obj.N;
            %determine delta theta
            dtheta = (2*pi)/obj.N;
            %put the locations in their place
            for i = 1:N
                %find the z locations of each chi location
                Chi = exp(1i*dtheta*i);
                z = chitoz(obj, Chi);
                obj.locations(i) = z;
                obj.chilocations(i) = Chi;
                %let's throw in the exp(imtheta) part, too, while
                %we're here.
                for j = 0:obj.degree
                    c = exp(1i*dtheta*i*j);
                    obj.values(i,j+1)=c;
                end
            end
        end %done
        function Chi = ztoChi(obj, z)
            %chi is a variable in a conformally mapped plane. On this
            %plane, the surface of the half space maps to a unit circle.
            Z = z/obj.d;
            Chi = (Z+1i)/(Z-1i);
        end %done
        function z = chitoz(obj, Chi)
            %chi is a variable in a conformally mapped plane. On this
            %plane, the surface of the half space maps to a unit circle.
            z = (1i*obj.d*(1+Chi))/(Chi-1);
        end %done
        function obj = solveforconstants(obj, tau11all, tau12all)%update this as needed for other elements
            %change inputs when I know what info I need from other elements
            %This function returns alphas and betas
            
            %check to see if eimthetas has something in it. if not, go
            %ahead and populate it.
            
            dtheta = (2*pi)/obj.N;
            %now start calculating each betam!
            betams = zeros(1,obj.degree+1);
            alphams = zeros(1,obj.degree+1);
            
            %sum the values of betam and alpham over the number of degrees.
            for deg = 0:obj.degree
                betam = 0;
                alpham = 0;
                
                %sum the components of betam and alpham over 2pi.
                for i = 1:length(obj.locations)
                    tau11others = tau11all(i);
                    tau12others = tau12all(i);
                    T = .5*1i*(tau11others-tau12others); %traction due to all other elements on a horizontal plane.
                    tothers1 = real(T); %The t1 values are the real parts of the traction
                    tothers2 = -imag(T); %the t2 values are the imaginary parts of the traction.
                    
                    betam = betam + (tothers2 * obj.values(i,deg+1));
                    alpham = alpham + (tothers1 /obj.values(i,deg+1));
                end
                
                betam = dtheta * betam*1i/pi;
                betams(deg+1) = betam;
                
                alpham = dtheta * alpham / pi;
                alphams(deg+1) = alpham;

            end
            %when we have all the alphams and betams in two nice little
            %arrays, make them into the constants property.
            
            obj.constants = vertcat(alphams, betams);
        end
        
        function Psi = calcPsi(obj,z)
        end
        
        function PhiBar = calcPhiBar(obj,z)
        end

        function Tau11 = calcTau11(obj,z,anglealpha)
            Tau11 = (z-conj(z))*calcPhiPrimePrime(obj,z)+calcPhiPrime(obj,z)+calclittlepsiprime(obj,z);
            Tau11 = Tau11*exp(-1i*anglealpha);
        end%done        
        function PhiPrimePrime = calcPhiPrimePrime(obj,z) %done
            PhiPrimePrime = 0.5*(calcPPrime(obj,z)+calcQPrime(obj,z));
        end
        function Tau12 = calcTau12(obj,z,anglealpha)
            Tau12 = calcPhiPrime(obj,z)+calcPhiPrimeBar(obj,z);
        end %done
        function PhiPrimeBar = calcPhiPrimeBar(obj,z) %done
            p = calcPhiPrime(obj,conj(z));
            PhiPrimeBar = conj(p);
        end
        function B1 = calcB1(obj,z)
            B1 = 0;
        end%done
        function PhiPrime = calcPhiPrime(obj,z)
            PhiPrime = 0.5*(Pz(obj,z)+Qz(obj,z));
        end %done
        function lpsiprime = calclittlepsiprime(obj,z)
            lpsiprime = 0.5*(Pz(obj,z)-Qz(obj,z));
        end %done
        function dChidz = dChidz(obj,Z)
            dChidz = (-(Z+i)/(Z-i)^2+(Z-1)^-1)*1/obj.d;
        end %done
        function P = Pz(obj, z)
            P = 0;
            for i = 1:obj.degree+1
                chin = ztoChi(obj,z);
                alphan = obj.constants(1,i);
                P = P + alphan * chin^(i-1);
            end
        end %done
        function Q = Qz(obj,z)
            Q = 0;
            for i = 1:obj.degree+1
                chin = ztoChi(obj,z);
                betan = obj.constants(2,i);
                Q = Q + betan * chin^(i-1);
            end
        end%done
        function PPrime = calcPPrime(obj,z)
            PPrime = 0;
            for i = 1:obj.degree+1
                chin = ztoChi(obj,z);
                alphan = obj.constants(1,i);
                PPrime = PPrime + alphan*(i-1)*chin^(i-1)*dChidz(obj,z);
            end
        end%done
        function QPrime = calcQPrime(obj,z)
            QPrime = 0;
            for i = 1:obj.degree+1
                chin = ztoChi(obj,z);
                betan = obj.constants(2,i);
                QPrime = QPrime + betan*(i-1)*chin^(i-1)*dChidz(obj,z);
            end
        end%done
    end
    
end

