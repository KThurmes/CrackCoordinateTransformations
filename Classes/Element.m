classdef Element
     %Abstract class for creating components of a world
    properties (Abstract = true)
        constants
        locations
        beta
    end
    
    methods (Abstract = true)
        obj = solveforconstants(obj, tau11all, tau12all )%update this as needed for other elements
        %change inputs when I know what info i need from other elements
        Psi = calcPsi(obj,z)
        PhiPrime = calcPhiPrime(obj,z)
        PhiBar = calcPhiBar(obj,z)
        Tau11 = calcTau11(obj,z, anglealpha)
        Tau12 = calcTau12(obj,z, anglealpha)
        B1 = calcB1(obj,z)
        %not entirely sure what to do with this one. Will have to see.
    end
    
end

