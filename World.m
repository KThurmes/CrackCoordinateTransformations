classdef World
    %WORLD A world containing a set of elements
    
    properties
        gravity = 9.81; %m/s^2
        G %Define as something reasonable
        rho %define as the density of some rock
        nu
        elements = {}
        kappa
    end
    
    methods
        %constructor
        function obj = World(G,rho,nu)
            obj.G = G;
            obj.rho = rho;
            obj.nu = nu;
            obj.kappa = 3-4*nu;
        end
        function obj = addelement(obj,value)
            obj.elements{end+1} = value;
        end
        function obj = solveworld(obj)
            numelements = length(obj.elements);
            change = 1;
            while change >= 1
                change = 0;
                for i = 1:numelements
                    %designate one element to be "active"
                    activeelement = obj.elements{i};
                    
                    %create an array of "inactive elements"
                    inactiveelements = {};
                    for j = 1:numelements
                        if j ~= i
                            inactiveelements{end+1} = obj.elements{j};
                        end
                    end
                    
                    %get an array of "locations" from the active element
                    activelocations = activeelement.locations;
                    activebeta = activeelement.beta;
                    numlocations = length(activelocations);
                    
                    %ARRAYS OF RELEVANT CONTRIBUTIONS. WILL NEED MORE!
                    tau11s = zeros(numlocations);
                    tau12s = zeros(numlocations);
                    
                    for locindex = 1:numlocations
                        
                        loc = activelocations(locindex);
                        
                        [tau11, tau12] = calctausum(inactiveelements, loc, activebeta);
                        
                        tau11s(locindex) = tau11;
                        tau12s(locindex) = tau12;
                        
                    end
                    
                    %save a snapshot of the old constants
                    oldconstants = activeelement.constants;
                    
                    %pass those arrays in to the active element's "solve"
                    %function
                    activeelement = solveforconstants(activeelement, tau11s, tau12s);
                    obj.elements{i}=activeelement;
                    %look at new constants
                    newconstants = activeelement.constants;
                    
                    %compare the new constants to the old constants. Have
                    %some sort of threshold for determining whether it's
                    %converged or not.
                    %iterate through each element in constants to see if
                    %they differ by a set amount
                    for g = 1:size(newconstants,1)
                        for h = 1:size(newconstants,2)
                            
                            %****add a part for if there are two rows of constants.
                            oldc = oldconstants(g,h);
                            newc = newconstants(g,h);
                            
                            if abs(oldc-newc)>0.0001
                                change = change+1;
                            end
                        end
                    end
                end
            end
        end %function
        function abstaus = calcabstaus(obj,elementlist, z)
            abstaus = z;
        end
        function w = calcw(obj, elementlist, z)
            w = z;
        end
        function T = traction(obj, z, anglealpha)
            [tau11s, tau12s] = calctausum(obj.elements, z, anglealpha);
            T = .5*1i*(tau11s*exp(2*1i*anglealpha)-tau12s);
        end
    end
    methods (Static)
        function [ tau11s, tau12s ] = calctausum( elementlist, z, anglealpha)
            %CALCTAUSUM will calculate the taus for all elements in a list of elements
            
            numelements = length(elementlist);
            tau12s = 0;
            tau11s = 0;
            for i = 1:numelements
                activeelement = elementlist{i};
                tau12indiv = calcTau12(activeelement, z, anglealpha);
                tau12s = tau12s + tau12indiv;
                tau11indiv = calcTau11(activeelement, z, anglealpha);
                tau11s = tau11s + tau11indiv;
            end
            
        end
    end
end %class



