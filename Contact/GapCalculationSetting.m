function Gapfunc = GapCalculationSetting(GapCalculation, n)
        
    addpath("Contact\GapFunctions")

    if GapCalculation == "nodes"

        Gapfunc = @GapfuncNodes;

    elseif GapCalculation == "Gauss"

        Gapfunc = @(ContactBody,TargetBody) GapfuncGauss(ContactBody,TargetBody,n);

    elseif GapCalculation == "LinSpace"
        
        if n < 2 
           warning("Number of points is not enough, it is set to 2 ")
           n = 2; 
        end 

        Gapfunc = @(ContactBody,TargetBody) GapfuncLinSpace(ContactBody,TargetBody,n);

    else

        warning("Gap-function is set to node-based calculations")
        Gapfunc  = @GapfuncNodes;

        
    end    