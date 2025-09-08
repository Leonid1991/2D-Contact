function Gapfunc = GapCalculationSetting(GapCalculation, n)

    if GapCalculation == "nodes"
        Gapfunc = @GapfuncNodes;

    elseif GapCalculation == "Gauss"

        Gapfunc = @(ContactBody,TargetBody) GapfuncGauss(ContactBody,TargetBody,n);

    end    