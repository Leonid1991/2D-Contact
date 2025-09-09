function Name = approachName(approach)

    names = ["None", "Penalty", "Nitsche Linear", ...
             "Nitsche for Solids", "Nitsche", ...
             "Lagrange multiplier", "Penalty (simplified)"];
        if approach >= 0 && approach <= numel(names)-1
            Name = names(approach+1);
        else
            Name = "NaN";
        end
    
    end    