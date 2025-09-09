function Body = SaveResults(Body, step, something)

    % Validate input
    if ~(isnumeric(something) || (isstring(something) && (something == "last" || something == "all")))
        error("Provide correct identifier: must be numeric, 'last', or 'all'");
    end

    % Initialize or reset results
    if ~isfield(Body, 'results') || (isstring(something) && something == "last")
        Body.results = [];
    end

    % Default: do not save
    flag = false;

    % Save logic
    if isstring(something) && something == "all"
        flag = true;
    elseif isnumeric(something) && mod(step, something) == 0
        flag = true;
    elseif isstring(something) && something == "last"
        flag = true;  % this only works if called at the last step externally
    end

    % Save data
    if flag
        NumberOfIDs = length(Body.Fext.nodalid);
        DofsAtNode = Body.DofsAtNode;            
        ux = Body.u(xlocChosen(DofsAtNode, Body.Fext.nodalid, 1));
        uy = Body.u(xlocChosen(DofsAtNode, Body.Fext.nodalid, 2));
        uxavg = sum(ux) / NumberOfIDs;
        uyavg = sum(uy) / NumberOfIDs;
        Body.results = [Body.results; Body.nElems.x, Body.nElems.y, Body.ndof, uxavg, uyavg];
    end

