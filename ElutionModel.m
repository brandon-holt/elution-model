function [time, bound] = ElutionModel(nCaptureAntibodies, nProtein, nSignalAntibodies, nCompetingAntibodies, nExcessProtein, kon, koff)
    
    if numel(kon) == 1; kon = repmat(kon, 1, 3); end
    global Kon; Kon = kon; % [Kcapture, Ksignal, Kcompeting]
    if numel(koff) == 1; koff = repmat(koff, 1, 3); end
    global Koff; Koff = koff; % [Kcapture, Ksignal, Kcompeting]
    
    % preallocate output arrays
    global tStep; tStep = .1;
    t = 0; tMax = 30;
    time = (0:tStep:tMax - tStep)';
    bound = zeros(numel(time),1);
    step = 1;
    
    % protein -> 0 = unbound, 1 = bound to capture Ab, 2 = bound to signal, 3 = bound to competing
    global protein; valency = 60; % note requirement: nProtein >= valency
    protein = zeros(nProtein + nExcessProtein, valency); % add excess protein
    protein(1:nProtein, 1) = 1; % all bound to one capture Ab
    protein(1:nProtein, 2) = 2; % all bound to one signal Ab

    % captureAntibodies -> 0 = unbound, i = linear index of protein bound to
    captureAntibodies = zeros(nCaptureAntibodies,3);
    captureAntibodies(1:min(nProtein, nCaptureAntibodies),1) = 1:min(nProtein, nCaptureAntibodies);
    captureAntibodies(:,2) = -1;
    % signalAntibodies -> 0 = unbound, i = linear index of protein bound to
    signalAntibodies = zeros(nSignalAntibodies,3);
    signalAntibodies(1:min(nProtein, nSignalAntibodies),1) = length(protein) + (1:min(nProtein, nSignalAntibodies));
    signalAntibodies(:,2) = -1;
    % competingAntibodies -> 0 = unbound, i = linear index of protein bound to
    competingAntibodies = zeros(nCompetingAntibodies,3);
    competingAntibodies(:,2) = -1;
    
    while t < tMax
        
        % record number of bound signalAntibodies
        if step > numel(bound); break; end
        bound(step) = sum(sum(protein == 2,2) .* any(protein == 1,2));
        
        % update state of all antibodies
        for i = randperm(max([length(signalAntibodies), length(competingAntibodies), length(captureAntibodies)]))
            captureAntibodies = UpdateAntibodyState(captureAntibodies, i, Antibodies.Capture);
            signalAntibodies = UpdateAntibodyState(signalAntibodies, i, Antibodies.Signal);
            competingAntibodies = UpdateAntibodyState(competingAntibodies, i, Antibodies.Competing);
        end

        % increment time
        t = t + tStep;
        step = step + 1;
        
    end

end

function [ab] = UpdateAntibodyState(ab, i, i_ab)
    
    if i > length(ab); return; end
    
    global protein; global tStep;
    global Kon; global Koff;
    current = 1; goal = 2; time = 3;
    
    if ab(i,goal) ~= -1 % if action is in progress continue

        ab(i,time) = ab(i,time) - tStep;
        if ab(i,time) <= 0 % finished, complete and reset
            if ab(i,current) > 0; protein(ab(i,current)) = 0; end % unbind the protein as well
            ab(i,current) = ab(i,goal);
            ab(i,goal) = -1;
        end

    elseif ab(i,current) == 0 % if currently unbound

        % bind antibody
        free_proteins = find(protein == 0);
        if i_ab == Antibodies.Capture; free_proteins = free_proteins(free_proteins <= length(protein)); % capture abs only bind to first position, ensure there is only one capture Ab per protein
        else; free_proteins = free_proteins(free_proteins > length(protein)); end % all other abs bind to not first position
        if numel(free_proteins) == 0; return; end
        free_protein = free_proteins(randi([1, length(free_proteins)]));
        ab(i,goal) = free_protein;
        ab(i,time) = exprnd(1 / Kon(i_ab));
        % bind protein
        protein(free_protein) = i_ab;

    elseif ab(i,current) > 0 % if currently bound, begin unbinding

        ab(i,goal) = 0;
        ab(i,time) = exprnd(1 / Koff(i_ab));

    end

end