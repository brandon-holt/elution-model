function [time, bound] = ElutionModel(nCaptureAntibodies, nProtein, nSignalAntibodies, nCompetingAntibodies, nExcessProtein, kon, koff)

    % element-1 = target protein + capture Ab
    % element-2 = target protein + signal Ab
    % element-3 = target protein + competing Ab
    global Kon; Kon = kon;
    global Koff; Koff = koff;
    
    % preallocate outputs
    global tStep; tStep = 0.01;
    t = 0; tMax = 30;
    time = (0:tStep:tMax - tStep)';
    bound = zeros(numel(time),1);
    step = 1;
    
    % starting at moment where all nProtein are bound to signal and capture Ab
    % 0 = unbound, 1 = bound to capture Ab, 2 = bound to signal, 3 = bound to competing
    global protein; valency = 2;
    protein = zeros(nProtein + nExcessProtein, valency);
    protein(1:nProtein, 1) = 1;
    protein(1:nProtein, 2) = 2;
    % for antibodies 3 columns: 1 = current state, 2 = goal state, 3 = time remaining
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
    
    if ab(i,2) ~= -1 % if action is in progress continue

        ab(i,3) = ab(i,3) - tStep;
        if ab(i,3) <= 0 % finished, complete and reset
            if ab(i,1) > 0; protein(ab(i,1)) = 0; end % unbind the protein as well
            ab(i,1) = ab(i,2);
            ab(i,2) = -1;
        end

    elseif ab(i,1) == 0 % if currently unbound

        % bind antibody
        free_proteins = find(protein == 0);
        if i_ab == Antibodies.Capture; free_proteins = free_proteins(free_proteins <= length(protein)); end % capture abs only bind to first position
        if numel(free_proteins) == 0; return; end
        free_protein = free_proteins(randi([1, length(free_proteins)]));
        ab(i,2) = free_protein;
        ab(i,3) = exprnd(1 / Kon(i_ab));
        % bind protein
        protein(free_protein) = i_ab;

    elseif ab(i,1) > 0 % if currently bound, begin unbinding

        ab(i,2) = 0;
        ab(i,3) = exprnd(1 / Koff(i_ab));

    end

end