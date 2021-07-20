function [time, bound] = ElutionModel(nCaptureAntibodies, nProtein, nSignalAntibodies, nCompetingAntibodies, nExcessProtein, valency, kon, koff)
    
    if numel(kon) == 1; kon = repmat(kon, 1, 3); end
    global Kon; Kon = kon; % [Kcapture, Ksignal, Kcompeting]
    if numel(koff) == 1; koff = repmat(koff, 1, 3); end
    global Koff; Koff = koff; % [Kcapture, Ksignal, Kcompeting]
    global pjump; pjump = [.8, .2, .2];
    
    % preallocate output arrays
    thresh = 0.1; % less than xx fraction of values will be less than tStep
    global tStep; tStep = expinv(thresh, 1/max(max(kon, koff))); 
    t = 0; tMax = 4000 * tStep;
    time = (0:tStep:tMax - tStep)';
    bound = zeros(numel(time),1);
    step = 1;
    
    % protein -> 0 = unbound, 1 = bound to capture Ab, 2 = bound to signal, 3 = bound to competing
    global protein; valency = round(sqrt(valency))^2; % note requirement: nProtein >= valency
    protein = zeros(nProtein + nExcessProtein, sqrt(valency), max(2,sqrt(valency))); % add excess protein
    protein(1:nProtein, 1, 1) = 1; % all bound to one capture Ab
    protein(1:nProtein, 1, end) = 2; % all bound to one signal Ab

    % captureAntibodies -> 0 = unbound, i = linear index of protein bound to
    captureAntibodies = zeros(nCaptureAntibodies,3);
    captureAntibodies(1:numel(find(protein == Antibodies.Capture)),1) = find(protein == Antibodies.Capture);
    captureAntibodies(:,2) = -1;
    % signalAntibodies -> 0 = unbound, i = linear index of protein bound to
    signalAntibodies = zeros(nSignalAntibodies,3);
    signalAntibodies(1:numel(find(protein == Antibodies.Signal)),1) = find(protein == Antibodies.Signal);
    signalAntibodies(:,2) = -1;
    % competingAntibodies -> 0 = unbound, i = linear index of protein bound to
    competingAntibodies = zeros(nCompetingAntibodies,3);
    competingAntibodies(:,2) = -1;
    
    while t < tMax
        
        % record number of bound signalAntibodies
        if step > numel(bound); break; end
        bound(step) = sum(sum(protein == 2,[2 3]) .* any(protein == 1,[2 3]));
        
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
    global Kon; global Koff; global pjump;
    current = 1; goal = 2; time = 3;
    
    if ab(i,goal) ~= -1 % if action is in progress continue

        ab(i,time) = ab(i,time) - tStep;
        if ab(i,time) <= 0 % finished, complete and reset
            
            if ab(i,current) > 0; protein(ab(i,current)) = 0; end % unbind the protein as well
            if ab(i,current) > 0 && ab(i,goal) == 0 && rand < pjump(i_ab)
                [l, w, h] = size(protein);
                [x, y, z] = ind2sub([l,w,h], ab(i,current));
                wrapN = @(x, N) (1 + mod(x-1, N));
                for j = -1:1
                    for k = -1:1
                        dy = wrapN(y+j,w); dz = wrapN(z+k,h);
                        if protein(x, dy, dz) == 0
                            protein(x, dy, dz) = i_ab;
                            ab(i,current) = sub2ind([l,w,h], x, dy, dz);
                            ab(i,goal) = -1;
                            return;
                        end
                    end
                end
            else
                ab(i,current) = ab(i,goal);
                ab(i,goal) = -1;
            end
            
        end

    elseif ab(i,current) == 0 % if currently unbound

        % bind antibody
        free_proteins = find(protein == 0);
        cap_bind_thresh = round(numel(protein)/2);
        if i_ab == Antibodies.Capture; free_proteins = free_proteins(free_proteins <= cap_bind_thresh); % capture abs only bind to first position, ensure there is only one capture Ab per protein
        else; free_proteins = free_proteins(free_proteins > cap_bind_thresh); end % all other abs bind to not first position
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