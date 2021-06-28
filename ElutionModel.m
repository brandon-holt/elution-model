function [time, bound] = ElutionModel(method, nSignalAntibodies, nProtein, nSolution, kon, koff)

    % kon/koff = 1x2, element-1 = for signal antibody, element-2 = for
    % competing antibody. can be 1x1 for Methods.ExcessProtein
    global Kon; Kon = kon;
    global Koff; Koff = koff;
    global tStep; tStep = 0.01;
    t = 0;
    tMax = 100;
    
    % starting at moment where all signal antibodies are bound to surface
    % protein. 3 columns: 1 = current state, 2 = goal state, 3 = time
    % remaining
    if nProtein < nSignalAntibodies; nProtein = nSignalAntibodies; end
    global antibodies;
    antibodies = zeros(nSignalAntibodies,3);
    antibodies(:,1) = 1;
    antibodies(:,2) = -1;
    global surface;
    surface = zeros(nProtein,1);
    surface(1:nSignalAntibodies, 1) = 1;
    global solution;
    solution = zeros(nSolution,3);
    solution(:,2) = -1;
    
    % preallocate outputs
    time = (0:tStep:tMax - tStep)';
    bound = zeros(numel(time),1);
    step = 1;
    
    while t < tMax
        
        % record number of bound antibodies
        bound(step) = sum(antibodies(:,1) == 1);
        
        % update state based on method
        if method == Methods.ExcessProtein
            ExcessProteinMethod();
        elseif method == Methods.CompetingAntibody
            CompetingAntibodyMethod();
        end
        
        % increment time
        t = t + tStep;
        step = step + 1;
        
    end

end

function ExcessProteinMethod()

    global Kon;
    global Koff; 
    global tStep; 
    global antibodies;
    global surface;
    global solution;

    for i = 1:length(antibodies)

        if antibodies(i,2) ~= -1 % if action is in progress continue
            
            antibodies(i,3) = antibodies(i,3) - tStep;
            if antibodies(i,3) <= 0 % finished, complete and reset
                if antibodies(i,2) == 0 % unbind protein             
                    if antibodies(i,1) == 1; surface(find(surface(:,1) == 1, 1, 'first'),1) = 0; end
                    if antibodies(i,1) == 2; solution(find(solution(:,1) == 1, 1, 'first'),1) = 0; end
                end
                antibodies(i,1) = antibodies(i,2);
                antibodies(i,2) = -1;
            end

        elseif antibodies(i,1) == 0 % if currently unbound
            
            % bind antibody
            unbound_1 = sum(surface(:,1) == 0);
            unbound_2 = sum(solution(:,1) == 0);
            p_bind_1 = unbound_1 / (unbound_1 + unbound_2);
            p_bind_2 = unbound_2 / (unbound_1 + unbound_2);
            cum_dist = cumsum([p_bind_1, p_bind_2]);
            rand_protein = find(rand < cum_dist, 1, 'first');
            if numel(rand_protein) == 0; continue; end
            antibodies(i,2) = rand_protein;
            antibodies(i,3) = exprnd(1 / Kon(1));
            % bind protein
            if rand_protein == 1; surface(find(surface(:,1) == 0, 1, 'first'),1) = 1; end
            if rand_protein == 2; solution(find(solution(:,1) == 0, 1, 'first'),1) = 1; end

        elseif antibodies(i,1) == 1 || antibodies(i,1) == 2 % if currently bound
            
            % begin unbinding antibody
            antibodies(i,2) = 0;
            antibodies(i,3) = exprnd(1 / Koff(1));

        end

    end
    
end

function CompetingAntibodyMethod()

    global Kon;
    global Koff; 
    global antibodies;
    global solution;
    
    max_count = max(length(antibodies), length(solution));
    
    for i = 1:max_count
        
        if i < length(antibodies)
            antibodies = ProcessAntibodyOneTarget(antibodies, i, Kon(1), Koff(1));
        end
        
        if i < length(solution)
            solution = ProcessAntibodyOneTarget(solution, i, Kon(2), Koff(2));
        end
        
    end

end

function [array] = ProcessAntibodyOneTarget(array, i, Kon, Koff)
    
    global surface;
    global tStep; 
    
    if array(i,2) ~= -1 % if action is in progress continue

        array(i,3) = array(i,3) - tStep;
        if array(i,3) <= 0 % finished, complete and reset
            array(i,1) = array(i,2);
            array(i,2) = -1;
            if array(i,1) == 0 % unbind protein
                surface(find(surface(:,1) == 1, 1, 'first'),1) = 0;
            end
        end

    elseif array(i,1) == 0 % if currently unbound

        % bind antibody
        free_protein = find(surface(:, 1) == 0, 1, 'first');
        if numel(free_protein) == 0; return; end
        array(i,2) = 1;
        array(i,3) = exprnd(1 / Kon);
        % bind protein
        surface(free_protein,1) = 1;

    elseif array(i,1) == 1 % if currently bound

        % begin unbinding antibody
        array(i,2) = 0;
        array(i,3) = exprnd(1 / Koff);

    end

end