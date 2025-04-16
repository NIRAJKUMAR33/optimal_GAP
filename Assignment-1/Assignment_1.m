function optimal_large_gap_modified()
    num_files = 12;
    all_results = cell(num_files, 1);
    all_objective_values = [];

    % Iterate through gap1 to gap12
    for g = 1:num_files
        filename = sprintf('gap%d.txt', g);
        fid = fopen(filename, 'r');
        if fid == -1
            error('Error opening file %s.', filename);
        end

        % Read the number of problem sets
        num_problems = fscanf(fid, '%d', 1);
        results = cell(num_problems, 1);

        for p = 1:num_problems
            % Read problem parameters
            m = fscanf(fid, '%d', 1);
            n = fscanf(fid, '%d', 1);

            % Read cost and resource matrices
            c = fscanf(fid, '%d', [n, m])';
            r = fscanf(fid, '%d', [n, m])';

            % Read server capacities
            b = fscanf(fid, '%d', [m, 1]);

            % Solve the problem
            x_matrix = solve_gap_max(m, n, c, r, b);
            objective_value = sum(sum(c .* x_matrix));

            results{p} = sprintf('c%d-%d\t%d', m*100 + n, p, round(objective_value));
            all_objective_values = [all_objective_values; objective_value];
        end

        fclose(fid);
        all_results{g} = results;
    end

    % Display results side by side
    files_per_row = 4;
    for start_idx = 1:files_per_row:num_files
        end_idx = min(start_idx + files_per_row - 1, num_files);

        for g = start_idx:end_idx
            fprintf('gap%d\t\t', g);
        end
        fprintf('\n');

        max_problems = max(cellfun(@length, all_results(start_idx:end_idx)));

        for p = 1:max_problems
            for g = start_idx:end_idx
                if p <= length(all_results{g})
                    fprintf('%s\t', all_results{g}{p});
                else
                    fprintf('\t\t');
                end
            end
            fprintf('\n');
        end
        fprintf('\n');
    end

    % Plot the optimal fitness values
    figure;
    plot(1:length(all_objective_values), all_objective_values, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
    title('Optimal Fitness Value per Problem Instance');
    xlabel('Instance Index');
    ylabel('Optimal Fitness Value');
    grid on;
end

% ✅ This is the sub-function that solves the GAP problem
function x_matrix = solve_gap_max(m, n, c, r, b)
    f = -c(:); % Convert to column vector for maximization

    % Constraint 1: Each user assigned exactly once
    Aeq_jobs = kron(eye(n), ones(1, m));
    beq_jobs = ones(n, 1);

    % Constraint 2: Server resource constraints
    Aineq_agents = zeros(m, m * n);
    for i = 1:m
        for j = 1:n
            Aineq_agents(i, (j-1)*m + i) = r(i,j);
        end
    end
    bineq_agents = b;

    % Define variable bounds (binary decision variables)
    lb = zeros(m * n, 1);
    ub = ones(m * n, 1);
    intcon = 1:(m*n);

    % Solve using intlinprog
    options = optimoptions('intlinprog', 'Display', 'off');
    x = intlinprog(f, intcon, Aineq_agents, bineq_agents, Aeq_jobs, beq_jobs, lb, ub, options);

    % Reshape into m × n matrix
    x_matrix = reshape(x, [m, n]);
end
optimal_large_gap_modified();
