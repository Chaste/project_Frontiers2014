function print_timesteps_for_latex
% This is just to have a look at the models and timesteps we used

close all
clear all

d = importdata('required_steps.txt');

% Index in this vector (minus 1) is the solver code in the results files.
solver_mapping = {'CVODE (analytic Jacobian)',...
    'CVODE (numerical Jacobian)',...
    'Forward Euler',...
    'Backward Euler',...
    'Runge-Kutta (2nd order)',...
    'Runge-Kutta (4th order)',...
    'Rush-Larsen',...
    'Generalised Rush-Larsen 1',...
    'Generalised Rush-Larsen 2'};

print_for_optimised(d, false)

% You can run this to get a printout of this table too,
% But you will also get a report of the differences in a file called
% different_steps.txt
% There were only three entries (two on CVODE tolerances) so we forget
% about this for now.
% fprintf('\n\n')
% 
% print_for_optimised(d, true)

end

function print_for_optimised(d, optimised_or_not)

% Load up the difficulty ordering from the timing_analysis.m script.
% Makes a variable called 'ordering'.
load('../results/difficulty_ordering.mat');

% Get the raw data out
model = d.textdata;
solver = d.data(:,1);
optimised = d.data(:,2);
timesteps = d.data(:,3);
error_metrics = d.data(:,4:11);
converged = d.data(:,12);

assert(length(model)==length(solver));
assert(length(model)==length(optimised));
assert(length(model)==length(timesteps));

model_list = unique(model);
solver_list = unique(solver);
optimised_list = unique(optimised);
converged_list = unique(converged);

% Check these are binary
assert(length(optimised_list)==2)
assert(length(converged_list)==2)

% Only look at converged answers
converged_indices = find(converged==1);

converged_error_metrics = error_metrics(converged_indices,2:end-1);
disp('Mean absolute error metrics')
mean(abs(converged_error_metrics))

fid = fopen('different_steps.txt', 'at');

% Reorder the model list
model_list = model_list(ordering);

for m = 1:length(model_list)
    model_indices = find(strcmp(model,model_list(m)));
    % Always consider the converged cases only...
    model_indices = intersect(model_indices,converged_indices);
    
    fprintf('%s',model_list{m})
    for s = 1:length(solver_list)
        solver_indices = find(solver==solver_list(s));
        model_and_solver = intersect(model_indices,solver_indices);
        
        % With or without lookup tables
        lookup_indices = find(optimised==optimised_or_not);
        index_of_interest = intersect(model_and_solver,lookup_indices);
        if (~isempty(index_of_interest))
            assert(length(index_of_interest)==1)
            if (s<=2)
                fprintf('&\t%g',10^(-(timesteps(index_of_interest)+2)))
            else
                fprintf('&\t%g',timesteps(index_of_interest))
            end
        else
            fprintf('&\t--')
        end
        
        if (optimised_or_not)
            % Make a fuss if the timestep is different for optimised
            % But only run this analysis for one call of this method
            timesteps_for_with_and_without_lookup_tables = timesteps(model_and_solver);

            fid = fopen('different_steps.txt', 'at');
            if (length(timesteps_for_with_and_without_lookup_tables)==2)            
                if (timesteps_for_with_and_without_lookup_tables(1)~=timesteps_for_with_and_without_lookup_tables(2))
                    fprintf(fid, '%i %i %g %g\n', m, s, timesteps_for_with_and_without_lookup_tables);
                end
            else
                fprintf(fid, '%i %i %g %g\n', m, s, timesteps_for_with_and_without_lookup_tables);
            end
            
        end
    end
    fprintf('\\\\\n')
end

fclose(fid);
end