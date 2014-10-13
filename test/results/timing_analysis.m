close all 
clear all

build_types = {'IntelProductionCvode',...
               'IntelProduction',...
               'Intel',...
               'GccOpt',...
               'GccOptNative',...
               'debug'};
           
% NB indexed as 0, 1, ..., 8 in the results file
% so we need to add one to those numbers to get the
% correct name!
solvers = {'CVODE AJ', 'CVODE NJ', 'F. Euler', ...         
        'B. Euler','RK2','RK4','Rush Larsen',...            
        'GRL1','GRL2'};
    
print_latex_table = false;
look_at_fake_pde_step_timings = false;

% Compile all the results into a table.
all_results = [];

  
for b=1:length(build_types)
    d = importdata([build_types{b} '_timings.txt']);

    % Get the raw data out
    model = d.textdata;
    solver = d.data(:,1);
    optimised = d.data(:,2);
    times = d.data(:,3);
    
    assert(length(model)==length(solver));
    assert(length(model)==length(optimised));
    assert(length(model)==length(times));
    
    clear d

    if b==1 % This file recently updated to include all the options, this 'if' could go if all updated to the same.
        model_list = unique(model);
        
        % Find out how many ODEs each model has from their summary files
        for i=1:length(model_list)
            d = importdata(['..' filesep 'data/reference_traces/' model_list{i} '.summary']);
            model_list_ODEs(i) = d.data(1);
        end
        
        solver_list = unique(solver);
        optimised_list = unique(optimised);
    end
    
    for model_idx = 1:length(model_list)
        indices_this_model = find(strcmp(model,model_list{model_idx}));
        for solver_idx = 1:length(solver_list)
            indices_this_solver = find(solver==solver_list(solver_idx));
            index_this_combination = intersect(indices_this_model,indices_this_solver);
            for optimised_idx = 1:length(optimised_list)
                indices_this_optimisation = find(optimised==optimised_list(optimised_idx));
                index_complete_combination = intersect(index_this_combination, indices_this_optimisation);
            
                if (~isempty(index_complete_combination))
                    assert(length(index_complete_combination)==1)
                    all_results(model_idx, solver_idx, b, optimised_idx) = times(index_complete_combination);
                else
                    all_results(model_idx, solver_idx, b, optimised_idx) = -1;
                end
            end
        end
    end            
        
    if print_latex_table
        % Now write out a table in latex format

        for optimised_idx = 1:length(optimised_list)
        
            max_time = max(max(all_results(b,:,:,optimised_idx))); 
            colours = colormap(autumn(100));

            % Print another summary table for O'hara colour coded.
            fprintf('%% This table was autogenerated by timing_analysis.m, so don''t edit it here, edit that instead!\n')
            fprintf('%% Please remember to add \\use{multirow} to your document preamble in order to support multirow cells\n');
            fprintf('%% Booktabs require to add \\usepackage{booktabs} to your document preamble.\n')
            fprintf('\\begin{table}[htb]\n')
            fprintf('\\begin{center}\n')
            fprintf('\\textbf{\\refstepcounter{table}\\label{Tab:Timings%s} Table \\arabic{table}.}{ Wall time (s) for 10s of simulation for the %s build }\\\\\n',...
                build_types{b}, put_underscores_in_latex(build_types{b}))
            fprintf('\\vspace{0.2cm}\n')
            fprintf('\\processtable{ }\n')
            fprintf('{\\begin{tabular}{@{}clccccccccccc@{}}\n')
            fprintf('\\toprule\n')
            fprintf('Model & ')
            for solver_idx=1:1:length(solver_list)
                fprintf(' %s ',solvers{solver_list(solver_idx)+1})
                if solver_idx < length(solver_list)
                    fprintf('&')
                end
            end
            fprintf(' \\\\\n \\midrule\n')
            for model_idx=1:length(model_list)
                % Print model name
                fprintf('%s & ',put_underscores_in_latex(model_list{model_idx}))
                % Print timings
                for solver_idx=1:1:length(solver_list)
                    if (all_results(model_idx, solver_idx,b) > 0)
                        % Real result
                        %hts = all_results(b, model_idx, solver_idx)./max_time;
                        %rgb = colours(floor(hts)+1,:);
                        %fprintf('\\cellcolor[rgb]{%f,%f,%f} %3.0f ',rgb(1),rgb(2),rgb(3),hts)
                        fprintf('%4.3f ',all_results(model_idx, solver_idx,b, optimised_idx))
                    else
                        % Wasn't run
                        %fprintf('\\cellcolor[rgb]{1,1,1} -- ')
                        fprintf(' --- ')
                    end            
                    if solver_idx < length(solver_list)
                        fprintf('&')
                    end
                end
                if model_idx < length(model_list)
                    fprintf('\\\\ \\midrule\n')
                end
            end
            fprintf('\\\\ \n\\botrule\n')
            fprintf('\\end{tabular}}{}\n')
            fprintf('\\end{center}\n')
            fprintf('\\end{table}\n')  
        end
    end
end

% Work out ranking of models in terms of number of ODEs
%[~,ordering] = sort(model_list_ODEs);

% % Rank the models in terms of how fast they are
% % in our best case - intel with intel cvode, for cvode numeric J
[~, ordering] = sort(all_results(:, 2, 1, 1));

save('difficulty_ordering.mat','ordering','-mat')

ordered_models = model_list(ordering);
for i=1:length(ordering)
    fprintf('%i Model:%s\n',i,ordered_models{i})
end

analytic_result_rows = find(all_results(ordering, 1, 1, 1)>0);

for s = 1:length(solvers)
    figure
    colorOrder = get(gca, 'ColorOrder');
    subplot(1,2,1)
    for b=1:length(build_types)
        semilogy(all_results(ordering, s, b, 1), '.-','Color',colorOrder(b,:))
        hold all
    end
    xlabel('Model index')
    %xlabel('Model index, ranked in terms of number of ODEs')
    ylabel('Wall time taken to simulate 1 second (s)')
    ylim([1e-4 1e3])
    title(solvers{s})
    legend(build_types,'Location','NorthWest')
    
    % Have a look how much faster the different compilers are
    subplot(1,2,2)
    for b=1:length(build_types)-1
        %plot(all_results(ordering, 2, b, 1),'.-')
        good_indices = find(all_results(ordering, s, 6, 1) > 0); 
        plot(good_indices, all_results(ordering(good_indices), s, 6, 1)./all_results(ordering(good_indices), s, b, 1),'.','MarkerFaceColor',colorOrder(b,:),'MarkerEdgeColor',colorOrder(b,:))
        hold all
        for i=1:length(good_indices)-1
            if good_indices(i)+1 == good_indices(i+1)
                h = plot([good_indices(i) good_indices(i+1)], all_results(ordering(good_indices([i i+1])), s, 6, 1)./all_results(ordering(good_indices([i i+1])), s, b, 1), '-','Color',colorOrder(b,:));
                if (i==1)
                    dont_show_in_legend(h, false);
                else
                    dont_show_in_legend(h, true);
                end
            end
        end    
    end
    title(solvers{s})
    xlabel('Model index')
    xlim([0 length(ordering)+1])
    ylabel('Compiler Speedup relative to gcc debug (x)')
    %legend(build_types(1:5),'Location','EastOutside')
end

figure
color_idx = 1;
legend_idx = 0;
colorOrder = get(gca, 'ColorOrder');

% semilogy(analytic_result_rows,all_results(ordering(analytic_result_rows), 1, 1, 1), '.-')
% legend_idx = legend_idx +1;
% solvers_legend{legend_idx} = solvers{1};
% hold all
% colorOrder = get(gca, 'ColorOrder');
% semilogy(analytic_result_rows,all_results(ordering(analytic_result_rows), 1, 1, 2), ...
%     'Color', colorOrder(color_idx,:), 'Marker','.','LineStyle','--')
% legend_idx = legend_idx +1;
% solvers_legend{legend_idx} = [solvers{1} ' Opt'];
% xlabel('Model indices, ordered by time taken using CVODE NJ')
% ylabel('Wall time taken for 10 paces (s)')

legend_idx = 0;
for i=1:length(solver_list)
    color_idx = color_idx+1;
    if i<8 
        linestyle = '-';
    else
        linestyle = '--';
    end    
    order = find(all_results(ordering, i, 1, 1)>0);  
    if (~isempty(order))
        semilogy(1:length(order), all_results(ordering(order), i, 1, 1), ...
            'Color',colorOrder(mod(color_idx,7)+1,:),...
            'Marker','.','LineStyle','-')
        hold all
        legend_idx = legend_idx + 1;
        solvers_legend{legend_idx} = solvers{solver_list(i)+1};
    end
    order = find(all_results(ordering, i, 1, 2)>0);
    if (~isempty(order))
        semilogy(1:length(order), all_results(ordering(order), i, 1, 2), ...
            'Color',colorOrder(mod(color_idx,7)+1,:),...
            'Marker','.','LineStyle','--')
        hold all
        legend_idx = legend_idx + 1;
        solvers_legend{legend_idx} = [solvers{solver_list(i)+1} ' Opt'];
    end
end
title('Solver benchmarking')
legend(solvers_legend,'Location','EastOutside')
xlabel('Model indices, ordered by time taken using CVODE NJ')
ylabel('Wall time taken to simulate 1 second (s)')

figure
for i=1:length(solver_list)
    valid_results = find(all_results(:, i, 1, 1)>0);
    [~, reordering] = sort(all_results(valid_results, i, 1, 1));
    if i<8 
        linestyle = '-';
    else
        linestyle = '--';
    end
    semilogy(all_results(valid_results(reordering), i, 1,1), ['.' linestyle])
    hold all
end
title('Solver benchmarking')
xlabel('Model indices ordered by time taken for each solver')
ylabel('Wall time taken to simulate 1 second (s)')
legend(solvers{solver_list+1},'Location','EastOutside')

% Have a look how much speed up we get with Lookup tables.
% For this analysis use CVODE AJ,NJ and IntelProductionCvode
figure

for s=1:length(solver_list)
    % The indices in all_results here are model, solver, build, lookuptables (off, on)
    proportions_time = all_results(ordering, s, 1, 2)./all_results(ordering, s, 1, 1);
    good_idx = intersect(find(proportions_time > 0), find(all_results(ordering, s, 1, 1) > 0));
    h = plot([0 length(ordering)+1], [1 1], 'k--');
    dont_show_in_legend(h, true);
    hold on
    
    color_this_plot = colorOrder(mod(s-1,size(colorOrder,1)-1)+1,:);
    
    % Plot all the points
    h = plot(good_idx, 1.0./proportions_time(good_idx), '.', 'MarkerFaceColor',...
        color_this_plot, 'MarkerEdgeColor', color_this_plot, 'MarkerSize', 10.0);
    
    fprintf('Solver %i: median speedup w. lookup tables = %g\n',s-1,median(1.0./proportions_time(good_idx)))
    
    if s<=6 
        linestyle = '.-';
    else
        linestyle = '.--';
    end

    % Uncomment these to get the 'lines' versions...
    dont_show_in_legend(h, true);
    for i=1:length(good_idx)-1
        % If there's no model missing, join the dots
        if good_idx(i)+1 == good_idx(i+1)
            h = plot([good_idx(i) good_idx(i+1)], 1.0./proportions_time(good_idx([i i+1])), linestyle, 'Color', color_this_plot);
            if (i==1)
                dont_show_in_legend(h, false);
            else
                dont_show_in_legend(h, true);
            end
        end
    end
end
xlabel('Model indices, ordered by time taken using CVODE NJ, no lookup tables')
xlim([0 length(ordering)+1])
title('Relative speed using Lookup Tables')
ylabel('Speed with Lookup Tables relative to without (x)')
legend(solvers,'Location','SouthWest')

% Some asserts in here, if they fail, some text in paper will need
% updating!

% Just look to see how many failures there are... for CVODE AJ without lookup
tmp1 = find(all_results(ordering, 1, 1, 1) < 0);
assert(length(tmp1) == 7);

% and with lookup tables
tmp2 = find(all_results(ordering, 1, 1, 2) < 0);
assert(length(tmp2) == 9);
% And check that all the ones that fell over (or weren't run) without lookup tables, 
% still don't work with lookup tables...
assert(length(intersect(tmp1,tmp2))==length(tmp1));

% Just look to see how many failures there are... for CVODE NJ without lookup
tmp = find(all_results(ordering, 2, 1, 1) < 0);
assert(isempty(tmp)) % No failures

% and with lookup tables
tmp = find(all_results(ordering, 2, 1, 2) < 0);
assert(length(tmp) == 6);

% analytic_result_rows_opt = find(all_results(ordering, 1, 1, 2)>0);
% figure
% semilogy(analytic_result_rows_opt,all_results(ordering(analytic_result_rows_opt), 1, 1, 2), '.-')
% xlabel('Model indices ordered by time taken for each Opt solver')
% ylabel('Wall time taken for 10 paces (s)')
% hold all
% for i=2:8
%     [~, ordering] = sort(all_results(:, i, 1, 2));
%     if i<8 
%         linestyle = '-';
%     else
%         linestyle = '--';
%     end
%     semilogy(all_results(ordering, i, 1, 2), ['.' linestyle])
% end
% title('Solver benchmarking')
% legend(solvers{solver_list+1},'Location','EastOutside')
% xlim([1 67]) % Include Clancy-Rudy again.



