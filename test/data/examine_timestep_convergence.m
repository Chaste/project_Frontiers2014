% You can run this after "CalculateTimesteps" to look at the resulting
% action potential traces relative to the definitive reference traces.

close all
clear all

file_listing = dir(['reference_traces' filesep '*.dat']);

for i=1:length(file_listing)
    
    % Skip the files that are to do with tissue sims
    if (~isempty(strfind(file_listing(i).name, 'tissue')))
        continue
    end
    
    file_listing(i).name
    
    d = importdata(['reference_traces' filesep file_listing(i).name]);
    if (isstruct(d))
        data = d.data;
    else
        data = d;
    end

    figure(i)
    plot(data(:,1), data(:,2),'b-')
    hold on
    title(file_listing(i).name)    
    xlabel('Time (ms)')
    xlabel('Voltage (mV)')   
end

% 
% 
% 
% file_listing = dir(['reference_traces' filesep '*_GccOpt_tissue_*.dat']);
% 
% % Pull out the model names
% for i=1:length(file_listing)
%     file = file_listing(i).name;
%     k = strfind(file, '_GccOpt_tissue_');
%     model{i} = file(1:k-1);
% end
% 
% model_list = unique(model);
% 
% % Let's look at summary files first for each model
% for m=1:length(model_list)
%     figure
%     fprintf('Processing model: %s\n', model_list{m});
%     for compiler=1:3
%     
%         if compiler==1 % IntelProduction on travis without CVODE reset
%             prefix = '_tissue_';
%             suffix = '';
%         elseif compiler==2 % GccOptNative on skip without CVODE reset
%             prefix = '_GccOpt_tissue_';
%             suffix = '';
%         elseif compiler==3 % IntelProduction on travis with CVODE reset
%             prefix = '_tissue_';
%             suffix = '_Reset';
%         end
%         file_listing = dir(['reference_traces' filesep model_list{m} prefix '*' suffix '.summary']);
% 
%         % Look at the two mesh resolutions separately.    
%         for h = [0.01 0.001]
%             pde_timesteps = [];
%             % Find all the traces to do with this h
%             for i=1:length(file_listing)
%                 file = file_listing(i).name;
%                 j = strfind(file, '_pde_');
%                 k = strfind(file, ['_h_' num2str(h)]);
%                 if (isempty(k))
%                     continue
%                 end
%                 pde_timesteps = [pde_timesteps; str2num(file(j+5:k-1))];
%             end
% 
%             % Right, now we have all the available combinations, let's plot the
%             % convergence of the summary metrics
%             for pde_idx = 1:length(pde_timesteps)
%                 pde_step = pde_timesteps(pde_idx);
%                 file_of_interest = ['reference_traces' filesep model_list{m} ...
%                     prefix 'pde_' num2str(pde_step)  '_h_' num2str(h) suffix '.summary'];
%                 d = importdata(file_of_interest);
%                 metrics(:,pde_idx) = d.data;
%             end
% 
%             num_metrics = size(metrics,1);
% 
%             for i=1:num_metrics
%                 assert(num_metrics==6)
%                 subplot(2,3,i)
%                 if compiler == 1
%                     linestyle = '.-';
%                 elseif compiler == 2
%                     linestyle = '.--';
%                 elseif compiler == 3
%                     linestyle = 'o-';
%                 end
%                 semilogx(pde_timesteps,metrics(i,:),linestyle)
%                 hold all
%                 if i == 2
%                     title(strrep(model_list{m}, '_', ' '))
%                 end
%                 xlim([1e-3 1])
%                 xlabel('PDE Timestep (ms)')
%                 set(gca,'XTick',[0.001 0.01 0.1 1])
%                 ylabel(strrep(d.textdata{i}, '_', ' '))
%             end        
%             clear metrics
%         end
%     end
%     legend('h = 0.01 cm Intel no Reset', 'h = 0.001 cm Intel no Reset',...
%            'h = 0.01 cm GccOpt no Reset', 'h = 0.001 cm GccOpt no Reset',...
%            'h = 0.01 cm Intel Reset', 'h = 0.001 cm Intel Reset', 'Location', 'EastOutside')
% end
