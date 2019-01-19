function finalplots_variance(f_dout, f_nsubjsubset, f_nsim, f_str_nsubjsubset, f_formatSpec_file, ...
                                         f_nbins, f_param_diff_subjave, f_strparam, f_strparam_short)
    % Plot the histogram curves of the difference between the recovered and fixed value according to subject subsets 
    %      for each parameter (prev & beta) as well as their normalised estimation.
    % Input:  'f_dout'                        folder for data (dout)
    %         'f_nsubjsubset'                 subsets of participant (nsubjsubset)
    %         'f_nsim'                        number of simulations (nsim)
    %         'f_str_nsubjsubset'             string version of nsubjsubset (str_nsubjsubset)
    %         'f_formatSpec_file'             string format for saved figure ('formatSpec_file)
    %         'f_nbins'                       number of bins for raw histogram ('nbins')
    %         'f_param_diff_subjave'          parameter average across subject for each simulation ...
    %                                         ('prev/beta_diff_subjave')
    %         'f_strparam'                    full name of the parameter ...
    %                                         ('probability of reversal' or 'choice variability')        
    %         'f_strparam_short'              shortened name of the parameter ('prev' or 'beta')
    % Output: two figures saved (one for each parameter)


    f_formatSpec_title = '%sifference recovered-fixed %s (2x%d simulations)';

    lgd = cell(length(f_nsubjsubset),1); % to create an 'automatic' legend corresponding with each subset
    for irevindexsubset = 1:length(f_nsubjsubset)
       lgd{irevindexsubset}= sprintf('%d subjects',f_nsubjsubset(irevindexsubset));
    end

    h = figure; % y = number of stimulations
    hold on
    name_title = sprintf(f_formatSpec_title, 'D', f_strparam, f_nsim);
    title(name_title);

    [n, x] = hist(f_param_diff_subjave,f_nbins);  
    plot(x, n, 'LineWidth',2)

    ylabel('Number of simulations');
    xlabel(['Difference recovered-fixed, ',f_strparam]);
    set(gca,'XGrid','on','YGrid','off')
    legend(lgd);
    hold off


    name_fig = sprintf(f_formatSpec_file, f_dout, ['fig_diff_', f_strparam_short, '_recov-fix'], ...
                                          f_nsim, f_str_nsubjsubset, '.png');
    saveas(h, name_fig) 




    clear n y % to avoid dimension issue with just-used variables
    h = figure; % with AUC = 1
    hold on
    name_title = sprintf(f_formatSpec_title, 'Probability density function of the d', f_strparam, f_nsim);
    title(name_title);

    nbins_pdf = length(f_param_diff_subjave);
    [~, x_norm] = hist(f_param_diff_subjave,nbins_pdf);
    y = zeros(nbins_pdf, length(f_nsubjsubset));
    for irevindexsubset = 1:length(f_nsubjsubset)
       pd = fitdist(f_param_diff_subjave(:,irevindexsubset), 'Normal');
       y(:,irevindexsubset) = pdf(pd,x_norm);
    end
    plot(x_norm, y,'LineWidth',2)

    set(gca,'ColorOrderIndex',1) % to reset the default color order
    [n, x] = hist(f_param_diff_subjave,f_nbins);
    empirical_pd = zeros(f_nbins, length(f_nsubjsubset));
    for irevindexsubset = 1:length(f_nsubjsubset)
       empirical_pd(:,irevindexsubset) = n(:,irevindexsubset)/trapz(x,n(:,irevindexsubset)); % so that AUC = 1
    end
    plot(x, empirical_pd,':','LineWidth',1)

    ylabel('Probability density');
    xlabel(['Difference recovered-fixed, ',f_strparam]);
    set(gca,'XGrid','on','YGrid','off')
    legend(lgd);
    hold off


    name_fig = sprintf(f_formatSpec_file, f_dout, ['fig_pdf_diff_', f_strparam_short, '_recov-fix'], ...
                                          f_nsim, f_str_nsubjsubset, '.png');
    saveas(h, name_fig)

end
