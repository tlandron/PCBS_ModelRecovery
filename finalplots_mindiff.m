function finalplots_mindiff (f_dout, f_nsubjsubset, f_nsim, f_nttest, f_ndiff, f_ndiffdim,      ...
                             f_diffparam_fit, f_param_diff_index, f_param_diffsigni_acrossdiff, ... 
                             f_strparam, f_strparam_short, f_strparam_cte)
    % Plots the number of significant ttests according to the difference tested.
    % Input:  'f_dout'                        folder for data (dout)
    %         'f_nsubjsubset'                 subsets of participant (nsubjsubset)
    %         'f_nsim'                        number of simulations (nsim)
    %         'f_nttest'                      number of ttests (nttest)
    %         'f_ndiff'                       differences to be tested (ndiff)
    %         'f_ndiffdim'                    1 for prev, 2 for beta
    %         'f_diffparam_fit'               difference between the twotasks in the actual dataset ...
    %                                         (prev/beta_diff_tasks)
    %         'f_param_diff_index'            index of the revelant difference for each ...
    %                                         parameter (prev/beta_diff_index)
    %         'f_param_diffsigni_acrossdiff'  results of ttests for the two parameters for all ...
    %                                         differences to be tested (prev/beta_diffsigni_acrossdiff)
    %         'f_strparam'                    full name of the parameter ...
    %                                         ('probability of reversal' or 'choice variability')        
    %         'f_strparam_short'              shortened name of the parameter ('prev' or 'beta')
    %         'f_strparam_cte'                constant parameter ('beta' or 'prev')
    % Output: two figures saved (one for each parameter)

    f_formatSpec_title = ['Significant t-tests as a function of the difference tested for %s ', 10, ...
                        '(%s constant, for %d t-tests on %d sim each, %d subjects)']; % ',10,' for new line
    f_formatSpec_fig = '%sfig_%s_diffsigni_across%ddiff_%dsubj_%dx%dsim%s';

   
    f_paramonly_diffsigni_acrossdiff = f_param_diffsigni_acrossdiff(:, :, :, f_param_diff_index); % to only select the differences 
                                                                                                  % for either prev or beta

    y = zeros(length(f_param_diff_index), 1); % preallocation of percentage of significant ttests
    for revindexsubset = 1:length(f_nsubjsubset) % loop to obtain the percentage of significant ttests
        for i = 1:length(f_param_diff_index)
            div = f_paramonly_diffsigni_acrossdiff(:, :, revindexsubset, i);
            y(i) = 100 * div(:, 1) ./ div(:, 2);
        end
        x = f_ndiff(f_param_diff_index, f_ndiffdim);

        h = figure;
        hold on
        name_title = sprintf(f_formatSpec_title, f_strparam, f_strparam_cte, f_nttest, f_nsim, f_nsubjsubset(revindexsubset));
        title(name_title);
        plot(x, y, 'o-', 'LineWidth', 2)
        
        xcst = [f_diffparam_fit f_diffparam_fit];             % to draw a vertical line for the difference ...
        ylimline = get(gca,'ylim');                           % between the two tasks in the actual dataset    
        line(xcst, ylimline, 'Color','red','LineStyle','--'); % also works with xline function in MATLAB 2018

        set(gca, 'XGrid', 'on', 'YGrid', 'off')
        ylabel('Percentage of significant t-tests');
        xlabel('Difference tested')
        hold off

        name_fig = sprintf(f_formatSpec_fig, f_dout, f_strparam_short, length(f_param_diff_index), ...
                           f_nsubjsubset(revindexsubset), f_nttest, f_nsim,'.png');
        saveas(h, name_fig);   
    end
end

