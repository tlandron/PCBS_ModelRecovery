function finalplots_mindiff (f_dout, f_nsubjsubset, f_nsim, f_nttest, f_ndiff, f_ndiffdim, f_param_diff_index, ...
                             f_param_diffsigni_acrossdiff, f_strparam, f_strparam_short, f_strparam_cte)


    f_formatSpec_title = ['Significant t-tests as a function of the difference tested for %s ', 10, ...
                        '(%s constant, for %d t-tests on %d sim each, %d subjects)'];
    f_formatSpec_fig = '%s%s_diffsigni_acrossdiff_%dsubj_%dx%dsim';

   
    f_paramonly_diffsigni_acrossdiff = f_param_diffsigni_acrossdiff(:,:,:,f_param_diff_index); % to only select the differences for either prev or beta

    y = zeros(length(f_param_diff_index), 1); % preallocation
    for revindexsubset = 1:length(f_nsubjsubset) % loop to obtain the percentage of significant ttests
        for i = 1:length(f_param_diff_index)
            div = f_paramonly_diffsigni_acrossdiff(:,:,revindexsubset,i);
            y(i) = 100 * div(:,1) ./ div(:,2);
        end
        x = f_ndiff(f_param_diff_index,f_ndiffdim);

        h = figure;
        hold on
        name_title = sprintf(f_formatSpec_title, f_strparam,f_strparam_cte,f_nttest,f_nsim,f_nsubjsubset(revindexsubset));
        title(name_title);
        plot(x,y,'o-','LineWidth',2)
        ylabel('Percentage of significant t-tests');
        hold off

        name_fig = sprintf(f_formatSpec_fig, f_dout,f_strparam_short,f_nsubjsubset(revindexsubset),f_nttest,f_nsim,'.png');
        saveas(h, name_fig);   
    end
end

