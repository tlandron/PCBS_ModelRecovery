function finalplots_variance(f_dout, f_nsubjsubset, f_nsim, f_str_nsubjsubset, f_formatSpec_file, f_param_diff_subjave, f_strparam, f_strparam_short)




f_formatSpec_title = 'Difference recovered-fixed %s (2x%d simulations)';

lgd = cell(length(f_nsubjsubset),1);
for irevindexsubset = 1:length(f_nsubjsubset)
   lgd{irevindexsubset}= sprintf('%d subjects',f_nsubjsubset(irevindexsubset));
end





h = figure; % y = number of stimulations
hold on
name_title = sprintf(f_formatSpec_title, f_strparam, f_nsim);
title(name_title);

nbins = 50;
[n, x] = hist(f_param_diff_subjave,nbins);  
plot(x, n, 'LineWidth',2)

ylabel('Number of simulations');
xlabel(['Difference recovered-fixed, ',f_strparam]);
legend(lgd);
hold off


name_fig = sprintf(f_formatSpec_file, f_dout,['diff_',f_strparam_short,'_recov-fix'],f_nsim,f_str_nsubjsubset,'.png');
saveas(h, name_fig) 




clear n y
h = figure; % with normalised distribution, AUC = 1
hold on
name_title = sprintf(f_formatSpec_title, f_strparam, f_nsim);
title(name_title);

nbins = 1000;
[~, x_norm] = hist(f_param_diff_subjave,nbins);
y = zeros(nbins, length(f_nsubjsubset));
for irevindexsubset = 1:length(f_nsubjsubset)
   pd = fitdist(f_param_diff_subjave(:,irevindexsubset), 'Normal');
   y(:,irevindexsubset) = pdf(pd,x_norm);
end
plot(x_norm, y,'LineWidth',2)

set(gca,'ColorOrderIndex',1) % to reset the default color order
nbins = 50;
[n, x] = hist(f_param_diff_subjave,nbins);
empirical_pd = zeros(nbins, length(f_nsubjsubset));
for irevindexsubset = 1:length(f_nsubjsubset)
   empirical_pd(:,irevindexsubset) = n(:,irevindexsubset)/trapz(x,n(:,irevindexsubset)); % so that AUC = 1
end
plot(x, empirical_pd,':','LineWidth',1)

ylabel('???');
xlabel(['Difference recovered-fixed, ',f_strparam]);
legend(lgd);
hold off


name_fig = sprintf(f_formatSpec_file, f_dout,['diff_',f_strparam_short,'_recov-fix'],f_nsim,f_str_nsubjsubset,'.png');
saveas(h, name_fig)

end
