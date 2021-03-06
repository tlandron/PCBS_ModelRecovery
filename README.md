# Development of a model recovery procedure for the study of reversal learning
## Validating project for PCBS course (Cogmaster, 2018-2019).
[github page](https://tlandron.github.io/PCBS_ModelRecovery/) - [github repository](https://github.com/tlandron/PCBS_ModelRecovery)

The original context for this recovery procedure can be in found [here](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/PCBS_ModelRecovery_Plan_TLANDRON.pdf).

Briefly, the goal of this project is to recover a model and its fitting procedure for a 2-task reversal learning paradigm. First, the recovered parameters of interest (the probability of reversal 'prev' and choice variability 'beta') will be studied across runs of simulations & fits, especially focusing on the mean difference between the recovered parameter and a reference value (from the actual dataset) to reveal any under- or overestimation.
Second, the sensitivity of the fitting procedure will be studied by testing different values for each task within parameters.


CAUTION: the web version of git was used almost until the very end of this project, resulting in all the files being uploaded in the main directory. However, for the scripts to run correctly, the initial folder organisation has to be respected (as well as having the data). The whole 'PBCS_ModelRecovery' was therefore commited and pushed using the terminal, the 'initially commited files' (i.e., those initially commited in the main directory) were nevertheless kept in the main directory so as to keep their commit history -- as asked -- the (script and figure) links in the following document refer to those 'initially commited files'. 

/!\ To test the scripts, clone the whole repository and run the PCBS_ModelRecovery/**Code/**..'script'.mlx (not the PCBS_ModelRecovery/..'script'.mlx)

Also, due to figure issues, reading README.md from the git repository is strongly advised.

---

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [Given scripts](#given-scripts)
- [Initial fits](#initial-fits)
- Part 1: Study of the difference between the parameter recovered and reference values
    - [Main script (part 1)](#main-script-part-1)
    - [Final plots script (part 1)](#final-plots-script-part-1)
    - [Conclusion (part 1)](#conclusion-part-1)
- Part 2 : Study of the sensitivity of the fitting procedure
    - [Main script (part 2)](#main-script-part-2)
    - [Final plots script (part 2)](#final-plots-script-part-2)
    - [Conclusion (part 2)](#conclusion-part-2)
- [General conclusion](#general-conclusion)
- [License](#license)

<!-- markdown-toc end -->

---


## Given scripts
- to load the data [loading script](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/loaddata.m):
    - creates a massive data array 'dat' storing the data of all subjects and conditions;
- simulation procedure [simulation script](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/sim_model_softmax.m):
    - produces simulated datasets from actual data, stores them in data array (cfg);
- fit procedure [fit script](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/fit_model_softmax.m):
    - given datasets (cfg), computes the two parameters of interest: the probability of reversal (prev) and choice variability (beta) for each subject and condition.

## Initial fits
A [script to run initial fits](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/initialfit.m), in order to obtain a reference value fixed for each variable of interest (prev_fit, beta_fit) was a adapted from the [fit script](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/fit_model_softmax.m) as a function that runs for 'nrep' repetitions of the fitting procedure using the actual dataset.


    function [outf_prev_fit, outf_beta_fit] = initialfit(f_nrep, f_nsubj, f_dat, f_subjlist)
        % Run 'nrep' repetitions of the fit procedure.
        % Input:  'f_nrep'      the number of repetitions
        %         'f_nsubj'     the number of subjects (integer output from loaddata.m)
        %         'f_dat'       the data (structure output from loaddata.m)
        %         'f_subjlist'  the list of subjects (array output from loaddata.m)
        % Output: 'outf_prev_fit outf_beta_fit' two 3D matrix,
        %             One matrix/parameter, structured as follow: 'outf_param_fit(isubj, itask, irep)', 
        %                 first dimension 'isubj' corresponds to each subject, 
        %                 second dimension 'itask' corresponds to each task,
        %                 third dimension 'irep' corresponds to each repetition.
    
        outf_prev_fit = zeros(f_nsubj, 2, f_nrep);
        outf_beta_fit = zeros(f_nsubj, 2, f_nrep);

        % given script -->
        for irep = 1:f_nrep
            for isubj = 1:f_nsubj
                for itask = 1:2

                    % filter trials of interest
                    ifilt = f_dat.subj == f_subjlist(isubj) & f_dat.task == itask;

                    % fit model
                    cfg_fit        = [];
                    cfg_fit.seqllr = f_dat.evid(ifilt);
                    cfg_fit.seqind = f_dat.sind(ifilt);
                    cfg_fit.rbef   = f_dat.rbef(ifilt);
                    cfg_fit.raft   = f_dat.raft(ifilt);
                    cfg_fit.brep   = 0; % force no repetition bias
                    cfg_fit.epsi   = 0; % force no lapses

                    out_fit = fit_model_softmax(cfg_fit);

                    % store model parameters
                    outf_prev_fit(isubj, itask, irep) = out_fit.prev;
                    outf_beta_fit(isubj, itask, irep) = out_fit.beta;

                end
            end
        end
        % <-- given script
    end




## Part 1: Study of the difference between the parameter recovered and reference values
### [Main script (part 1)](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/script_recovery_ACTOBS_GTS_TL.mlx)
Data loading.
    
    loaddata


Inital fits (nrep number of fits) of the actual data to obtain an estimate of the two variables of interest, the probability of reversal and the choice probability, that will serve as fixed refence values in the next analyses (prev_fix, beta_fix). 

    nrep = 10;

    clear prev_fit beta_fit % to prevent any dimension error when using MATLAB LiveScript environment 
    [prev_fit, beta_fit] = initialfit(nrep, nsubj, dat, subjlist);

    prev_fix = mean(prev_fit(:)) % same fixed value across the tasks     % prev_fix = 0.2087
    beta_fix = mean(beta_fit(:))                                         % beta_fix = 1.2439


In this first part, the two tasks are not differentiated: the parameters are the same across the 2 tasks, resulting in 2x'nsim' simulated datasets.

    nsim              = 5000;         % number of simulations
    nsubjsubset       = [12 24 48 96] % subset(s) of participants (w/out consideration for condition, 96 = all)

    % preallocation, to save data for each subset
    prev_recov_subjave = zeros(2*nsim, length(nsubjsubset)); % recovered value of each parameter averaged ...
    beta_recov_subjave = zeros(2*nsim, length(nsubjsubset)); % across subjects, for each subset

    prev_diff_subjave = zeros(2*nsim, length(nsubjsubset)); % difference (recovered-fixed value) of each parameter ...
    beta_diff_subjave = zeros(2*nsim, length(nsubjsubset)); % averaged across subjects, for each subset

    meansd_prev_diff_subjave = zeros(length(nsubjsubset), 2); % mean & sd of the difference, for each subset
    meansd_beta_diff_subjave = zeros(length(nsubjsubset), 2);

    mean_corr_prev_beta = zeros(length(nsubjsubset), 1); % mean correlation coefficient between prev and beta, for each subset
    shared_var          = zeros(length(nsubjsubset), 1); % shared variance between prev and beta, for each subset


    for isubjsubset = nsubjsubset % for each subset of subjects

        prev_recov = zeros(isubjsubset, 2, nsim); % to save the parameters for each subject
        beta_recov = zeros(isubjsubset, 2, nsim);

        isubjsubset % to display the current subset running in MATLAB LiveScript environment 
        seq_isubjsubsest = randperm(nsubj, isubjsubset) % to avoid selecting deterministically the same subjects across subsets
        for isubj = seq_isubjsubsest
        % given script -->
            for itask = 1:2

                % filter trials of interest
                ifilt = dat.subj == subjlist(isubj) & dat.task == itask;

                % simulate model
                cfg_sim        = [];
                cfg_sim.nsim   = nsim;
                cfg_sim.seqllr = dat.evid(ifilt);
                cfg_sim.seqind = dat.sind(ifilt);
                cfg_sim.prev   = prev_fix;
                cfg_sim.beta   = beta_fix;
                cfg_sim.brep   = 0; % force no repetition bias
                cfg_sim.epsi   = 0; % force no lapses

                out_sim        = sim_model_softmax(cfg_sim);

   Then each subject is fit using the function "fit_model_softmax".
   
                % fit simulated decisions
                for isim = 1:nsim
                    cfg_fit        = [];
                    cfg_fit.seqllr = out_sim.cfg.seqllr;
                    cfg_fit.seqind = out_sim.cfg.seqind;
                    cfg_fit.raft   = out_sim.raft(:, isim);
                    cfg_fit.rbef   = out_sim.rbef(:, isim);
                    cfg_fit.brep   = out_sim.cfg.brep; % force no repetition bias
                    cfg_fit.epsi   = out_sim.cfg.epsi; % force no lapses

                    out_recov = fit_model_softmax(cfg_fit);
        % <-- given script
        
   Estimated values of prev and beta from each simulation are stored.

                    % store model parameters
                    revindexsubj = find(seq_isubjsubsest == isubj); % isubj goes over the index limit due to the 
                                                                    % use of subsets --> reverse indexing
                    prev_recov(revindexsubj, itask, isim) = out_recov.prev;
                    beta_recov(revindexsubj, itask, isim) = out_recov.beta;
                end    
            end
        end
        
Variables to plot histograms + histograms of the recovered parameters 
(if figures wanted, uncomment;  CAUTION: variables stored across subsets needed for final saving).

        prev_recov_cattask = squeeze(cat(3, prev_recov(:, 1, :), prev_recov(:, 2, :)))'; % concatenation of the data for the two ...
        beta_recov_cattask = squeeze(cat(3, beta_recov(:, 1, :), beta_recov(:, 2, :)))'; % tasks into one array (as if only one task)
                                                                                         % transposed matrix for dimension reason
        revindexsubset = find(nsubjsubset == isubjsubset); 
        prev_recov_subjave(:, revindexsubset) = mean(prev_recov_cattask, 2);
        beta_recov_subjave(:, revindexsubset) = mean(beta_recov_cattask, 2);

    %     % figure to be played in MATLAB LiveScript environment, but not saved.
    %     figure; 
    %     hold on
    %     title(["Probability of reversal (from n = ", num2str(isubjsubset), " subjects)"]);
    %     hist(prev_recov_subjave(:, revindexsubset));
    %     hold off
    %     
    %     figure;
    %     hold on
    %     title(["Choice probability (from n = ", num2str(isubjsubset), " subjects)"]);
    %     hist(beta_recov_subjave(:, revindexsubset));
    %     hold off

Variables to plot histograms + histograms of the difference between the recovered-fixed for each parameter 
(if figures wanted, uncomment; CAUTION: variables stored across subsets needed for final saving & final plots).

        prev_diff = prev_recov - prev_fix;
        beta_diff = beta_recov - beta_fix;

        prev_diff_cattask = squeeze(cat(3, prev_diff(:, 1, :), prev_diff(:, 2, :)))'; % concatenation of the data for the two ...
        beta_diff_cattask = squeeze(cat(3, beta_diff(:, 1, :), beta_diff(:, 2, :)))'; % tasks into one array (as if only one task)
                                                                                      % transposed matrix for dimension reason

        prev_diff_subjave(:, revindexsubset) = mean(prev_diff_cattask, 2);
        beta_diff_subjave(:, revindexsubset) = mean(beta_diff_cattask, 2);

    %     figure;
    %     hold on
    %     title(["Difference recovered-fixed probability of reversal (from n = ", num2str(isubjsubset), " subjects)"]);
    %     hist(prev_diff_subjave(:, revindexsubset));
    %     hold off
    %     
    %     figure;
    %     hold on
    %     title(["Difference recovered-fixed choice probability (from n = ", num2str(isubjsubset), " subjects)"]);
    %     hist(beta_diff_subjave(:, revindexsubset));
    %     hold off

Mean & variance of the distribution of the difference recovered-fixed for each parameter across simulations and for each subset of subjects.
        
        meansd_prev_diff_subjave(revindexsubset, 1) = mean(prev_diff_subjave(:, revindexsubset));
        meansd_prev_diff_subjave(revindexsubset, 2) = std(prev_diff_subjave(:, revindexsubset));

        meansd_beta_diff_subjave(revindexsubset, 1) = mean(beta_diff_subjave(:, revindexsubset));
        meansd_beta_diff_subjave(revindexsubset, 2) = std(beta_diff_subjave(:, revindexsubset));

Correlation coefficient and shared variance between the two parameters between subject is computed for each simulation.  
        
        corr_prev_beta = zeros(2*nsim, 1);
        for isim = 1:2*nsim
            corr_prev_beta(isim) = corr(prev_recov_cattask(isim,: )', beta_recov_cattask(isim, :)'); % transposed to use function 'corr'
        end

        mean_corr_prev_beta(revindexsubset) = mean(corr_prev_beta)
        shared_var(revindexsubset)          = mean_corr_prev_beta(revindexsubset).^2;

    end



Saving workspace. 
NB: As git only allows files < 25 Mb, the commented code was used to split the variable of interest into two separate files ([prev_variables](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/prev_recov_2x5000sim_12-24-48-96subj.mat) & [beta_variables](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/beta_recov_2x5000sim_12-24-48-96subj.mat)).
    
    dout = '../Out/'; % folder for produced data

    str_nsubjsubset = ''; % to nicely convert nsubjsubset into a string, e.g., "12-24-48-96"
    for isubjsubset = nsubjsubset
        str_nsubjsubset = [str_nsubjsubset, num2str(isubjsubset), '-'];
    end
    str_nsubjsubset = str_nsubjsubset(1:end-1) % to delete the last useless '-'

    formatSpec_file = '%s%s_%dsim_%ssubj%s';

    name_file = sprintf(formatSpec_file, dout, 'recov', nsim, str_nsubjsubset, '.mat');
    save(name_file) % to save the whole workspace
    
    % % Or use the following to save only specific variables:
    % name_file = sprintf(formatSpec_file, dout, 'prev_recov', nsim, str_nsubjsubset, '.mat');
    % save(name_file, 'nsim', 'nsubjsubset', 'prev_fix', 'prev_recov_subjave', 'prev_diff_subjave',    ...
    %                 'meansd_prev_diff_subjave', 'mean_corr_prev_beta', 'shared_var')
    % 
    % name_file = sprintf(formatSpec_file, dout, 'beta_recov', nsim, str_nsubjsubset, '.mat');
    % save(name_file, 'nsim', 'nsubjsubset', 'beta_fix', 'beta_recov_subjave', 'beta_diff_subjave',    ...
    %                 'meansd_beta_diff_subjave', 'mean_corr_prev_beta', 'shared_var')


Histograms of the paramterers according the subset they were estimated from (in this case, n = 12, 24, 48 or 96).

    nbins = 100; % can be manually changed to the 'apparent noise' on the raw curves
            
    finalplots_variance(dout, nsubjsubset, nsim, str_nsubjsubset, formatSpec_file, ...
                              nbins, prev_diff_subjave, 'probability of reversal', 'prev')

    finalplots_variance(dout, nsubjsubset, nsim, str_nsubjsubset, formatSpec_file, ...
                              nbins, beta_diff_subjave, 'choice variability', 'beta')
                      
### [Final plots script (part 1)](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/finalplots_variance.m)

    function finalplots_variance(f_dout, f_nsubjsubset, f_nsim, f_str_nsubjsubset, f_formatSpec_file, ...
                                         f_nbins, f_param_diff_subjave, f_strparam, f_strparam_short)
        % Plot the histogram curves of the difference between the recovered and fixed value according to subject subsets 
        %      for each parameter (prev & beta) as well as their normalised estimation.
        % Input:  'f_dout'                        folder for data (dout)
        %         'f_nsubjsubset'                 subsets of participant (nsubjsubset)
        %         'f_nsim'                        number of simulations (nsim)
        %         'f_str_nsubjsubset'             string version of nsubjsubset (str_nsubjsubset)
        %         'f_formatSpec_file'             string format for saved figure (formatSpec_file)
        %         'f_nbins'                       number of bins for raw histogram (nbins)
        %         'f_param_diff_subjave'          parameter average across subject for each simulation ...
        %                                         (prev/beta_diff_subjave)
        %         'f_strparam'                    full name of the parameter ...
        %                                         ('probability of reversal' or 'choice variability')        
        %         'f_strparam_short'              shortened name of the parameter ('prev' or 'beta')
        % Output: two figures saved (one for each parameter)


        f_formatSpec_title = '%sifference recovered-fixed %s (2x%d simulations)';

        lgd = cell(length(f_nsubjsubset), 1); % to create an 'automatic' legend corresponding with each subset
        for irevindexsubset = 1:length(f_nsubjsubset)
           lgd{irevindexsubset} = sprintf('%d subjects', f_nsubjsubset(irevindexsubset));
        end


        h = figure; % y = number of stimulations
        hold on
        name_title = sprintf(f_formatSpec_title, 'D', f_strparam, f_nsim);
        title(name_title);

        [n, x] = hist(f_param_diff_subjave, f_nbins);  
        plot(x, n, 'LineWidth', 2)

        ylabel('Number of simulations');
        xlabel(['Difference recovered-fixed, ', f_strparam]);
        set(gca, 'XGrid', 'on', 'YGrid', 'off')
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
        [~, x_norm] = hist(f_param_diff_subjave, nbins_pdf);
        y = zeros(nbins_pdf, length(f_nsubjsubset));
        for irevindexsubset = 1:length(f_nsubjsubset)
           pd = fitdist(f_param_diff_subjave(:, irevindexsubset), 'Normal');
           y(:, irevindexsubset) = pdf(pd, x_norm);
        end
        plot(x_norm, y, 'LineWidth', 2)

        set(gca, 'ColorOrderIndex', 1) % to reset the default color order
        [n, x] = hist(f_param_diff_subjave, f_nbins);
        empirical_pd = zeros(f_nbins, length(f_nsubjsubset));
        for irevindexsubset = 1:length(f_nsubjsubset)
           empirical_pd(:, irevindexsubset) = n(:, irevindexsubset)/trapz(x, n(:, irevindexsubset)); % so that AUC = 1
        end
        plot(x, empirical_pd, ':', 'LineWidth', 1)

        ylabel('Probability density');
        xlabel(['Difference recovered-fixed, ', f_strparam]);
        set(gca, 'XGrid', 'on', 'YGrid', 'off')
        legend(lgd);
        hold off

        name_fig = sprintf(f_formatSpec_file, f_dout, ['fig_pdf_diff_', f_strparam_short, '_recov-fix'], ...
                                              f_nsim, f_str_nsubjsubset, '.png');
        saveas(h, name_fig)

    end

### Conclusion (part 1)

<img src="./fig_diff_prev_recov-fix_2x5000sim_12-24-48-96subj.png" width="415"> <img src="./fig_diff_beta_recov-fix_2x5000sim_12-24-48-96subj.png" width="415">

<img src="./fig_pdf_diff_prev_recov-fix_2x5000sim_12-24-48-96subj.png" width="415"> <img src=".//fig_pdf_diff_beta_recov-fix_2x5000sim_12-24-48-96subj.png" width="415">

The plots show a slight surestimation in the recovered values (e.g., around +0.0050 for the probability of reversal and around +0.0526 for choice variability, both for 24 subjects -- can be extracted from 'meansd_prev_diff_subjave' and 'meansd_beta_diff_subjave'). Unsurprisingly, the more subjects, the smaller the variance. 

From 'shared_var', the shared variance between the two parameter is extracted to suggest 12-13% of shared variance between the probability of reversal and choice variability, whatever the subject subset size. 


 
## Part 2 : Study of the sensitivity of the fitting procedure
### [Main script (part 2)](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/script_recovery_ACTOBS_GTS_TL(mindiff).mlx)

Data loading.

    loaddata

    dout = '../Out/'; % folder for produced data
    
Inital fits (nrep number of fits) of the actual data to obtain estimates of the two variables of interest, the probability of reversal and the choice probability.

    nrep = 10;

    clear prev_fit beta_fit % to prevent any dimension error when using MATLAB LiveScript environment 
    [prev_fit, beta_fit] = initialfit(nrep, nsubj, dat, subjlist);

    prev_fit_tasks = mean(mean(prev_fit, 1),3) % to display prev mean for each task     % [0.2411 0.1763]
    beta_fit_tasks = mean(mean(beta_fit, 1),3) % to display beta mean for each task     % [1.3808 1.1071]

    mean_prev_fit = mean(prev_fit(:))
    mean_beta_fit = mean(beta_fit(:))


Using a function "sim_model_softmax", datasets are simulated using different values as fixed parameters, in order to study the ability of the fitting procedure to recover a certain different within a parameter while the other is kept constant (see 'ndiff'; nsim simulated datasets).

    nsim              = 30;                % number of simulations
    nttest            = 100;               % number of ttests
    ndiff             = [0      0   ; ...  % couples of differences to be tested
                         0.01   0   ; ...
                         0.02   0   ; ...
                         0.03   0   ; ...
                         0.04   0   ; ...
                         0.05   0   ; ...
                         0.1    0   ; ...
                         0      0.1 ; ...
                         0      0.2 ; ...
                         0      0.3 ; ...
                         0      0.4 ; ...
                         0      0.5 ; ...
                         0      1   ]     
    nsubjsubset       = [24]              % subset(s) of participants (w/out consideration for medical condition, 96 = all)

    % preallocation
    prev_diffsigni_acrossdiff = zeros(1, 2, length(nsubjsubset), 1); % significant ttest counter across the differences tested
    beta_diffsigni_acrossdiff = zeros(1, 2, length(nsubjsubset), 1);

    revindexdiff = 0
    for idiff = ndiff'
        revindexdiff = revindexdiff + 1; % to facilitate counting of the difference outside the loop

        idiff % to display the difference being computed in MATLAB LiveScript environment 
        prev_fix = [mean_prev_fit-(idiff(1)/2) mean_prev_fit+(idiff(1)/2)]
        beta_fix = [mean_beta_fit-(idiff(2)/2) mean_beta_fit+(idiff(2)/2)]

        % preallocation
        prev_ttest = zeros(nttest, 5, length(nsubjsubset)); % ttest results saved
        beta_ttest = zeros(nttest, 5, length(nsubjsubset));

        prev_diffsigni = zeros(1, 2, length(nsubjsubset)); % significant ttest counter
        beta_diffsigni = zeros(1, 2, length(nsubjsubset));

        for isubjsubset = nsubjsubset
                % preallocation, to save the parameters for each subject
                prev_recov = zeros(isubjsubset, 2, nsim);
                beta_recov = zeros(isubjsubset, 2, nsim);

                isubjsubset % to display the current subset running in MATLAB LiveScript environment 
                seq_isubjsubsest = randperm(nsubj, isubjsubset) % to avoid selecting deterministically the same subjects across subsets
            for ittest = 1:nttest
                for isubj = seq_isubjsubsest
                % given script -->
                    for itask = 1:2

                        % filter trials of interest
                        ifilt = dat.subj == subjlist(isubj) & dat.task == itask;

                        % simulate model
                        cfg_sim        = [];
                        cfg_sim.nsim   = nsim;
                        cfg_sim.seqllr = dat.evid(ifilt);
                        cfg_sim.seqind = dat.sind(ifilt);
                        cfg_sim.prev   = prev_fix(itask);
                        cfg_sim.beta   = beta_fix(itask);
                        cfg_sim.brep   = 0; % force no repetition bias
                        cfg_sim.epsi   = 0; % force no lapses

                        out_sim        = sim_model_softmax(cfg_sim);

Then each subject is fit using the function "fit_model_softmax".
                        
                        % fit simulated decisions
                        for isim = 1:nsim
                            cfg_fit        = [];
                            cfg_fit.seqllr = out_sim.cfg.seqllr;
                            cfg_fit.seqind = out_sim.cfg.seqind;
                            cfg_fit.raft   = out_sim.raft(:, isim);
                            cfg_fit.rbef   = out_sim.rbef(:, isim);
                            cfg_fit.brep   = out_sim.cfg.brep; % force no repetition bias
                            cfg_fit.epsi   = out_sim.cfg.epsi; % force no lapses

                            out_recov = fit_model_softmax(cfg_fit);
                % <-- given script 
                
Estimated values of prev and beta from each simulation are stored.

                            % store model parameters
                            revindexsubj = find(seq_isubjsubsest == isubj); % isubj goes over the index limit due to the 
                                                                            % use of subsets --> reverse indexing
                            prev_recov(revindexsubj, itask, isim) = out_recov.prev;
                            beta_recov(revindexsubj, itask, isim) = out_recov.beta;
                        end    
                    end
                end

Paired-sample t-test + data saving for each given difference (ndiff).
                
                revindexsubset = find(nsubjsubset == isubjsubset);    

                [h,p,ci,stats] = ttest(prev_recov(1, 1, :), prev_recov(1, 2, :));
                prev_ttest(ittest, 1, revindexsubset) = h;
                prev_ttest(ittest, 2, revindexsubset) = p;
                prev_ttest(ittest, 3, revindexsubset) = ci(:, :, 1);
                prev_ttest(ittest, 4, revindexsubset) = ci(:, :, 2);
                prev_ttest(ittest, 5, revindexsubset) = stats.sd;

                [h,p,ci,stats] = ttest(beta_recov(1,1,:), beta_recov(1,2,:));
                beta_ttest(ittest, 1, revindexsubset) = h;
                beta_ttest(ittest, 2, revindexsubset) = p;
                beta_ttest(ittest, 3, revindexsubset) = ci(:, :, 1);
                beta_ttest(ittest, 4, revindexsubset) = ci(:, :, 2);
                beta_ttest(ittest, 5, revindexsubset) = stats.sd;

            end

            prev_diffsigni(:, :, revindexsubset) = [length(find(prev_ttest(:, 1, revindexsubset) == 1)) nttest]; %ttest counter
            beta_diffsigni(:, :, revindexsubset) = [length(find(beta_ttest(:, 1, revindexsubset) == 1)) nttest];


        end

        prev_diffsigni_acrossdiff(:, :, :, revindexdiff) = prev_diffsigni % to save data across differences tested
        beta_diffsigni_acrossdiff(:, :, :, revindexdiff) = beta_diffsigni % displayed in MATLAB LiveScript environment

        str_nsubjsubset = ''; % to nicely convert nsubjsubset into a string, e.g., "12-24-48-96"
        for isubjsubset = nsubjsubset
            str_nsubjsubset = [str_nsubjsubset, num2str(isubjsubset), '-'];
        end
        str_nsubjsubset = str_nsubjsubset(1:end-1) % to delete the last (useless) '-'

        formatSpec_file = '%s%s_(%g-%g)diff_%ssubj_%dx%dsim.mat';
    
        name_file = sprintf(formatSpec_file, dout, 'prev_ttest', idiff(1), idiff(2), str_nsubjsubset, nttest, nsim);
        save(name_file, 'dout', 'ndiff', 'idiff','nsubjsubset','nsim','nttest','prev_ttest','prev_diffsigni')

        name_file = sprintf(formatSpec_file, dout, 'beta_ttest', idiff(1), idiff(2), str_nsubjsubset, nttest, nsim);
        save(name_file, 'dout', 'ndiff', 'idiff', 'nsubjsubset', 'nsim', 'nttest', 'beta_ttest', 'beta_diffsigni')

    end



Plot of significant t-tests as a function of the tested difference for each parameter (while the value of the other is kept constant).
    
    prev_diff_tasks = prev_fit_tasks(1)-prev_fit_tasks(2) % computes the each parameter difference ...
    beta_diff_tasks = beta_fit_tasks(1)-beta_fit_tasks(2) % between the two tasks in the actual dataset 

    prevonly_diff_index = [1 2 3 4 5 6 7]    ; % to enter manually according to respective index in 'ndiff'
    finalplots_mindiff(dout, nsubjsubset, nsim, nttest, ndiff, 1, prev_diff_tasks, ...
                             prevonly_diff_index, prev_diffsigni_acrossdiff,       ...
                             'probability of reversal', 'prev', 'choice variability')

    betaonly_diff_index = [1 8 9 10 11 12 13]; % to enter manually according to respective index in 'ndiff'
    finalplots_mindiff(dout, nsubjsubset, nsim, nttest, ndiff, 2, beta_diff_tasks, ...
                              betaonly_diff_index, beta_diffsigni_acrossdiff,      ...
                              'choice variability', 'beta', 'probability of reversal') 
 
Final workscape saving (e.g., [here](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/all_diffsigni_across7diff_24subj_100x30sim.mat)).

    formatSpec_file = '%s%s_diffsigni_across%ddiff_%ssubj_%dx%dsim%s';

    name_file = sprintf(formatSpec_file, dout, 'all', length(prevonly_diff_index), str_nsubjsubset, nttest, nsim, '.mat');
    save(name_file) % to save the whole workspace

    % % To save only the variables needed for plotting:
    % name_file = sprintf(formatSpec_file, dout, 'prev', length(prevonly_diff_index), str_nsubjsubset, nttest, nsim, '.mat');
    % save(name_file, 'dout', 'nsubjsubset', 'nttest', 'ndiff', 'nsim','nttest', 'prev_diff_tasks', ...
    %                         'prev_diffsigni_acrossdiff', 'prevonly_diff_index')
    % 
    % name_file = sprintf(formatSpec_file, dout, 'beta', length(prevonly_diff_index), str_nsubjsubset, nttest, nsim, '.mat');
    % save(name_file, 'dout', 'nsubjsubset', 'nttest', 'ndiff', 'nsim','nttest', 'beta_diff_tasks', ...
    %                         'beta_diffsigni_acrossdiff', 'betaonly_diff_index')


                             
### [Final plots script (part 2)](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/finalplots_mindiff.m)

    function finalplots_mindiff (f_dout, f_nsubjsubset, f_nsim, f_nttest, f_ndiff, f_ndiffdim,      ...
                                 f_diffparam_fit, f_param_diff_index, f_param_diffsigni_acrossdiff, ... 
                                 f_strparam, f_strparam_short, f_strparam_cte)
        % Plot the number of significant ttests according to the difference tested.
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



### Conclusion (part 2)

<img src="./fig_prev_diffsigni_across7diff_24subj_100x30sim.png" width="415"> <img src="./fig_beta_diffsigni_across7diff_24subj_100x30sim.png" width="415">

Arbitrarily, a subset of 24 subjects was considered in the second part.
The plots confirm that the fitting procedure false positive rate is around 5% (when there is no difference between the simulated parameter and the fitted one) and suggest that the fitting procedure is able to acknowledge a difference of 0.05 in the probablity of reversal up to more than 90%, and a difference of 0.3 in the choice variability (close to 100%). 
Such values are the same order as the difference within parameters observed between the two tasks in the actual data (*red dotted line*, 0.2411 - 0.1763 = 0.0648 for prev and 1.3808 - 1.1071 = 0.2737 for beta).


## General conclusion
The model recovery supports the fitting procedure being reliable. 
The first part focused on the difference between a reference value for each of the two parameters of interest and their recovered values from 10 000 simulations, in four different-sized subsets (12, 24, 48 and 96 subjects). While the fitting procedure slightly overestimates both parameters, such overestimations only represent around 0.005/0.209 = 2% of the mean probability of reversal fixed value (average across the two tasks) and 0.05/1.24 = 4% of the mean choice variability fixed value (average across the two tasks).
The second part studied the sensitivity of the fitting procedure to discriminate a difference in the parameters values between the two tasks, in a subset of 24 subjects and confirmed the fitting procedure sensitivy above 90% for each parameter.

How to improve the recovery procedure:
- An easy way would be to increase the number of simulations, especially in the second part: more t-tests as well as more differences tested would sharpen the measure of sensitivity.

What this project has brought me:
- I started MATLAB language almost from scratch, I now have strong bases of it (e.g., MATLAB environment, 'MATLAB' mindset - quite different from python, the figures, etc.).

What this course has brought me:
- This course gave the tools to learn scientific python (class) and MATLAB (project), I now have the resources to build simple experiments (in python) or to conduct simple data analyses (MATLAB). Using my knowledge from one language to the other will be easier from now on.

How to improve the course:
- I would suggest to conduct the course as a tutorial, with each session being a subpart of an entire experiment. Each session would focus on different part of the script: stimuli, tasks, data collection, data analysis, and so on.

# License
All documents in this repository are distributed under the [Creative Commons Attribution-ShareAlike 4](https://creativecommons.org/licenses/by-sa/4.0/). The code is distributed under the [GPL v3 LICENSE](https://www.gnu.org/licenses/gpl-3.0.en.html).
