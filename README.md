# Programming for Cognitive and Brain Sciences: Development of model and parameter recovery procedures for the study of reversal learning
## Validating project for PCBS course (Cogmaster, 2018-2019).

The context for this recovery procedure can be in found [here](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/PCBS_ModelRecovery_Plan_TLANDRON.pdf).

## Given scripts
- for loading the data [loading script](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/finalplots.m),
    - ROLE OF THE SCRIPT, main variable out = dat;
- simulation procedure [simulation script](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/sim_model_softmax.m),
    - ROLE OF THE SCRIPT, main variable out = cfg ? ;
- fit procedure [fit script](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/fit_model_softmax.m),
    - ROLE OF THE SCRIPT, main variable out = prev, beta.

## Initial fits
A [script to run initial fits](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/initialfit.m), in order to obtain a reference value for each variable of interest (prev_fit, beta_fit) was a adapted from the [fit script](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/fit_model_softmax.m) as a function that runs for 'nrep' repetitions of the fitting procedure using the actual dataset.


    function [outf_prev_fit, outf_beta_fit] = initialfit(f_nrep, f_nsubj, f_dat, f_subjlist)
    % Run 'nrep' repetitions of the fit procedure.
    % Input:  'f_nrep'      the number of repetitions
    %         'f_nsubj'     the number of subject (integer output from loaddata.m)
    %         'f_dat'       the data (structure output from loaddata.m)
    %         'f_subjlist'  the list of subjects (array output from loaddata.m)
    % Output: 'outf_prev_fit outf_beta_fit' two 3D matrix,
    %             One matrix/parameter, structured as follow: 'outf_param_fit(isubj,itask,irep)', 
    %                 first dimension 'isubj' corresponds to each subject, 
    %                 second dimension 'itask' corresponds to each task,
    %                 third dimension 'irep' corresponds to each repetition.
    
        outf_prev_fit = zeros(f_nsubj,2,f_nrep);
        outf_beta_fit = zeros(f_nsubj,2,f_nrep);

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
                    outf_prev_fit(isubj,itask,irep) = out_fit.prev;
                    outf_beta_fit(isubj,itask,irep) = out_fit.beta;

                end
            end
        end
    end

## Part 1: Study of the mean and variance difference between the recovered values and the reference "fix" values of the parameters.

Load the data.
    
    loaddata


Inital fits ('nrep' number of fits) of the actual data to obtain an estimate of the two variables of interest, the probability of reversal (prev) and the choice probability (beta), that will serve fixed refence values to the next analyses.

    nrep = 10;

    clear prev_fit beta_fit % to prevent any dimension error when using MATLAB LiveScript environment 
    [prev_fit, beta_fit] = initialfit(nrep, nsubj, dat, subjlist);

    prev_fix = mean(prev_fit(:)) % same fixed value across the tasks
    beta_fix = mean(beta_fit(:))


In this first part, the two different task are not differentiated: the parameters are the same across the 2 tasks, resulting in (2 x nsim) simulated datasets.

    nsim              = 5000;         % number of simulations
    nsubjsubset       = [12 24 48 96] % subset(s) of participants (w/out consideration for condition, 96 = all)

    % preallocation, to save data across the subsets
    prev_recov_subjave = zeros(2*nsim,length(nsubjsubset));
    beta_recov_subjave = zeros(2*nsim,length(nsubjsubset));

    prev_diff_subjave = zeros(2*nsim,length(nsubjsubset));
    beta_diff_subjave = zeros(2*nsim,length(nsubjsubset));

    meansd_prev_diff_subjave = zeros(length(nsubjsubset),2);
    meansd_beta_diff_subjave = zeros(length(nsubjsubset),2);

    mean_corr_prev_beta = zeros(length(nsubjsubset),1);
    shared_var          = zeros(length(nsubjsubset),1);


    for isubjsubset = nsubjsubset % for each subset of subjects

        prev_recov = zeros(isubjsubset,2,nsim); % to save data across subjects
        beta_recov = zeros(isubjsubset,2,nsim);

        isubjsubset % to display the current subset running in MATLAB LiveScript environment 
        seq_isubjsubsest = randperm(nsubj,isubjsubset) % to avoid selecting deterministically the same subjects across subsets
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

            Then each datasets is fit using the (given) function "fit_model_softmax".
                % fit simulated decisions
                for isim = 1:nsim
                    cfg_fit        = [];
                    cfg_fit.seqllr = out_sim.cfg.seqllr;
                    cfg_fit.seqind = out_sim.cfg.seqind;
                    cfg_fit.raft   = out_sim.raft(:,isim);
                    cfg_fit.rbef   = out_sim.rbef(:,isim);
                    cfg_fit.brep   = out_sim.cfg.brep; % force no repetition bias
                    cfg_fit.epsi   = out_sim.cfg.epsi; % force no lapses

                    out_recov = fit_model_softmax(cfg_fit);
        % <-- given script             
            Estimated value of prev and beta from each simulation are stored and their average value is computed.
                    % store model parameters
                    revindexsubj = find(seq_isubjsubsest==isubj); % isubj goes over the index limit due to the 
                                                                  % use of subsets --> reverse indexing
                    prev_recov(revindexsubj,itask,isim) = out_recov.prev;
                    beta_recov(revindexsubj,itask,isim) = out_recov.beta;
                end    
            end
        end
        
Variables to plot histograms & histograms of the recovered parameters (if figures wanted, uncomment;  CAUTION: variables stored across subsets needed for final saving).

        prev_recov_cattask = squeeze(cat(3,prev_recov(:,1,:),prev_recov(:,2,:)))'; % concatenation of the data for the ...
        beta_recov_cattask = squeeze(cat(3,beta_recov(:,1,:),beta_recov(:,2,:)))'; % two tasks into one array (as if one task)
                                                                                   % transposed matrix for dimension reason
        revindexsubset = find(nsubjsubset==isubjsubset); 
        prev_recov_subjave(:,revindexsubset) = mean(prev_recov_cattask,2);
        beta_recov_subjave(:,revindexsubset) = mean(beta_recov_cattask,2);

    %     % figure to be played in MATLAB LiveScript environment, but not saved.
    %     figure; 
    %     hold on
    %     title(["Probability of reversal (from n = ",num2str(isubjsubset)," subjects)"]);
    %     hist(prev_recov_subjave(:,revindexsubset));
    %     hold off
    %     
    %     figure;
    %     hold on
    %     title(["Choice probability (from n = ",num2str(isubjsubset)," subjects)"]);
    %     hist(beta_recov_subjave(:,revindexsubset));
    %     hold off


Variables to plot histograms & histograms of the difference between the recovered-fixed parameters (if figures wanted, uncomment; CAUTION: variables stored across subsets needed for final saving & final plots).

        prev_diff = prev_recov - prev_fix;
        beta_diff = beta_recov - beta_fix;

        prev_diff_cattask = squeeze(cat(3, prev_diff(:, 1, :), prev_diff(:,2,:)))'; % concatenation of the data for the ...
        beta_diff_cattask = squeeze(cat(3, beta_diff(:, 1, :), beta_diff(:,2,:)))'; % two tasks into one array (as if one task)
                                                                                    % transposed matrix for dimension reason

        prev_diff_subjave(:,revindexsubset) = mean(prev_diff_cattask, 2);
        beta_diff_subjave(:,revindexsubset) = mean(beta_diff_cattask, 2);

    %     figure;
    %     hold on
    %     title(["Difference recovered-fixed probability of reversal (from n = ",num2str(isubjsubset)," subjects)"]);
    %     hist(prev_diff_subjave(:,revindexsubset));
    %     hold off
    %     
    %     figure;
    %     hold on
    %     title(["Difference recovered-fixed choice probability (from n = ",num2str(isubjsubset)," subjects)"]);
    %     hist(beta_diff_subjave(:,revindexsubset));
    %     hold off

Mean & variance of the distribution of recovered parameters across simulations for each subset of participants.
        
        meansd_prev_diff_subjave(revindexsubset,1) = mean(prev_diff_subjave(:,revindexsubset));
        meansd_prev_diff_subjave(revindexsubset,2) = std(prev_diff_subjave(:,revindexsubset));

        meansd_beta_diff_subjave(revindexsubset,1) = mean(beta_diff_subjave(:,revindexsubset));
        meansd_beta_diff_subjave(revindexsubset,2) = std(beta_diff_subjave(:,revindexsubset));

Correlation between the two parameters between subject is computed for each simulation.  
        
        prev_recov_cattask = squeeze(cat(3,prev_recov(:,1,:),prev_recov(:,2,:))); % TO CORRECT!!!
        beta_recov_cattask = squeeze(cat(3,beta_recov(:,1,:),beta_recov(:,2,:)));

        corr_prev_beta = zeros(2*nsim,1);
        for isim = 1:2*nsim
            corr_prev_beta(isim) = corr(prev_recov_cattask(:,isim), beta_recov_cattask(:,isim));
        end

        mean_corr_prev_beta(revindexsubset) = mean(corr_prev_beta)
        shared_var(revindexsubset)          = mean_corr_prev_beta(revindexsubset).^2;

    end


Saving data of interest
    
    dout = '../Out/' % folder for produced data

    str_nsubjsubset = ''; % to nicely convert nsubjsubset into a string, e.g., "12-24-48-96"
    for isubjsubset = nsubjsubset
        str_nsubjsubset = [str_nsubjsubset, num2str(isubjsubset),'-'];
    end
    str_nsubjsubset = str_nsubjsubset(1:end-1) % to delete the last useless '-'

    formatSpec_file = '%s%s_%dsim_%ssubj%s';

    name_file = sprintf(formatSpec_file, dout,'prev_recov',nsim,str_nsubjsubset,'.mat');
    save(name_file, 'nsim','prev_fix','prev_recov_subjave','prev_diff_subjave',    ...
                    'meansd_prev_diff_subjave','mean_corr_prev_beta','shared_var')

    name_file = sprintf(formatSpec_file, dout,'beta_recov',nsim,str_nsubjsubset,'.mat');
    save(name_file, 'nsim','beta_fix','beta_recov_subjave','beta_diff_subjave',    ...
                    'meansd_beta_diff_subjave','mean_corr_prev_beta','shared_var')


Histograms of the paramterer according the number of subject they were estimated from (n = 12, 24, 48 or 96).
    
    finalplots_variance(dout, nsubjsubset, nsim, str_nsubjsubset, formatSpec_file, ...
                              prev_diff_subjave, 'probability of reversal', 'prev')

    finalplots_variance(dout, nsubjsubset, nsim, str_nsubjsubset, formatSpec_file, ...
                              beta_diff_subjave, 'choice variability', 'beta')


ADD FINALS PLOTS AND CCL
 
## Part 2 : Study of the minimal difference discriminated by the fitting procedure

Load the data (given script).

    loaddata

    dout = '../Out/' % folder for produced data
    

Inital fits ('nrep' number of fits) of the actual data to obtain an estimate of the two variables of interest, the probability of reversal (prev) and the choice variability (beta), that will serve fixed refence values to the next analyses.

    nrep = 10;

    clear prev_fit beta_fit % to prevent any dimension error when using MATLAB LiveScript environment 
    [prev_fit, beta_fit] = initialfit(nrep, nsubj, dat, subjlist);

    mean_prev_fit = mean(prev_fit(:))
    mean_beta_fit = mean(beta_fit(:))


Using a (given) function "sim_model_softmax", dataset are simulated using the fixed parameters (nsim simulated datasets).

    nsim              = 30;                % number of simulations
    nttest            = 100;               % number of ttests
    ndiff             = [0      0   ; ...
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
                         0      1   ]     % couples of differences to be tested
    nsubjsubset       = [24]              % subset(s) of participants (w/out consideration for medical condition, 96 = all)

    % preallocation
    prev_diffsigni_acrossdiff = zeros(1, 2, length(nsubjsubset), 1); % significant ttest counter across the differences tested
    beta_diffsigni_acrossdiff = zeros(1, 2, length(nsubjsubset), 1);

    revindexdiff = 0
    for idiff = ndiff'
        revindexdiff = revindexdiff+1; % to facilitate counting of the difference outside the loop

        idiff % to display the difference being computed in MATLAB LiveScript environment 
        prev_fix = [mean_prev_fit-(idiff(1)/2) mean_prev_fit+(idiff(1)/2)]
        beta_fix = [mean_beta_fit-(idiff(2)/2) mean_beta_fit+(idiff(2)/2)]

        % preallocation
        prev_ttest = zeros(nttest, 5, length(nsubjsubset)); % ttest results saved
        beta_ttest = zeros(nttest, 5, length(nsubjsubset));

        prev_diffsigni = zeros(1, 2, length(nsubjsubset)); % significant ttest counter
        beta_diffsigni = zeros(1, 2, length(nsubjsubset));

        for isubjsubset = nsubjsubset
                % preallocation
                prev_recov = zeros(isubjsubset,2,nsim);
                beta_recov = zeros(isubjsubset,2,nsim);

                isubjsubset % to display the current subset running in MATLAB LiveScript environment 
                seq_isubjsubsest = randperm(nsubj,isubjsubset) % to avoid selecting deterministically the same subjects across subsets
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

Then each datasets is fit using the (given) function "fit_model_softmax".
                        
                        % fit simulated decisions
                        for isim = 1:nsim
                            cfg_fit        = [];
                            cfg_fit.seqllr = out_sim.cfg.seqllr;
                            cfg_fit.seqind = out_sim.cfg.seqind;
                            cfg_fit.raft   = out_sim.raft(:,isim);
                            cfg_fit.rbef   = out_sim.rbef(:,isim);
                            cfg_fit.brep   = out_sim.cfg.brep; % force no repetition bias
                            cfg_fit.epsi   = out_sim.cfg.epsi; % force no lapses

                            out_recov = fit_model_softmax(cfg_fit);
                % <-- given script 
                
Estimated value of prev and beta from each simulation are stored and their average value is computed.

                            % store model parameters
                            revindexsubj = find(seq_isubjsubsest==isubj); % isubj goes over the index limit due to the 
                                                                          % use of subsets --> reverse indexing
                            prev_recov(revindexsubj,itask,isim) = out_recov.prev;
                            beta_recov(revindexsubj,itask,isim) = out_recov.beta;
                        end    
                    end
                end


Paired-sample t-test + data saving for each given difference (ndiff)
                
                revindexsubset = find(nsubjsubset==isubjsubset);    

                [h,p,ci,stats] = ttest(prev_recov(1,1,:),prev_recov(1,2,:));
                prev_ttest(ittest, 1, revindexsubset) = h;
                prev_ttest(ittest, 2, revindexsubset) = p;
                prev_ttest(ittest, 3, revindexsubset) = ci(:,:,1);
                prev_ttest(ittest, 4, revindexsubset) = ci(:,:,2);
                prev_ttest(ittest, 5, revindexsubset) = stats.sd;

                [h,p,ci,stats] = ttest(beta_recov(1,1,:),beta_recov(1,2,:));
                beta_ttest(ittest, 1, revindexsubset) = h;
                beta_ttest(ittest, 2, revindexsubset) = p;
                beta_ttest(ittest, 3, revindexsubset) = ci(:,:,1);
                beta_ttest(ittest, 4, revindexsubset) = ci(:,:,2);
                beta_ttest(ittest, 5, revindexsubset) = stats.sd;

            end

            prev_diffsigni(:,:,revindexsubset) = [length(find(prev_ttest(:,1,revindexsubset) == 1)) nttest]; %ttest counter
            beta_diffsigni(:,:,revindexsubset) = [length(find(beta_ttest(:,1,revindexsubset) == 1)) nttest];


        end

        prev_diffsigni_acrossdiff(:,:,:,revindexdiff) = prev_diffsigni % to save data across differences tested
        beta_diffsigni_acrossdiff(:,:,:,revindexdiff) = beta_diffsigni % displayed in MATLAB LiveScript environment

        str_nsubjsubset = ''; % to nicely convert nsubjsubset into a string, e.g., "12-24-48-96"
        for isubjsubset = nsubjsubset
            str_nsubjsubset = [str_nsubjsubset, num2str(isubjsubset),'-'];
        end
        str_nsubjsubset = str_nsubjsubset(1:end-1) % to delete the last (useless) '-'

        formatSpec_file = '%s%s_(%g-%g)diff_%ssubj_%dx%dsim.mat';

        name_prev_ttest = sprintf(formatSpec_file, dout,'prev_ttest',idiff(1),idiff(2),str_nsubjsubset,nttest,nsim);
        save(name_prev_ttest, 'idiff','nsubjsubset','nsim','nttest','prev_ttest','prev_diffsigni')

        name_beta_ttest = sprintf(formatSpec_file, dout,'beta_ttest',idiff(1),idiff(2),str_nsubjsubset,nttest,nsim);
        save(name_beta_ttest, 'idiff','nsubjsubset','nsim','nttest','beta_ttest','beta_diffsigni')

    end



Plot of significant t-tests as a function of the difference in value for each parameter (only works while the value of the other is kept constant)
    
    prevonly_diff_index = [1 2 3 4 5 6 7]    ; % to enter manually according to the chosen 'ndiff'
    finalplots_mindiff(dout, nsubjsubset, nsim, nttest, ndiff, 1, prevonly_diff_index,     ...
                             prev_diffsigni_acrossdiff, 'probability of reversal', 'prev', ...
                             'choice variability')

    betaonly_diff_index = [1 8 9 10 11 12 13]; % to enter manually according to the chosen 'ndiff'
    finalplots_mindiff(dout, nsubjsubset, nsim, nttest, ndiff, 2, betaonly_diff_index,     ...
                             beta_diffsigni_acrossdiff, 'choice variability', 'beta',      ...
                             'probability of reversal')
                             
ADD PLOT AND CLL

## Files of interest
(different files uploaded, with different number of sim)

## Conclusion
Ccl of the model recovery

What this project has brought to me:
- starting (almost) from scrath MATLAB language
- a first hand-on on the data I'll be using for the long internship during the second semester
