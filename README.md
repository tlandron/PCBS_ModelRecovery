# Programming for Cognitive and Brain Sciences: Development of model and parameter recovery procedures for the study of reversal learning
## Validating project for PCBS course (Cogmaster, 2018-2019).

The background for this recovery procedure can be in found [here](https://github.com/tlandron/PCBS_ModelRecovery/blob/master/PCBS_ModelRecovery_Plan_TLANDRON.pdf).

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

##Part 1: Study of the mean and variance difference between the recovered values and the reference "fix" values of the parameters.

Loading of the data
    
    loaddata

Ten initial fits using the actual data to obtain an estimate of the two variables of interest, the probability of reversal (prev) and the choice probability (beta), that will serve fixed refence values to the next analyses.

    nrep = 10;

    clear prev_fit beta_fit  % to prevent any dimension error when using the LiveScript environment in MATLAB
    [prev_fit, beta_fit] = initialfit(nrep, nsubj, dat, subjlist);

In this first part, the two different task will not be differentiated: the parameters are the same across the 2 tasks, resulting in (2 x nsim) simulated datasets.

    prev_fix = mean(prev_fit(:))
    beta_fix = mean(beta_fit(:))


Using a (given) function "sim_model_softmax", dataset are simulated using the fixed parameters).
    
    nsim              = 5000; % number of simulations
    nsubjsubset       = [12 24 48 96] % number of participants (w/out consideration for condition, 96 = all)

    % preallocation
    prev_recov_subjave = zeros(2*nsim,length(nsubjsubset));
    beta_recov_subjave = zeros(2*nsim,length(nsubjsubset));

    prev_diff_subjave = zeros(2*nsim,length(nsubjsubset));
    beta_diff_subjave = zeros(2*nsim,length(nsubjsubset));

    meansd_prev_diff_subjave = zeros(length(nsubjsubset),2);
    meansd_beta_diff_subjave = zeros(length(nsubjsubset),2);

    mean_corr_prev_beta = zeros(length(nsubjsubset),1);
    shared_var          = zeros(length(nsubjsubset),1);


    for isubjsubset = nsubjsubset
        
        % preallocation
        prev_recov = zeros(isubjsubset,2,nsim);
        beta_recov = zeros(isubjsubset,2,nsim);

        isubjsubset % to track 
        seq_isubjsubsest = randperm(nsubj,isubjsubset) % to avoid selecting deterministically the same subjects across subsets
        for isubj = seq_isubjsubsest
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

Then each dataset is fit using the (given) function "fit_model_softmax".

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

 Estimated value of prev and beta from each simulation are stored and their average value is computed.
 
                    % store model parameters
                    revindexsubj = find(seq_isubjsubsest==isubj); % isubj goes over the index limit due to subsets --> reverse indexing 
                    prev_recov(revindexsubj,itask,isim) = out_recov.prev;
                    beta_recov(revindexsubj,itask,isim) = out_recov.beta;
                end    
            end
        end
            
Histogram of the recovered parameters and their difference from the fix one.

        revindexsubset = find(nsubjsubset==isubjsubset);

        prev_recov_cattask = squeeze(cat(3,prev_recov(:,1,:),prev_recov(:,2,:)))';
        beta_recov_cattask = squeeze(cat(3,beta_recov(:,1,:),beta_recov(:,2,:)))';

        prev_recov_subjave(:,revindexsubset) = mean(prev_recov_cattask,2);
        beta_recov_subjave(:,revindexsubset) = mean(beta_recov_cattask,2);

        figure;
        hold on
        title(["Probability of reversal (from n = ",num2str(isubjsubset)," subjects)"]);
        hist(prev_recov_subjave(:,revindexsubset));
        hold off

        figure;
        hold on
        title(["Choice probability (from n = ",num2str(isubjsubset)," subjects)"]);
        hist(beta_recov_subjave(:,revindexsubset));
        hold off


        prev_diff = prev_recov - prev_fix;
        beta_diff = beta_recov - beta_fix;

        prev_diff_cattask = squeeze(cat(3, prev_diff(:, 1, :), prev_diff(:,2,:)))';
        beta_diff_cattask = squeeze(cat(3, beta_diff(:, 1, :), beta_diff(:,2,:)))';

                Parameters for each stimulation (average across subjects) is stored outside the loop.

        prev_diff_subjave(:,revindexsubset) = mean(prev_diff_cattask, 2);
        beta_diff_subjave(:,revindexsubset) = mean(beta_diff_cattask, 2);

        figure;
        hold on
        title(["Difference recovered-fixed probability of reversal (from n = ",num2str(isubjsubset)," subjects)"]);
        hist(prev_diff_subjave(:,revindexsubset));
        hold off

        figure;
        hold on
        title(["Difference recovered-fixed choice probability (from n = ",num2str(isubjsubset)," subjects)"]);
        hist(beta_diff_subjave(:,revindexsubset));
        hold off
        
           Mean & variance of the distribution of recovered parameters across simulations for each subset of participants
        meansd_prev_diff_subjave(revindexsubset,1) = mean(prev_diff_subjave(:,revindexsubset));
        meansd_prev_diff_subjave(revindexsubset,2) = std(prev_diff_subjave(:,revindexsubset));

        meansd_beta_diff_subjave(revindexsubset,1) = mean(beta_diff_subjave(:,revindexsubset));
        meansd_beta_diff_subjave(revindexsubset,2) = std(beta_diff_subjave(:,revindexsubset));
            Correlation between the two parameters between subject is computed for each simulation.  
        prev_recov_cattask = squeeze(cat(3,prev_recov(:,1,:),prev_recov(:,2,:)));
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

    str_nsubjsubset = '';
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

    Histogram of the paramterer according the number of subject they were estimated from (n = 12, 24, 48 or 96).
    finalplots_variance(dout, nsubjsubset, nsim, str_nsubjsubset, formatSpec_file, ...
               prev_diff_subjave, 'probability of reversal', 'prev')

    finalplots_variance(dout, nsubjsubset, nsim, str_nsubjsubset, formatSpec_file, ...
               beta_diff_subjave, 'choice variability', 'beta')

 
