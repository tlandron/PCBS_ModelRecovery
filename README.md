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


 
