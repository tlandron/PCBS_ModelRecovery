function [out] = fit_model_softmax(cfg)

% get configuration parameters
seqllr = cfg.seqllr(:); % sequence evidence
seqind = cfg.seqind(:); % sequence index within block
rbef   = cfg.rbef(:); % response preceding sequence #i
raft   = cfg.raft(:); % response following sequence #i

% set response lookup indices
nseq = numel(seqllr);
nobs = nseq;
ipost = sub2ind([nseq,2],(1:nseq)',raft);

% define model parameters
npar = 0;
pnam = {}; % parameter name
psiz = []; % parameter size
pini = []; % parameter initialization value
pmin = []; % parameter minimum value
pmax = []; % parameter maximum value
% reversal probability
npar = npar+1;
pnam{npar,1} = 'prev';
psiz(npar,1) = 1;
pini(npar,1) = 0.5;
pmin(npar,1) = 0.001;
pmax(npar,1) = 0.999;
% repetition bias
npar = npar+1;
pnam{npar,1} = 'brep';
psiz(npar,1) = 1;
pini(npar,1) = 0;
pmin(npar,1) = -10;
pmax(npar,1) = +10;
% softmax inverse temperature
npar = npar+1;
pnam{npar,1} = 'beta';
psiz(npar,1) = 1;
pini(npar,1) = 0;
pmin(npar,1) = -10;
pmax(npar,1) = +10;
% lapse rate
npar = npar+1;
pnam{npar,1} = 'epsi';
psiz(npar,1) = 1;
pini(npar,1) = 0.01;
pmin(npar,1) = 0;
pmax(npar,1) = 0.5;

% define fixed parameters
pfix = cell(npar,1);
for i = 1:npar
    if isfield(cfg,pnam{i})
        pfix{i} = reshape(cfg.(pnam{i}),[psiz(i),1]);
    end
end

% define free parameters to fit
ifit = cell(npar,1);
pfit_ini = [];
pfit_min = [];
pfit_max = [];
n = 1;
for i = 1:npar
    if isempty(pfix{i}) % free parameter
        ifit{i} = n+[1:psiz(i)]-1;
        pfit_ini = cat(1,pfit_ini,pini(i)*ones(psiz(i),1));
        pfit_min = cat(1,pfit_min,pmin(i)*ones(psiz(i),1));
        pfit_max = cat(1,pfit_max,pmax(i)*ones(psiz(i),1));
        n = n+psiz(i);
    elseif numel(pfix{i}) ~= psiz(i)
        error('wrong size for fixed parameter %s!',pnam{i});
    elseif any(isnan(pfix{i})) % partly fixed parameter
        nfit = nnz(isnan(pfix{i}));
        ifit{i} = n+[1:nfit]-1;
        pfit_ini = cat(1,pfit_ini,pini(i)*ones(nfit,1));
        pfit_min = cat(1,pfit_min,pmin(i)*ones(nfit,1));
        pfit_max = cat(1,pfit_max,pmax(i)*ones(nfit,1));
        n = n+nfit;
    end
end
nfit = length(pfit_ini);

% fit free parameters
if nfit > 0
    pval = fmincon(@fmin,pfit_ini,[],[],[],[],pfit_min,pfit_max,[], ...
        optimset('Display','notify','FunValCheck','on','Algorithm','interior-point','TolX',1e-20,'MaxFunEvals',1e6));
    [~,phat] = fmin(pval);
else
    phat = pfix;
end

% get best-fitting parameters
out = cell2struct(phat,pnam);

% save configuration structure
fnames = fieldnames(out);
fnames = cat(1,{'cfg'},fnames);
out.cfg = cfg;
out = orderfields(out,fnames);

% get quality of fit
out.nfit = nfit; % number of fitted parameters
out.nobs = nobs; % number of observations
out.llh  = get_llh(phat{:}); % model log-likelihood
out.aic  = -2*out.llh+2*nfit+2*nfit*(nfit+1)/(nobs-nfit+1); % model AIC
out.bic  = -2*out.llh+nfit*log(nobs); % model BIC

% get best-fitting posterior probabilities
[out.ppost,out.loglr] = get_ppost(phat{:});

    function [f,pfit] = fmin(p)
    % create list of model parameters
    pfit = cell(npar,1);
    for i = 1:npar
        if isempty(pfix{i}) % free parameter
            pfit{i} = p(ifit{i});
        elseif any(isnan(pfix{i})) % partly fixed parameter
            pfit{i} = pfix{i};
            pfit{i}(isnan(pfit{i})) = p(ifit{i});
        else % fully fixed parameter
            pfit{i} = pfix{i};
        end
    end
    % return negative model log-likelihood
    f = -get_llh(pfit{:});
    end

    function [llh,ppost] = get_llh(varargin)
    % get model posterior probabilities
    ppost = get_ppost(varargin{:});
    ppost = [ppost,1-ppost];
    % get model log-likelihood
    llh = sum(log(max(ppost(ipost),eps)));
    end

    function [ppost,loglr] = get_ppost(prev,brep,beta,epsi)
    loglr = nan(nseq,1); % log-belief following sequence #i
    ppost = nan(nseq,1); % posterior choice probability following sequence #i
    for iseq = 1:nseq
        loglr(iseq) = seqllr(iseq);
        if seqind(iseq) > 1
            loglr(iseq) = loglr(iseq)+update_prior(loglr(iseq-1),prev);
        end
        qval = beta*loglr(iseq);
        if seqind(iseq) > 1
            qval = qval+brep*(3-2*rbef(iseq));
        end
        ppost(iseq) = (1-epsi)./(1+exp(-qval))+0.5*epsi;
    end
    end

end

function [lo] = update_prior(li,h)
% update prior belief according to hazard rate
lo = li+log((1-h)./h+exp(-li))-log((1-h)./h+exp(+li));
end

function [p] = softmax(q,b,e)
% non-greedy softmax function
p = (1-e)./(1+exp(-b.*q))+0.5*e;
end
