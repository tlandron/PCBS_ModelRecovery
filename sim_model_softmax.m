function [out] = sim_model_softmax(cfg)

% get simulation parameters
nsim = cfg.nsim; % number of simulations
prev = cfg.prev; % probability of reversal = hazard rate
brep = cfg.brep; % repetition bias
beta = cfg.beta; % softmax inverse temperature
epsi = cfg.epsi; % lapse rate

% get experiment parameters
seqllr = cfg.seqllr(:); % sequence evidence
seqind = cfg.seqind(:); % sequence index within block

nseq = length(seqllr); % number of sequences
loglr = nan(nseq,nsim); % log-belief after sequence #i
resp = nan(nseq+1,nsim); % response before sequence #i
resp(1,:) = ceil(2*rand(1,nsim)); % initial guess
for iseq = 1:nseq
    loglr(iseq,:) = seqllr(iseq);
    if seqind(iseq) > 1
        loglr(iseq,:) = loglr(iseq,:)+update_prior(loglr(iseq-1,:),prev);
    end
    presp = 1./(1+exp(-beta*loglr(iseq,:)-brep*(3-2*resp(iseq,:))));
    resp(iseq+1,:) = 1+(rand(1,nsim) > presp);
end
% add lapses
if epsi > 0
    i = randsample(numel(resp),round(epsi*numel(resp)));
    resp(i) = ceil(2*rand(size(i)));
end

% create output structure
out      = [];
out.cfg  = cfg;
out.rbef = resp(1:nseq,:);
out.raft = resp(2:nseq+1,:);

end

function [lo] = update_prior(li,h)
% update prior belief according to hazard rate
lo = li+log((1-h)./h+exp(-li))-log((1-h)./h+exp(+li));
end