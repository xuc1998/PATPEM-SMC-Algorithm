function theta_new = randw_redraw(theta_old, LB, UB, RWM)
%RANDW_REDRAW  Box-constrained RWM proposal using reject-and-redraw.
%
% Repeatedly draws a Gaussian random-walk proposal until it falls inside
% [LB,UB]. This is intended for boundary-treatment comparison against the
% reflection-based proposal.

    d = numel(theta_old);

    if ~isfield(RWM,'Cov') || isempty(RWM.Cov)
        RWM.Cov = 1e-4 * eye(d);
    end
    if ~isfield(RWM,'Jitter') || isempty(RWM.Jitter)
        RWM.Jitter = 1e-6;
    end

    Sigma = RWM.Cov;
    if ~ismatrix(Sigma) || any(size(Sigma) ~= d) || any(~isfinite(Sigma(:)))
        Sigma = 1e-4 * eye(d);
    end

    Sigma = 0.5 * (Sigma + Sigma.');

    baseScale = max(1e-10, mean(abs(diag(Sigma))));
    jit = max(RWM.Jitter, 1e-8);
    jit = max(jit, 1e-6 * baseScale);

    [R,p] = cholcov(Sigma + jit * eye(d), 0);
    tries = 0;
    while p ~= 0 && tries < 6
        jit = jit * 10;
        [R,p] = cholcov(Sigma + jit * eye(d), 0);
        tries = tries + 1;
    end

    if p ~= 0
        diagVar = abs(diag(Sigma));
        if all(diagVar <= 0)
            diagVar = ones(d,1) * 1e-4;
        end
        [R,~] = chol(diag(diagVar) + (jit + 1e-6) * eye(d));
    end

    while true
        delta = randn(1,d) * R;
        z = theta_old + delta;
        if all(z >= LB) && all(z <= UB)
            theta_new = z;
            return;
        end
    end
end