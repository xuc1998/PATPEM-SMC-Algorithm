function theta_new = Generatep_redraw(particles, k, LB, UB, DEMH)
%GENERATEP_REDRAW  DE-MH proposal using reject-and-redraw.
%
% Repeatedly generates a DE-MH proposal until it falls inside [LB,UB].

    [N, d] = size(particles);

    if ~isfield(DEMH,'Gamma') || isempty(DEMH.Gamma)
        DEMH.Gamma = 2.38 / sqrt(2*d);
    end
    if ~isfield(DEMH,'NoiseSD') || isempty(DEMH.NoiseSD)
        DEMH.NoiseSD = 1e-4;
    end

    while true
        B = setdiff(1:N, k);
        R = B(randperm(N-1, 2));

        prop = particles(k,:) ...
             + DEMH.Gamma   .* (particles(R(1),:) - particles(R(2),:)) ...
             + DEMH.NoiseSD .* randn(1,d);

        if all(prop >= LB) && all(prop <= UB)
            theta_new = prop;
            return;
        end
    end
end