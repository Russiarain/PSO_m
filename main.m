clear
rng default

% Define PSO parameters
m = 3;          % num of variables
n = 100;        % num of particles
wmax = 0.9;     % inertia weight
wmin = 0.4;     % inertia weight
c1 = 2;         % self
c2 = 2;         % social

%% PSO main
maxIteration = 1000;
maxRepeat = 10;
lb = [0 0 0];
ub = [10 10 10];
vmax = 0.1*(ub - lb);   % velocity limit
fminContainer = zeros(maxRepeat,maxIteration);
fiteContainer = zeros(maxRepeat,1);

for r = 1:maxRepeat
    % PSO Initialization
    x0 = ones(n,1)*lb + rand(n,m).*(ub-lb);
    
    x = x0;         % initial position
    v = 0.1 * x0;   % initial velocity

    % initial fitness function value
    f0 = zeros(n, 1);
    for p = 1:n
        f0(p) = object_fun(x0(p,:));
    end
    [fmin0, idx0] = min(f0);

    pbest = x0;         % Particle historical best position
    gbest = x0(idx0,:); % Global historical best position
    % PSO Initialization End

    % PSO Start
    ite = 1;
    tolerance = 1;
    while ite < maxIteration && tolerance > 1e-9
        % Update inertia weight (optional)
        % reducing inertia weight during iteration helps fine search in later iterations
        w = wmax - (wmax- wmin) * ite/maxIteration;

        % update velocity
        v = w*v + c1*rand(n,m).*(pbest-x) + c2*rand(n,m).*(gbest-x);

        % handle velocity boundary violations
        vlvio = v < -vmax;
        vuvio = v > vmax;
        for i = 1:n
            v(i,vlvio(i,:)) = -vmax(vlvio(i,:));
            v(i,vuvio(i,:)) = vmax(vuvio(i,:));
        end

        % update position
        x = x + v;

        % handle boundary violations
        lvio = x < lb;
        uvio = x > ub;
        for i = 1:n
            x(i,lvio(i,:)) = lb(lvio(i,:));
            x(i,uvio(i,:)) = ub(uvio(i,:));
        end

        % evaluate fitness
        f = zeros(n,1);
        for i = 1:n
            f(i) = object_fun(x(i,:));
        end

        % update pbest
        pbest(f<f0, :) = x(f<f0,:);
        f0(f<f0) = f(f<f0);

        % find the best particle
        [fmin, idx] = min(f0);
        % store best fitness
        fminContainer(r,ite) = fmin;
        % store iteration count
        fiteContainer(r) = ite;

        % update gbest
        if fmin < fmin0
            gbest = pbest(idx,:);
            fmin0 = fmin;
        end

        % update tolerance
%         if ite > 600
%             tolerance = abs(fminContainer(r,ite-1)-fmin0);  % question here
%         end

        ite = ite + 1;
    end

    fvalue = 10*(gbest(1)-1)^2+20*(gbest(2)-2)^2+(gbest(3)-3)^2;

    fprintf('min=%7.5f ;best=[%.5f %.5f %.5f]\n',fvalue, gbest(1), gbest(2), gbest(3))
end
