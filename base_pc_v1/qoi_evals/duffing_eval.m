function output = duffing_eval(inputs, eval_opt)
    n_runs = size(inputs,1);
    output = zeros(n_runs,1);
    t0 = 0;
    tend = 5;
    tstep = 500;
    for k = 1 : n_runs
        [T,x] = ode45(@(t,y)le_duffing(t,y, inputs(k,:)),linspace(t0,tend,tstep),[1; 0]');

        output(k) = interp1(T,x(:,1),4);
    end
end

function dxdt = le_duffing(~, x, inputs)
epsilon = inputs(1);
omega = inputs(2);
xi = inputs(3);

    dxdt = [x(2); -2*omega*xi*x(2) - omega^2*(x(1) - epsilon*x(1)^3)];
end



