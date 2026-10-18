

function main_Pendulum

    % Parameters
    g = 9.81;
    l = 1;
    m = 1;
    k = 0.5;

    % Initial conditions (theta; omega)
    theta0 = linspace(-2*pi, 2*pi, 50);
    omega0 = linspace(-7, 7, 50);   %Based on max posible omega = sqrt(4*g/l)
    
    % Simulation time
    tspan = [0 20];

    % Define stability check matrix
    stable = zeros(length(omega0), length(theta0));

    % Solve ODE
    for i = 1:length(theta0)
        for j = 1:length(omega0)

            theta0_loop = theta0(i);
            omega0_loop = omega0(j);
            x0 = [theta0_loop; omega0_loop];
            [t,x] = ode45(@rate_fun, tspan, x0);

            theta_final = x(end,1);
            omega_final = x(end,2);

            % Classify Convergence as within 0.05 (ode45 will not output exactly 0)
            stable(j,i) = abs(theta_final) < 0.05 && abs(omega_final) < 0.05;
        end
    end

    figure;
    imagesc(theta0, omega0, stable);    % Displays matrix values as a square colour map
    axis xy;
    clim([0 1]);
    xlabel('\theta_0 (rad)');
    ylabel('\omega_0 (rad/s)');
    title('Simulated Region of Attraction');
    hold on;

    x0_sample = [pi/2; 0];  % Sample trajectory (arbitrary)
    [t,x] = ode45(@rate_fun, tspan, x0_sample);
    
    [THETA, OMEGA] = meshgrid(theta0, omega0);  % Map theta/omega pairs onto grid for contour layout (specifies X/Y values)
    V_grid = 0.5*m*l^2.*OMEGA.^2 + m*g*l.*(1 - cos(THETA));

    c = 2*m*g*l;    % Plug in V (theta = pi, omega = 0) for value at upright equlibrium (perfectly upside down)
    contour(THETA, OMEGA, V_grid, [c c], 'r', 'LineWidth', 2);  % Draw contour over omega/theta grid where energy V is equal to C (upright equilibrium)
    hold off;

    % States
    theta = x(:,1);
    omega = x(:,2);

    % Plot theta
    figure;
    plot(t,theta,LineWidth=1.5);
    xlabel('Time (s)');
    ylabel('\theta (rad)');
    title('Pendulum Angle (x1)');

    % Plot omega
    figure;
    plot(t,omega,LineWidth=1.5);
    xlabel('Time (s)');
    ylabel('\omega (rad/s)');
    title('Pendulum Angular Velocity (x2)');

    % Phase portrait
    figure;
    plot(theta,omega,LineWidth=1.5);
    xlabel('\theta (rad)');
    ylabel('\omega (rad/s)');
    title('Phase Portrait');
    grid on;

    % Define Lyapunov function as sum of potential and kinetic energy
    V = (1/2) * m * l^2 .* omega.^2 + m * g * l .* (1-cos(theta));  % Use element-wise multiplication because omega and theta are in vector form
    figure;
    plot(t,V, LineWidth=1.5)
    xlabel('Time (s)')
    ylabel('Energy function V(x)')
    title('Lyapunov Function Over Time')

end

function dxdt = rate_fun(~,x)

    % Parameters
    g = 9.81;
    l = 1;
    m = 1;
    k = 0.5;

    % States
    x1 = x(1);   % theta
    x2 = x(2);   % omega

    % State equations
    dx1dt = x2;

    dx2dt = ...
        -(g/l)*sin(x1) ...
        -(k/m)*x2;

    dxdt = [dx1dt; dx2dt];

end

