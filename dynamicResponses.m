%%% CLEAN UP AND START
clear, clc
disp('Start')

%%% INITIAL CONDITIONS AND SETUP
wn = 1;
x0 = 1;
xdot0 = 1;
t = 0:0.02:15;
zeta = 0:0.01:2;
x = zeros(length(t), length(zeta));

%%% OBTAINING RESPONSES
for k = 1:length(zeta)
    if zeta(k) < 1
        wd = wn*sqrt(1-zeta(k)^2);
        C1 = x0;
        C2 = (xdot0+wn)/wd;
        x(:,k) = C1.*exp(-zeta(k).*wn.*t).*cos(wd.*t) + C2.*exp(-zeta(k).*wn.*t).*sin(wd.*t);
    elseif zeta(k) == 1
        C1 = x0;
        C2 = x0 + wn*x0;
        x(:,k) = (C1+C2.*t).*exp(-wn.*t);
    else
        wd = wn*sqrt(zeta(k)^2-1);
        r1 = -wn*zeta(k)+wd;
        r2 = -wn*zeta(k)-wd;
        det = r2-r1;
        C1 = (r2*x0 - xdot0)/det;
        C2 = (-r1*x0 + 1)/det;
        x(:,k) = C1.*exp(r1.*t) + C2.*exp(r2.*t);
    end
end

% PLOTTING DATA
% Figure and Color Scheme
figure('Name', 'Second Order Response',...
    'Position', [50 50 1000 800])
colormap winter;

% Dynamic Response Data
mesh(zeta, t, x)
hold on

% Critically Damped Response (in red)
indexZeta1 = find(zeta == 1);
plot3(ones(length(t), 1), t, x(:, indexZeta1),...
    'Color', 'r', 'LineWidth', 4)

%%% DECORATIONS
% User Interface
grid on
rotate3d on
view(70, 30)

% Labels
title('2nd Order Linear Diff Eq Response with Varying Zeta')
xlabel('\zeta')
ylabel('time (s)')
zlabel('x(t)')

% Annotations
underdampedLabel = text(0.5,9.5,0.25,'Underdamped');
set(underdampedLabel,'BackgroundColor','w','EdgeColor','k')
critdampedLabel = text(1,9.5,0.25,'Critically Damped');
set(critdampedLabel,'BackgroundColor','w','EdgeColor','k')
overdampedLabel = text(1.75,9.5,0.25,'Overdamped');
set(overdampedLabel,'BackgroundColor','w','EdgeColor','k')


%%% DONE
disp('Done')
