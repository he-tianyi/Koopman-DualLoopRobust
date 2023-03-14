
function f = DCMotor_Shifted( t,x,u )

La = 0.314; Ra = 12.345; km = 0.253; J = 0.00441; B = 0.00732; taul = 1.47; ua = 60;

f = [ -(Ra/La)*x(1,:) - (km/La)*x(2,:).*u + ((km*taul)/(La*B))*u; ...
    -(B/J)*x(2,:) + (km/J)*x(1,:).*u + ((km*ua)/(J*Ra))*u];
end

% Model
%   dx1 = -(Ra/La)*x(1) - (km/La)*x(2)*u + ua/La
%   dx2   = -(B/J)*x(2) + (km/J)*x(1)*u - taul/J;
%   y = x1

% x1 - current in [-10,10]
% x2 - angular velocity in [-100,100]
% u  - control input in [-4,4]