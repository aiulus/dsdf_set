function [Ad, Bd] = simple_discretize(A, B, dt)
    Ad = eye(size(A)) + dt * A;
    Bd = eye(size(B)) + dt * B;
end  