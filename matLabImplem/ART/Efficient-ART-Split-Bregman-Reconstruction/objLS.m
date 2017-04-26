function [cost] = objLS(x,d,A)

cost = 0.5*norm(A(x)-d)^2;