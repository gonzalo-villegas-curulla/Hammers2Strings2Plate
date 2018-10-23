function y = phi (K, alpha, eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function phi() returns the potential from the interpenetration 'eta'
% between hammer and string. K is the stiffness constant of the collision,
% and alpha is the nonlinear exponent of the collision.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


y = K/(alpha+1)*( eta/2 + abs(eta)/2  )^(alpha+1);

end
