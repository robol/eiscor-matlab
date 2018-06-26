function eiscor_roots(p, qz)
%EISCOR_ROOTS Compute the roots of a polynomial.
%
% X = EISCOR_ROOTS(P) computes the roots of the polynomial
%
%           P(Z) = P(1) * Z^N + ... + P(N) * Z + P(N+1)
%
%     expressed through the vectors of its coefficients in the
%     monomial basis. The algorithm used is described in [1,2].
%
% X = EISCOR_ROOTS(P, QZ) allows to select if the users wants
%     to use a QZ iteration on the companion pencil (QZ = 1),
%     or the QR iteration on the companion matrix of the monic
%     polynomial P(Z) / P(1). The latter choice if often preferable
%     (and is indeed the default), since it is faster and just
%     as stable as the QZ, at least for this very special problem.
%
% [X,Q,Z,S,T] = EISCOR_ROOTS(P, QZ) additionally returns the Schur
%     form of the companion pencil. Notice that, at the moment,
%     this is exposed only for debugging reasons, and it is
%     expected to be slow. 
%
% References
%
% [1] Aurentz, J. L., Mach, T., Vandebril, R., & Watkins, D. S. (2015). 
%     Fast and backward stable computation of roots of polynomials. 
%     SIAM Journal on Matrix Analysis and Applications, 36(3), 942-973.
% 
% [2] Aurentz, J. L., Mach, T., Robol, L., Vandebril, R., & Watkins, D. S. 
%     (2018). Fast and backward stable computation of roots of polynomials,
%     part II: backward error analysis; companion matrix and companion 
%     pencil. SIAM Journal on Matrix Analysis and Applications.

  warning('Could not file the MEX file for eiscor: trying to compile it');
  
  mex eiscor_roots.F90 ~/eiscor/lib/libeiscor.so.0.2.0
  
  fprintf([ '\n\nI tried to compile the MEX file for you. Assuming \n' ...
      'you have EISCOR installed in ~/eiscor it should have worked.\n' ...
      'If it did, you can now try to run eiscor again' ]);

end
  
