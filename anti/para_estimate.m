function lambda = para_estimate(error, tau)

diter = 2 * tau + 1;
error = error(:);
N = length(error);
N0 = length(find(error == 0));
N1 = N - N0;

if (N1 ~= 0) & ((N1/N) > 1e-3)
     S= sum( abs(error) );
     gamma= -(N0*diter)/(2*N*diter+4*S) + ...
                    sqrt((N0*diter)^2-(2*N1*diter-4*S)*(2*N*diter+4*S))/(2*N*diter+4*S);
     lambda= -(2/diter)*log(gamma);     
end;