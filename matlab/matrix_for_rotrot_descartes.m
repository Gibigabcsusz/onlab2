function M = matrix_for_rotrot_descartes(n,d)

% the vector z goes from zero to the thickness of the plate in n steps

Mdiag  = eye(n-2)*( 2/d^2);
Mleft  = eye(n-2)*(-1/d^2);
Mright = eye(n-2)*(-1/d^2);

M = [Mleft, zeros(n-2, 2)] + ...
    [zeros(n-2, 1), Mdiag, zeros(n-2, 1)] + ...
    [zeros(n-2, 2), Mright];
