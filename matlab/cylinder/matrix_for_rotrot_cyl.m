function M = matrix_for_rotrot_cyl(r)

% the vector r goes from zero to the radius of the wire in n steps
% with modifed difference scheme

n = length(r);

d = r(2)-r(1);

s = linspace(d/2, r(end)-d/2, n-1); % shifted grid

M = zeros(n-2,n);

for i = 1:n-2
    M(i,i)   = r(i)/s(i); 
    M(i,i+1) = -r(i+1)*(s(i)+s(i+1))/(s(i)*s(i+1));
    M(i,i+2) = r(i+2)/s(i+1);
end

M = M/d^2;