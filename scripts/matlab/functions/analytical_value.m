% T(r) = A*ln(r) + B

function phi = analytical_value(A, B, r)
    phi = log(r).*A + B;
end