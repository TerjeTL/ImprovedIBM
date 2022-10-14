% T(r) = A*ln(r) + B

function [A, B] = analytical_solver(r_inner, bc_inner, inner_is_neumann, ...
    r_outer, bc_outer, outer_is_neumann)

    syms A B
    
    if (inner_is_neumann && outer_is_neumann)
        error("Cannot solve for Neumann - Neumann Boundary Conditions")
    end
    
    if (inner_is_neumann)
        eqns = [A/r_inner == bc_inner, A*log(r_outer) + B == bc_outer];
    elseif (outer_is_neumann)
        eqns = [A/r_outer == bc_outer, A*log(r_inner) + B == bc_inner];
    else
        eqns = [A*log(r_inner) + B == bc_inner, A*log(r_outer) + B == bc_outer];
    end
    
    S = vpasolve(eqns, [A B]);
    
    A = S.A;
    B = S.B;
end