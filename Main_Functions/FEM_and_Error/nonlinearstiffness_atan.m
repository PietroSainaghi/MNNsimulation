function redidual = nonlinearstiffness_atan(F_forresidual, Ksub_linear, Ksub_nonlinear, MaxForce, U)

% define function within scope

redidual = -F_forresidual + Ksub_linear * U + (MaxForce * 2 / pi) * Ksub_nonlinear * atan(U);
% redidual = -F_forresidual + Ksub_linear * U + Ksub_nonlinear * atan(U);