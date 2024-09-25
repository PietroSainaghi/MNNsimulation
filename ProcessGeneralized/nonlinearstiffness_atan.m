function redidual = nonlinearstiffness_atan(F_forresidual, Ksub_linear, Ksub_nonlinear, U)

% define function within scope

redidual = -F_forresidual + Ksub_linear * U + Ksub_nonlinear * atan(U);