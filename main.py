from Output import Output

#-------------------------- GENERAL PARAMETERS -----------------------------#
T          = 230
P          = [2529983.93, 3550022.7, 4560030.3, 5570037.9]
tol        = 1e-10
n_ite      = 3000
save       = 0
beta_v     = 0.001
case       = 0
ftol       = 1e-8
system     = 6
complete   = 0

Output.results(P=P, T=T, system=system, beta_v=beta_v, n_ite=n_ite, tol=tol, save=save,
               ftol=ftol, case=case, complete=complete)