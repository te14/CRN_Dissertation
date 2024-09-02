# Examples 2.1--2.6:
print("Simplified dynamics with Hill-coefficient 2\n")

# variables
N,M,a,b,c,d,k1,k2,k3,v,u,z,fp,fpp = var('N,M,a,b,c,d,k1,k2,k3,v,u,z,fp,fpp')

# LLN
G0f = -2*k1*v*v*z*fp + (2*k2+k3)*(M-z)*fp 
G2g = -k1*v*v*z*a*fp + (k2+k3)*(M-z)*a*fp

coeff_fp(z) = (G0f+G2g).coefficient(fp)
coeff_fp_z = (coeff_fp(z) - coeff_fp(0))/z

sol_a = solve([coeff_fp_z==0], a, solution_dict=True)
a_sol =  (sol_a[0][a]).full_simplify()


mu_LLN(a) = coeff_fp(0).full_simplify()
mu_LLN_sol = mu_LLN(a_sol).full_simplify()

print("\nmu in LLN: ", mu_LLN_sol)


# CLT

L0f = 2*k1*v*v*z*fpp - 4*k1*v*u*z*fp + 2*k2*(M-z)*fpp + 1/2*k3*(M-z)*fpp
L1g = a * (-2*k1*v*v*z*(z-1)*fpp - 2*k1*v*u*z*fp + 2*k2*(M-z)*(z+1)*fpp + k3*(M-z)*(z+1)*fpp - mu_LLN(a) * z *fpp)
L2h = - k1*v*v*z * (b*fp + c*fpp + (z-1)*d*fpp) + (k2+k3)*(M-z)*(b*fp + c*fpp + z*d*fpp)

# Drift in CLT
coeff_fp(z) = (L0f+L1g+L2h).coefficient(fp)
a(z) = ((coeff_fp(z) - coeff_fp(0))/z).full_simplify()
coeff_fp_z = a(0)
sol_b = solve([coeff_fp_z==0], b, solution_dict=True)
b(a) = sol_b[0][b]
b_sol = b(a_sol)

# print("b = ", b_sol)
mu_CLT(a, b) = coeff_fp(0,0).full_simplify()
mu_CLT_sol = mu_CLT(a_sol, b_sol)
print("\nmu in CLT: ", mu_CLT(a_sol, b_sol).factor())	


# Diffusion in CLT
coeff_fpp(z) = (L0f+L1g+L2h).coefficient(fpp)
a(z) = ((coeff_fpp(z) - coeff_fpp(0))/z).full_simplify()
coeff_fpp_z = a(0)
coeff_fpp_zz = ((coeff_fpp(z) - z*a(0) - coeff_fpp(0))/z/z).full_simplify()

sol_cd = solve([coeff_fpp_z == 0, coeff_fpp_zz == 0], c, d, solution_dict=True)
c(a,b) =  (sol_cd[0][c]).full_simplify()
d(a,b) =  (sol_cd[0][d]).full_simplify()
c_sol=c(a_sol, b_sol)
d_sol=d(a_sol, b_sol)


sigma2_CLT(a, b, c, d) = coeff_fpp(0).full_simplify() 
sigma2_CLT_sol = sigma2_CLT(a_sol, b_sol, c_sol, d_sol).full_simplify().factor()
print("\n\nsigma2 in CLT: ", sigma2_CLT_sol)

# (sigma2_CLT_sol*(k1*v*v + k2 + k3)^3).full_simplify()





print("Section 3.1: Michaelis-Menten kinetics\n")

v,M,N,a,b,c,d,k1,k2,k3,p,u = var('v,M,N,a,b,c,d,k1,k2,k3,p,u')
a = -(k1*v + k2)/(k1*v+k2+k3)
d = a*a
p = k3*k1*v/(k1*v+k2+k3)

sol = solve(k1*u*(1+a) + b*(k1*v+k2+k3), b, solution_dict=True)
b = sol[0][b]
sol = solve(1/2*(k1*v-k2) + (k1*v+(M-1)*k2+M*p)*a - (k1*v+k2+k3)*c + (k1*v+M*(k2+k3))*d, c, solution_dict=True)
c = sol[0][c]

# Drift:
(M*(k2+k3)*b*(k1*v+k2+k3)^2).full_simplify()

# Diffusion:
(M*(k2*(1/2+a) + (k2+k3)*c) * 2 * (k1*v+k2+k3)^3 / k3).full_simplify()



print("Section 3.2: Dynamics with Hill-coefficient 2\n")

# variables
N,M,a1,a3,b1,b3,c1,c2,d1,d2,d3,k0,k1,k2,k3,k4,v,u,z1,z3,fp,fpp = var('N,M,a1,a3,b1,b3,c1,c2,d1,d2,d3,k0,k1,k2,k3,k4,v,u,z1,z3,fp,fpp')

# LLN
G0f = -k0*v*z1*fp + k1*(M-z1-z3)*fp - k2*v*(M-z1-z3)*fp + k3*z3*fp
G1g = -k0*v*z1*a1*fp + k1*(M-z1-z3)*a1*fp + k2*v*(M-z1-z3)*a3*fp - (k3+k4)*z3*a3*fp 

coeff_fp(z1,z3) = (G0f+G1g).coefficient(fp)
a(z1) = ((coeff_fp(z1,0) - coeff_fp(0,0))/z1).full_simplify()
coeff_fp_z1 = a(0)
a(z3) = ((coeff_fp(0,z3) - coeff_fp(0,0))/z3).full_simplify()
coeff_fp_z3 = a(0)

sol_a = solve([coeff_fp_z1==0, coeff_fp_z3==0], a1, a3, solution_dict=True)
a1_sol =  (sol_a[0][a1]).full_simplify()
a3_sol =  (sol_a[0][a3]).full_simplify()
# print("a1 = ", a1_sol)
# print("a3 = ", a3_sol)
mu_LLN(a1, a3) = coeff_fp(0,0).full_simplify()
mu_LLN_sol = mu_LLN(a1_sol, a3_sol)
print("\nmu in LLN: ", mu_LLN(a1_sol, a3_sol).factor())	
mu_LLN_fromk(k0,k1,k2,k3,k4) = mu_LLN(a1_sol, a3_sol)
print("...as N->infty: ", limit(mu_LLN_fromk(k0,N*k1,N*k2,k3,k4), N=oo))

# CLT

L0f = -k0*u*z1*fp + 1/2*k0*v*z1*fpp + 1/2*k1*(M-z1-z3)*fpp - k2*(M-z1-z3)*u*fp + 1/2 * k2*v*(M-z1-z3)*fpp + 1/2*k3*z3*fpp
L1g = -k0*z1*v*(a1*(z1-1) + a3*z3)*fpp - k0*z1*u*a1*fp + k1*(M-z1-z3)*(a1*(z1+1)+a3*z3)*fpp - k2*(M-z1-z3)*v*(a1*z1+a3*(z3+1))*fpp + k2*(M-z1-z3)*u*a3*fp + k3*z3*(a1*z1+a3*(z3-1))*fpp - mu_LLN(a1,a3)*(a1*z1+a3*z3)*fpp
L2h = -k0*z1*v*(b1*fp + (c1+(z1-1)*d1+(z1+z3-1)*d3)*fpp) + k1*(M-z1-z3)*(b1*fp + (c1+z1*d1+(z1+z3)*d3)*fpp) + k2*(M-z1-z3)*v*(b3*fp + (c2 + z3*d2 + (z1+z3)*d3)*fpp) - (k3+k4)*z3*(b3*fp + (c2+(z3-1)*d2+(z1+z3-1)*d3)*fpp)

# Drift in CLT
coeff_fp(z1,z3) = (L0f+L1g+L2h).coefficient(fp)
a(z1) = ((coeff_fp(z1,0) - coeff_fp(0,0))/z1).full_simplify()
coeff_fp_z1 = a(0)
a(z3) = ((coeff_fp(0,z3) - coeff_fp(0,0))/z3).full_simplify()
coeff_fp_z3 = a(0)
sol_b = solve([coeff_fp_z1==0, coeff_fp_z3==0], b1, b3, solution_dict=True)
b1(a1,a3) =  (sol_b[0][b1]).full_simplify()
b3(a1,a3) =  (sol_b[0][b3]).full_simplify()
b1_sol = b1(a1_sol, a3_sol)
b3_sol = b3(a1_sol, a3_sol)
# print("b1 = ", b1_sol)
# print("b3 = ", b3_sol, "\n")
mu_CLT(a1, a3, b1, b3) = coeff_fp(0,0).full_simplify()
mu_CLT_sol = mu_LLN(a1_sol, a3_sol, b1_sol, b3_sol)
print("\nmu in CLT: ", mu_CLT(a1_sol, a3_sol, b1_sol, b3_sol).factor())	
mu_CLT_fromk(k0,k1,k2,k3,k4) = mu_CLT(a1_sol, a3_sol, b1_sol, b3_sol).full_simplify()
print("...as N->infty: ", limit(mu_CLT_fromk(k0,N*k1,N*k2,k3,k4), N=oo).factor())


# Diffusion in CLT
coeff_fpp(z1,z3) = (L0f+L1g+L2h).coefficient(fpp)
a(z1) = ((coeff_fpp(z1,0) - coeff_fpp(0,0))/z1).full_simplify()
coeff_fpp_z1 = a(0)
coeff_fpp_z1z1 = ((coeff_fpp(z1,0) - z1*a(0) - coeff_fpp(0,0))/z1/z1).full_simplify()
a(z3) = ((coeff_fpp(0,z3) - coeff_fpp(0,0))/z3).full_simplify()
coeff_fpp_z3 = a(0)
coeff_fpp_z3z3 = ((coeff_fpp(0,z3) - z3*a(0) - coeff_fpp(0,0))/z3/z3).full_simplify()
coeff_fpp_z1z3 = ((coeff_fpp(z1,z3) - z1*coeff_fpp_z1 - z3*coeff_fpp_z3 - z1^2*coeff_fpp_z1z1 - z3^2*coeff_fpp_z3z3 - coeff_fpp(0,0))/z1/z3).full_simplify()

sol_c = solve([coeff_fpp_z1 == 0, coeff_fpp_z3 == 0, coeff_fpp_z1z1 == 0, coeff_fpp_z3z3 == 0, coeff_fpp_z1z3 == 0], c1, c2, d1, d2, d3, solution_dict=True)
c1(a1,a3,b1,b3) =  (sol_c[0][c1]).full_simplify()
c2(a1,a3,b1,b3) =  (sol_c[0][c2]).full_simplify()
d1(a1,a3,b1,b3) =  (sol_c[0][d1]).full_simplify()
d2(a1,a3,b1,b3) =  (sol_c[0][d2]).full_simplify()
d3(a1,a3,b1,b3) =  (sol_c[0][d3]).full_simplify()
c1_sol=c1(a1_sol, a3_sol, b1_sol, b3_sol)
c2_sol=c2(a1_sol, a3_sol, b1_sol, b3_sol)
d1_sol=d1(a1_sol, a3_sol, b1_sol, b3_sol)
d2_sol=d2(a1_sol, a3_sol, b1_sol, b3_sol)
d3_sol=d3(a1_sol, a3_sol, b1_sol, b3_sol)

sigma2_CLT(a1, a3, b1, b3, c1, c2, d1, d2, d3) = coeff_fpp(0,0).full_simplify() 
sigma2_CLT_sol = sigma2_CLT(a1_sol, a3_sol, b1_sol, b3_sol, c1_sol, c2_sol, d1_sol, d2_sol, d3_sol).full_simplify().factor()
print("\n\nsigma2 in CLT: ", sigma2_CLT_sol)
sigma2_CLT_fromk(k0,k1,k2,k3,k4) = sigma2_CLT_sol
print("...as N->infty: ", limit(sigma2_CLT_fromk(k0,N*k1,N*k2,k3,k4), N=oo).factor())




print("3.3 Extension of Section 4.1 to two fast species\n")
a1,a2,b1,b2,d1,d2,d3,k0,k1,k2,k3,k4,v,u = var('a1,a2,b1,b2,d1,d2,d3,k0,k1,k2,k3,k4,v,u')

sol_a = solve([(k1+k2)*a1 - k2*a2 == k1, k3*a1 - (k3+k4)*a2 == 0], a1, a2, solution_dict=True)
a1 = sol_a[0][a1]
a2 = sol_a[0][a2]

sol_b = solve([1/2*k1 - v*k0*a1*a1 - (k1+k2)*b1 + k2*b2 + v*k0*d1 + v*k0*d3 == 0, -v*k0*a1*a2 + k3*b1 - (k3+k4)*b2 + v*k0*d3 == 0, k1*a1 - (k1 + k2)*d1 - k1*d3 == 0, (k3+k4)*d2 + k4*d3 == 0, k1*a2 + k3*d1 + k2*d2 - (k1+k4)*d3 == 0], b1, b2, d1, d2, d3, solution_dict=True)
b1 = sol_b[0][b1]

print("b1 = ")
b1
print("\n\n LLN:\n\n")
(k0*v*(a1-1)).full_simplify()
print("\n\nCLT:\n\n")
print("Drift:\n")
(k0*u*(a1-1)).full_simplify()
print("\nDiffusion:\n")
(1/2 * k0*v*(1-2*a1+2*b1)).full_simplify()





print("Section 3.4: Main example from Cappelletti and Wiuf\n")

# variables
v1,v3,u1,u3,z1,z2,N,M,a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,c5,c6,d1,d2,d3,d4,d5,d6,d7,d8,d9,k1,k2,k3,k4,k5,k6,k7,k8,D1f,D3f,D11f,D13f,D33f = var('v1,v3,u1,u3,z1,z2,N,M,a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,c5,c6,d1,d2,d3,d4,d5,d6,d7,d8,d9,k1,k2,k3,k4,k5,k6,k7,k8,D1f,D3f,D11f,D13f,D33f')

# LLN

G0f = - v1*z1*(k1*D1f + k2*D1f) + k3*z2*D1f + k6*(M-z1-z2)*D3f + k8*v3*z1*(D1f-D3f)
g(z1,z2) = (a1*z1+a2*z2)*D1f + (a3*z1+a4*z2)*D3f
G1g = k1*v1*z1*(g(z1-1,z2+1)-g(z1,z2)) + k2*v1*z1*(g(z1-1,z2)-g(z1,z2)) + k3*z2*(g(z1+1,z2-1)-g(z1,z2)) + k4*z2*(g(z1+1,z2-1)-g(z1,z2)) + k5*z2*(g(z1,z2-1)-g(z1,z2)) + k6*(M-z1-z2)*(g(z1+1,z2)-g(z1,z2)) + k7*(M-z1-z2)*(g(z1,z2+1)-g(z1,z2))

coeff_LLN_D1f(z1,z2) = (G0f+G1g).coefficient(D1f)
a(z1) = ((coeff_LLN_D1f(z1,0) - coeff_LLN_D1f(0,0))/z1).full_simplify()
coeff_LLN_D1f_z1 = a(0)
a(z2) = ((coeff_LLN_D1f(0,z2) - coeff_LLN_D1f(0,0))/z2).full_simplify()
coeff_LLN_D1f_z2 = a(0)
coeff_LLN_D3f(z1,z2) = (G0f+G1g).coefficient(D3f)
a(z1) = ((coeff_LLN_D3f(z1,0) - coeff_LLN_D3f(0,0))/z1).full_simplify()
coeff_LLN_D3f_z1 = a(0)
a(z2) = ((coeff_LLN_D3f(0,z2) - coeff_LLN_D3f(0,0))/z2).full_simplify()
coeff_LLN_D3f_z2 = a(0)

sol_ab = solve([coeff_LLN_D1f_z1==0, coeff_LLN_D1f_z2==0, coeff_LLN_D3f_z1==0, coeff_LLN_D3f_z2==0], a1, a2, a3, a4, solution_dict=True)
a1_sol =  (sol_ab[0][a1]).full_simplify()
a2_sol =  (sol_ab[0][a2]).full_simplify()
a3_sol =  (sol_ab[0][a3]).full_simplify()
a4_sol =  (sol_ab[0][a4]).full_simplify()

mu_LLN_1(a1, a2, a3, a4) = coeff_LLN_D1f(0,0).full_simplify()
mu_LLN_1_sol = mu_LLN_1(a1_sol, a2_sol, a3_sol, a4_sol)
mu_LLN_3(a1, a2, a3, a4) = coeff_LLN_D3f(0,0).full_simplify()
mu_LLN_3_sol = mu_LLN_3(a1_sol, a2_sol, a3_sol, a4_sol)
mu_LLN_1_sol_fromk(k1,k2,k3,k4,k5,k6,k7,k8) = mu_LLN_1_sol
mu_LLN_3_sol_fromk(k1,k2,k3,k4,k5,k6,k7,k8) = mu_LLN_3_sol

print("LLN:")
print("coeff of D1f as N->infty: ", limit(mu_LLN_1_sol_fromk(k1,k2,k3,N^2*k4,N^2*k5,N*k6,N*k7,k8), N=oo).full_simplify())
print("coeff of D3f as N->infty: ", limit(mu_LLN_3_sol_fromk(k1,k2,k3,N^2*k4,N^2*k5,N*k6,N*k7,k8), N=oo).full_simplify())

lambda1,lambda2 = var('lambda1, lambda2')
lambda1 = (k1*k4*(k6+k7) + k2*k4*k7)/(k4*(k6+k7)+k5*k6)
lambda2 = (k2*k6*(k4+k5) + k1*k5*k6)/(k4*(k6+k7)+k5*k6)




# CLT
L0f = -k1*u1*z1*D1f + 1/2*k1*v1*z1*D11f - k2*u1*z1*D1f + 1/2*k2*v1*z1*D11f + 1/2*k3*z2*D11f + 1/2*k6*(M-z1-z2)*D33f + k8*u3*z1*(D1f - D3f) + 1/2*k8*v3*z1*(D11f - 2*D13f + D33f)
g(z1,z2) = (a1*z1+a2*z2)*D1f + (a3*z1+a4*z2)*D3f
D1g(z1,z2) = (a1*z1+a2*z2)*D11f + (a3*z1+a4*z2)*D13f
D3g(z1,z2) = (a1*z1+a2*z2)*D13f + (a3*z1+a4*z2)*D33f
L1g = -k1*v1*z1*D1g(z1-1,z2+1) + k1*u1*z1*(g(z1-1,z2+1)-g(z1,z2)) - k2*v1*z1*D1g(z1-1,z2) + k2*u1*z1*(g(z1-1,z2)-g(z1,z2)) + k3*z2*D1g(z1+1,z2-1) + k6*(M-z1-z2)*D3g(z1+1,z2) + k8*v3*z1*(D1g(z1,z2) - D3g(z1,z2)) - M*((-(lambda1+lambda2)*v1+k8*v3)*D1g(z1,z2) + (lambda2*v1-k8*v3)*D3g(z1,z2))

h(z1,z2)  = (b1*z1+b2*z2)*D1f + (b3*z1+b4*z2)*D3f + (c1*z1+c2*z2+d1*z1*(z1-1)/2 + d2*z2*(z2-1)/2 + d3*(z1+z2)*(z1+z2-1)/2) * D11f + (c3*z1+c4*z2+d4*z1*(z1-1)/2 + d5*z2*(z2-1)/2 + d6*(z1+z2)*(z1+z2-1)/2) * D33f + (c5*z1+c6*z2+d7*z1*(z1-1)/2 + d8*z2*(z2-1)/2 + d9*(z1+z2)*(z1+z2-1)/2) * D13f

L2h = k1*v1*z1*(h(z1-1,z2+1)-h(z1,z2)) + k2*v1*z1*(h(z1-1,z2)-h(z1,z2)) + k3*z2*(h(z1+1,z2-1)-h(z1,z2)) + k4*z2*(h(z1+1,z2-1)-h(z1,z2)) + k5*z2*(h(z1,z2-1)-h(z1,z2)) + k6*(M-z1-z2)*(h(z1+1,z2) - h(z1,z2)) + k7*(M-z1-z2)*(h(z1,z2+1)-h(z1,z2))



# Drift in CLT
coeff_CLT_D1f(z1,z2) = (L0f+L1g+L2h).coefficient(D1f)
coeff_CLT_D3f(z1,z2) = (L0f+L1g+L2h).coefficient(D3f)
a(z1) = ((coeff_CLT_D1f(z1,0) - coeff_CLT_D1f(0,0))/z1).full_simplify()
coeff_CLT_D1f_z1 = a(0)
a(z2) = ((coeff_CLT_D1f(0,z2) - coeff_CLT_D1f(0,0))/z2).full_simplify()
coeff_CLT_D1f_z2 = a(0)
a(z1) = ((coeff_CLT_D3f(z1,0) - coeff_CLT_D3f(0,0))/z1).full_simplify()
coeff_CLT_D3f_z1 = a(0)
a(z2) = ((coeff_CLT_D3f(0,z2) - coeff_CLT_D3f(0,0))/z2).full_simplify()
coeff_CLT_D3f_z2 = a(0)

sol_b = solve([coeff_CLT_D1f_z1==0, coeff_CLT_D1f_z2==0, coeff_CLT_D3f_z1==0, coeff_CLT_D3f_z2==0], b1, b2, b3, b4, solution_dict=True)
b1(a1,a2,a3,a4) =  (sol_b[0][b1]).full_simplify()
b2(a1,a2,a3,a4) =  (sol_b[0][b2]).full_simplify()
b3(a1,a2,a3,a4) =  (sol_b[0][b3]).full_simplify()
b4(a1,a2,a3,a4) =  (sol_b[0][b4]).full_simplify()
b1_sol = b1(a1_sol, a2_sol, a3_sol, a4_sol).full_simplify()
b2_sol = b2(a1_sol, a2_sol, a3_sol, a4_sol).full_simplify()
b3_sol = b3(a1_sol, a2_sol, a3_sol, a4_sol).full_simplify()
b4_sol = b4(a1_sol, a2_sol, a3_sol, a4_sol).full_simplify()

mu_CLT_D1f(a1, a2, a3, a4, b1, b2, b3, b4) = coeff_CLT_D1f(0,0).full_simplify()
mu_CLT_D1f_sol = mu_CLT_D1f(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol)
mu_CLT_D3f(a1, a2, a3, a4, b1, b2, b3, b4) = coeff_CLT_D3f(0,0).full_simplify()
mu_CLT_D3f_sol = mu_CLT_D3f(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol)
mu_CLT_D1f_sol_fromk(k1,k2,k3,k4,k5,k6,k7,k8) = mu_CLT_D1f_sol
mu_CLT_D3f_sol_fromk(k1,k2,k3,k4,k5,k6,k7,k8) = mu_CLT_D3f_sol

print("CLT:")
print("coeff of D1f as N->infty: ", limit(mu_CLT_D1f_sol_fromk(k1,k2,k3,N^2*k4,N^2*k5,N*k6,N*k7,k8), N=oo))
print("coeff of D3f as N->infty: ", limit(mu_CLT_D3f_sol_fromk(k1,k2,k3,N^2*k4,N^2*k5,N*k6,N*k7,k8), N=oo))


# Diffusiont in CLT
coeff_CLT_D11f(z1,z2) = (L0f+L1g+L2h).coefficient(D11f).full_simplify()
coeff_CLT_D33f(z1,z2) = (L0f+L1g+L2h).coefficient(D33f).full_simplify()
coeff_CLT_D13f(z1,z2) = (L0f+L1g+L2h).coefficient(D13f).full_simplify()

a(z1) = ((coeff_CLT_D11f(z1,0) - coeff_CLT_D11f(0,0))/z1).full_simplify()
coeff_CLT_D11f_z1 = a(0).full_simplify()
coeff_CLT_D11f_z1z1 = ((coeff_CLT_D11f(z1,0) - z1*a(0) - coeff_CLT_D11f(0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D11f(0,z2) - coeff_CLT_D11f(0,0))/z2).full_simplify()
coeff_CLT_D11f_z2 = a(0).full_simplify()
coeff_CLT_D11f_z2z2 = ((coeff_CLT_D11f(0,z2) - z2*a(0) - coeff_CLT_D11f(0,0))/z2/z2).full_simplify()
coeff_CLT_D11f_z1z2 = ((coeff_CLT_D11f(z1,z2) - z1*coeff_CLT_D11f_z1 - z2*coeff_CLT_D11f_z2 - z1^2*coeff_CLT_D11f_z1z1 - z2^2*coeff_CLT_D11f_z2z2 - coeff_CLT_D11f(0,0))/z1/z2).full_simplify()

a(z1) = ((coeff_CLT_D33f(z1,0) - coeff_CLT_D33f(0,0))/z1).full_simplify()
coeff_CLT_D33f_z1 = a(0).full_simplify()
coeff_CLT_D33f_z1z1 = ((coeff_CLT_D33f(z1,0) - z1*a(0) - coeff_CLT_D33f(0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D33f(0,z2) - coeff_CLT_D33f(0,0))/z2).full_simplify()
coeff_CLT_D33f_z2 = a(0).full_simplify()
coeff_CLT_D33f_z2z2 = ((coeff_CLT_D33f(0,z2) - z2*a(0) - coeff_CLT_D33f(0,0))/z2/z2).full_simplify()
coeff_CLT_D33f_z1z2 = ((coeff_CLT_D33f(z1,z2) - z1*coeff_CLT_D33f_z1 - z2*coeff_CLT_D33f_z2 - z1^2*coeff_CLT_D33f_z1z1 - z2^2*coeff_CLT_D33f_z2z2 - coeff_CLT_D33f(0,0))/z1/z2).full_simplify()

a(z1) = ((coeff_CLT_D13f(z1,0) - coeff_CLT_D13f(0,0))/z1).full_simplify()
coeff_CLT_D13f_z1 = a(0).full_simplify()
coeff_CLT_D13f_z1z1 = ((coeff_CLT_D13f(z1,0) - z1*a(0) - coeff_CLT_D13f(0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D13f(0,z2) - coeff_CLT_D13f(0,0))/z2).full_simplify()
coeff_CLT_D13f_z2 = a(0).full_simplify()
coeff_CLT_D13f_z2z2 = ((coeff_CLT_D13f(0,z2) - z2*a(0) - coeff_CLT_D13f(0,0))/z2/z2).full_simplify()
coeff_CLT_D13f_z1z2 = ((coeff_CLT_D13f(z1,z2) - z1*coeff_CLT_D13f_z1 - z2*coeff_CLT_D13f_z2 - z1^2*coeff_CLT_D13f_z1z1 - z2^2*coeff_CLT_D13f_z2z2 - coeff_CLT_D13f(0,0))/z1/z2).full_simplify()

sol_cA = solve([coeff_CLT_D11f_z1==0, coeff_CLT_D11f_z2==0, coeff_CLT_D11f_z1z1==0, coeff_CLT_D11f_z2z2==0, coeff_CLT_D11f_z1z2==0], c1, c2, d1, d2, d3, solution_dict=True)
c1(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cA[0][c1]).full_simplify()
c2(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cA[0][c2]).full_simplify()
d1(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cA[0][d1]).full_simplify()
d2(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cA[0][d2]).full_simplify()
d3(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cA[0][d3]).full_simplify()
c1_sol = c1(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
c2_sol = c2(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
d1_sol = d1(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
d2_sol = d2(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
d3_sol = d3(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()

sol_cB = solve([coeff_CLT_D33f_z1==0, coeff_CLT_D33f_z2==0, coeff_CLT_D33f_z1z1==0, coeff_CLT_D33f_z2z2==0, coeff_CLT_D33f_z1z2==0], c3, c4, d4, d5, d6, solution_dict=True)
c3(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cB[0][c3]).full_simplify()
c4(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cB[0][c4]).full_simplify()
d4(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cB[0][d4]).full_simplify()
d5(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cB[0][d5]).full_simplify()
d6(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cB[0][d6]).full_simplify()
c3_sol = c3(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
c4_sol = c4(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
d4_sol = d4(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
d5_sol = d5(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
d6_sol = d6(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()


sol_cC = solve([coeff_CLT_D13f_z1==0, coeff_CLT_D13f_z2==0, coeff_CLT_D13f_z1z1==0, coeff_CLT_D13f_z2z2==0, coeff_CLT_D13f_z1z2==0], c5, c6, d7, d8, d9, solution_dict=True)
c5(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cC[0][c5]).full_simplify()
c6(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cC[0][c6]).full_simplify()
d7(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cC[0][d7]).full_simplify()
d8(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cC[0][d8]).full_simplify()
d9(a1,a2,a3,a4,b1,b2,b3,b4) =  (sol_cC[0][d9]).full_simplify()
c5_sol = c5(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
c6_sol = c6(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
d7_sol = d7(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
d8_sol = d8(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()
d9_sol = d9(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol).full_simplify()

sigma2_CLT_D11f(a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,d1,d2,d3,c3,c4,d4,d5,d6,c5,c6,d7,d8,d9) = coeff_CLT_D11f(0,0).full_simplify()
sigma2_CLT_D11f_sol = sigma2_CLT_D11f(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol, c1_sol, c2_sol, d1_sol, d2_sol, d3_sol, c3_sol, c4_sol, d4_sol, d5_sol, d6_sol, c5_sol, c6_sol, d7_sol, d8_sol, d9_sol).full_simplify()
sigma2_CLT_D33f(a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,d1,d2,d3,c3,c4,d4,d5,d6,c5,c6,d7,d8,d9) = coeff_CLT_D33f(0,0).full_simplify()
sigma2_CLT_D33f_sol = sigma2_CLT_D33f(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol, c1_sol, c2_sol, d1_sol, d2_sol, d3_sol, c3_sol, c4_sol, d4_sol, d5_sol, d6_sol, c5_sol, c6_sol, d7_sol, d8_sol, d9_sol).full_simplify()
sigma2_CLT_D13f(a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,d1,d2,d3,c3,c4,d4,d5,d6,c5,c6,d7,d8,d9) = coeff_CLT_D13f(0,0).full_simplify()
sigma2_CLT_D13f_sol = sigma2_CLT_D13f(a1_sol, a2_sol, a3_sol, a4_sol, b1_sol, b2_sol, b3_sol, b4_sol, c1_sol, c2_sol, d1_sol, d2_sol, d3_sol, c3_sol, c4_sol, d4_sol, d5_sol, d6_sol, c5_sol, c6_sol, d7_sol, d8_sol, d9_sol).full_simplify()

sigma2_CLT_D11f_sol_fromk(k1,k2,k3,k4,k5,k6,k7,k8) = sigma2_CLT_D11f_sol
sigma2_CLT_D33f_sol_fromk(k1,k2,k3,k4,k5,k6,k7,k8) = sigma2_CLT_D33f_sol
sigma2_CLT_D13f_sol_fromk(k1,k2,k3,k4,k5,k6,k7,k8) = sigma2_CLT_D13f_sol

print("coeff of D11f as N->infty: ", limit(sigma2_CLT_D11f_sol_fromk(k1,k2,k3,N^2*k4,N^2*k5,N*k6,N*k7,k8), N=oo))
print("coeff of D33f as N->infty: ", limit(sigma2_CLT_D33f_sol_fromk(k1,k2,k3,N^2*k4,N^2*k5,N*k6,N*k7,k8), N=oo))
print("coeff of D13f as N->infty: ", limit(sigma2_CLT_D13f_sol_fromk(k1,k2,k3,N^2*k4,N^2*k5,N*k6,N*k7,k8), N=oo))



print("Section 3.7: Michaelis-Menten with 2 compartments, gamma>beta")

# variables
M,N,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,d5,d6,k11,k12,k13,k21,k22,k23,n11,n12,n21,n22,n31,n32,v1,v2,u2,z1,z2,z3,fp,fp2,fpp = var('M,N,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,d5,d6,k11,k12,k13,k21,k22,k23,n11,n12,n21,n22,n31,n32,v1,v2,u2,z1,z2,z3,fp,fp2,fpp')

k11=k13=k21=k22=k23=n12=n22=n32=n11=n21=n31=1
k12=2
# LLN(beta<gamma)
v1=v2*n32/(n31+n32)
G0f = (-v1*z1*k11+z2*k12-(v2-v1)*(z3-z1)*k21+(M-z2-z3)*k22)*fp
G2g = (-v1*z1*k11+z2*k12+z2*k13-(v2-v1)*(z3-z1)*k21+(M-z2-z3)*k22+(M-z2-z3)*k23)*(b1*fp)+(v1*z1*k11-z2*k12-z2*k13)*c1*fp

G2h1 = (v1*z1*k11*(-b1+c1)+z2*k12*(b1-c1)+z2*k13*(b1-c1))*fp*0
G2h2 = (z1*n11*(-a1)+(z3-z1)*n12*a1+z2*n21*(-a1)+(M-z2-z3)*n22*a1)*fp

coeff_fp(z1,z2,z3) = (G0f+G2g+G2h1+G2h2).coefficient(fp)
a(z1) = ((coeff_fp(z1,0,0) - coeff_fp(0,0,0))/z1).full_simplify()
coeff_fp_z1 = a(0)
a(z2) = ((coeff_fp(0,z2,0) - coeff_fp(0,0,0))/z2).full_simplify()
coeff_fp_z2 = a(0)
a(z3) = ((coeff_fp(0,0,z3) - coeff_fp(0,0,0))/z3).full_simplify()
coeff_fp_z3 = a(0)


sol_a = solve([coeff_fp_z1==0, coeff_fp_z2==0, coeff_fp_z3==0], a1, b1,
c1, solution_dict=True)
a1_sol =  (sol_a[0][a1]).full_simplify()
b1_sol =  (sol_a[0][b1]).full_simplify()
c1_sol =  (sol_a[0][c1]).full_simplify()
mu_LLN(a1, b1, c1) = coeff_fp(0,0,0).full_simplify()
mu_LLN_sol = mu_LLN(a1_sol, b1_sol, c1_sol)
print("\nmu in LLN: ", mu_LLN(a1_sol, b1_sol, c1_sol).factor())	
mu_LLN_fromk(k11,k12,k13,k21,k22,k23,n11,n12,n21,n22,n31,n32) = mu_LLN(a1_sol, b1_sol, c1_sol)
print("...as N->infty: ",
limit(mu_LLN_fromk(k11,k12,k13,k21,k22,k23,N^-0.5*n11,N^-0.5*n12,N^-0.5*n21,N^-0.5*n22,n31,n32),
N=oo).full_simplify())

k11=k13=k21=k22=k23=n12=n22=n32=n11=n21=n31=1
k12=2
# CLT(beta<gamma)
v1=n32/(n31+n32)*v2
L0f = 1/2*(v1*z1*k11+z2*k12+(v2-v1)*(z3-z1)*k21+(M-z2-z3)*k22)*fpp-u2*n31/(n31+n32)*(z3-z1)*k21*fp-u2*z1*k11*n32/(n31+n32)*fp
L1g = (-v1*z1*k11*(b1*(z3-1)+c1*(z2+1))+z2*k12*(b1*(z3+1)+c1*(z2-1))-(v2-v1)*(z3-z1)*k21*(b1*(z3-1)+c1*z2)+(M-z2-z3)*k22*(b1*(z3+1)+c1*z2))*fpp+u2*n31/(n31+n32)*(z3-z1)*k21*(-b1)*fp+u2*z1*k11*n32/(n31+n32)*(-b1+c1)*fp-mu_LLN(a1,b1,c1)*(b1*z3+c1*z2)*fpp
L2h1 = v1*z1*k11*((-b2+c2)*fp+(-b3+c3+d1*(-z1+1)+d2*(-z3+1)+d3*z2+d5*(-2*z1-2*z3+3))*fpp)+z2*k12*((b2-c2)*fp+(b3-c3+d1*z1+d2*z3+d3*(-z2+1)+d5*(2*z3+2*z1+1))*fpp)+z2*k13*((b2-c2)*fp+(b3-c3+d1*z1+d2*z3+d3*(-z2+1)+d5*(2*z3+2*z1+1))*fpp)+(v2-v1)*(z3-z1)*k21*(-b2*fp+(-b3+d2*(-z3+1)+d4*(-z2-z3+1)+d5*(-z3-z1+1))*fpp)+(M-z2-z3)*k22*(b2*fp+(b3+d2*z3+d4*(z2+z3)+d5*(z1+z3))*fpp)+(M-z2-z3)*k23*(b2*fp+(b3+d2*z3+d4*(z2+z3)+d5*(z1+z3))*fpp)
L2h2 = (z1*n11+z2*n21)*(-a2*fp+(-a3+d6*(-z1-z2+1))*fpp)+((z3-z1)*n12+(M-z2-z3)*n22)*(a2*fp+(a3+d6*(z1+z2))*fpp)

# Drift in CLT
coeff_fp(z1,z2,z3) = (L0f+L1g+L2h1+L2h2).coefficient(fp)
a(z1) = ((coeff_fp(z1,0,0) - coeff_fp(0,0,0))/z1).full_simplify()
coeff_fp_z1 = a(0)
a(z2) = ((coeff_fp(0,z2,0) - coeff_fp(0,0,0))/z2).full_simplify()
coeff_fp_z2 = a(0)
a(z3) = ((coeff_fp(0,0,z3) - coeff_fp(0,0,0))/z3).full_simplify()
coeff_fp_z3 = a(0)
sol_b = solve([coeff_fp_z1==0, coeff_fp_z2==0, coeff_fp_z3==0], a2, b2,
c2, solution_dict=True)
a2(a1,b1,c1) =  (sol_b[0][a2]).full_simplify()
b2(a1,b1,c1) =  (sol_b[0][b2]).full_simplify()
c2(a1,b1,c1) =  (sol_b[0][c2]).full_simplify()
a2_sol = a2(a1_sol, b1_sol, c1_sol)
b2_sol = b2(a1_sol, b1_sol, c1_sol)
c2_sol = c2(a1_sol, b1_sol, c1_sol)
# print("a2 = ", a2_sol)
# print("b2 = ", b2_sol, "\n")
# print("c2 = ", c2_sol, "\n")
mu_CLT(a1, b1, c1, a2, b2, c2) = coeff_fp(0,0,0).full_simplify()
mu_CLT_sol = mu_CLT(a1_sol, b1_sol, c1_sol, a2_sol, b2_sol, c2_sol)
print("\nmu in CLT: ", mu_CLT(a1_sol, b1_sol, c1_sol, a2_sol, b2_sol,
c2_sol).factor())	
mu_CLT_fromk(k11,k12,k13,k21,k22,k23,n11,n12,n21,n22,n31,n32) = mu_CLT(a1_sol,
b1_sol, c1_sol, a2_sol, b2_sol, c2_sol).full_simplify()
print("...as N->infty: ",
limit(mu_CLT_fromk(k11,k12,k13,k21,k22,k23,N*n11,N*n12,N*n21,N*n22,N*n31,N*n32),
N=oo).factor())

# Diffusion in CLT
coeff_fpp(z1,z2,z3) = (L0f+L1g+L2h1+L2h2).coefficient(fpp)
a(z1) = ((coeff_fpp(z1,0,0) - coeff_fpp(0,0,0))/z1).full_simplify()
coeff_fpp_z1 = a(0)
coeff_fpp_z1z1 = ((coeff_fpp(z1,0,0) - z1*a(0) -
coeff_fpp(0,0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_fpp(0,z2,0) - coeff_fpp(0,0,0))/z2).full_simplify()
coeff_fpp_z2 = a(0)
coeff_fpp_z2z2 = ((coeff_fpp(0,z2,0) - z2*a(0) -
coeff_fpp(0,0,0))/z2/z2).full_simplify()
a(z3) = ((coeff_fpp(0,0,z3) - coeff_fpp(0,0,0))/z3).full_simplify()
coeff_fpp_z3 = a(0)
coeff_fpp_z3z3 = ((coeff_fpp(0,0,z3) - z3*a(0) -
coeff_fpp(0,0,0))/z3/z3).full_simplify()
coeff_fpp_z1z3 = ((coeff_fpp(z1,0,z3) - z1*coeff_fpp_z1 -
z3*coeff_fpp_z3 - z1^2*coeff_fpp_z1z1 - z3^2*coeff_fpp_z3z3 -
coeff_fpp(0,0,0))/z1/z3).full_simplify()
coeff_fpp_z1z2 = ((coeff_fpp(z1,z2,0) - z1*coeff_fpp_z1 -
z2*coeff_fpp_z2 - z1^2*coeff_fpp_z1z1 - z2^2*coeff_fpp_z2z2 -
coeff_fpp(0,0,0))/z1/z2).full_simplify()
coeff_fpp_z2z3 = ((coeff_fpp(0,z2,z3) - z2*coeff_fpp_z2 -
z3*coeff_fpp_z3 - z2^2*coeff_fpp_z2z2 - z3^2*coeff_fpp_z3z3 -
coeff_fpp(0,0,0))/z2/z3).full_simplify()

sol_c = solve([coeff_fpp_z1 == 0, coeff_fpp_z2 == 0, coeff_fpp_z3 == 0,
coeff_fpp_z1z1 == 0, coeff_fpp_z2z2 == 0, coeff_fpp_z3z3 == 0,
coeff_fpp_z1z3 == 0, coeff_fpp_z1z2 == 0, coeff_fpp_z2z3 == 0], a3, b3,
c3, d1, d2, d3, d4, d5, d6, solution_dict=True)
a3(a1,a2,b1,b2,c1,c2) =  (sol_c[0][a3]).full_simplify()
b3(a1,a2,b1,b2,c1,c2) =  (sol_c[0][b3]).full_simplify()
c3(a1,a2,b1,b2,c1,c2) =  (sol_c[0][c3]).full_simplify()
d1(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d1]).full_simplify()
d2(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d2]).full_simplify()
d3(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d3]).full_simplify()
d4(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d4]).full_simplify()
d5(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d5]).full_simplify()
d6(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d6]).full_simplify()
a3_sol=a3(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
b3_sol=b3(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
c3_sol=c3(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d1_sol=d1(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d2_sol=d2(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d3_sol=d3(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d4_sol=d4(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d5_sol=d5(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d6_sol=d6(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)

sigma2_CLT(a1, a2, b1, b2, c1, c2, a3, b3, c3, d1, d2, d3, d4, d5, d6) = coeff_fpp(0,0,0).full_simplify()
sigma2_CLT_sol = sigma2_CLT(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol,
c2_sol, a3_sol, b3_sol, c3_sol, d1_sol, d2_sol, d3_sol, d4_sol, d5_sol,
d6_sol).full_simplify().factor()
print("\n\nsigma2 in CLT: ", sigma2_CLT_sol)
sigma2_CLT_fromk(k11,k12,k13,k21,k22,k23,n11,n12,n21,n22) = sigma2_CLT_sol
print("...as N->infty: ",
limit(sigma2_CLT_fromk(k11,k12,k13,k21,k22,k23,N^-0.5*n11,N^-0.5*n12,N^-0.5*n21,N^-0.5*n22,n31,n32),
N=oo).factor())


print("Section 3.7: Michaelis-Menten with 2 compartments, beta>gamma")

# variables
M,N,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,d5,d6,k11,k12,k13,k21,k22,k23,n11,n12,n21,n22,n31,n32,v1,v2,u2,z1,z2,z3,fp,fp2,fpp = var('M,N,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,d5,d6,k11,k12,k13,k21,k22,k23,n11,n12,n21,n22,n31,n32,v1,v2,u2,z1,z2,z3,fp,fp2,fpp')

k11=k13=k21=k22=k23=n12=n22=n32=n11=n21=n31=1
k12=2
# LLN(beta>gamma)
v1=n32/(n31+n32)*v2
G0f = (-v1*z1*k11+z2*k12-(v2-v1)*(z3-z1)*k21+(M-z2-z3)*k22)*fp
G2g = (-v1*z1*k11+z2*k12+z2*k13-(v2-v1)*(z3-z1)*k21+(M-z2-z3)*k22+(M-z2-z3)*k23)*(a1*fp)

G2h1 = (v1*z1*k11*(-b1+c1)+z2*k12*(b1-c1)+z2*k13*(b1-c1))*fp*0
G2h2 = (z1*n11*(-b1)+(z3-z1)*n12*b1+z2*n21*(-c1)+(M-z2-z3)*n22*c1)*fp

coeff_fp(z1,z2,z3) = (G0f+G2g+G2h1+G2h2).coefficient(fp)
a(z1) = ((coeff_fp(z1,0,0) - coeff_fp(0,0,0))/z1).full_simplify()
coeff_fp_z1 = a(0)
a(z2) = ((coeff_fp(0,z2,0) - coeff_fp(0,0,0))/z2).full_simplify()
coeff_fp_z2 = a(0)
a(z3) = ((coeff_fp(0,0,z3) - coeff_fp(0,0,0))/z3).full_simplify()
coeff_fp_z3 = a(0)


sol_a = solve([coeff_fp_z1==0, coeff_fp_z2==0, coeff_fp_z3==0], a1, b1,
c1, solution_dict=True)
a1_sol =  (sol_a[0][a1]).full_simplify()
b1_sol =  (sol_a[0][b1]).full_simplify()
c1_sol =  (sol_a[0][c1]).full_simplify()
mu_LLN(a1, b1, c1) = coeff_fp(0,0,0).full_simplify()
mu_LLN_sol = mu_LLN(a1_sol, b1_sol, c1_sol)
print("\nmu in LLN: ", mu_LLN(a1_sol, b1_sol, c1_sol).factor())	
mu_LLN_fromk(k11,k12,k13,k21,k22,k23,n11,n12,n21,n22,n31,n32) = mu_LLN(a1_sol, b1_sol, c1_sol)
print("...as N->infty: ",
limit(mu_LLN_fromk(k11,k12,k13,k21,k22,k23,N*n11,N*n12,N*n21,N*n22,N*n31,N*n32),
N=oo).full_simplify())

k11=k13=k21=k22=k23=n12=n22=n32=n11=n21=n31=1
k12=2
# CLT(beta>gamma)
v1=n32/(n31+n32)*v2
L0f = 1/2*(v1*z1*k11+z2*k12+(v2-v1)*(z3-z1)*k21+(M-z2-z3)*k22)*fpp-u2*n31/(n31+n32)*(z3-z1)*k21*fp-u2*z1*k11*n32/(n31+n32)*fp
L1g = (-v1*z1*k11*(a1*(z3-1))+z2*k12*(a1*(z3+1))-(v2-v1)*(z3-z1)*k21*(a1*(z3-1))+(M-z2-z3)*k22*(a1*(z3+1)))*fpp-u2*n31/(n31+n32)*(z3-z1)*k21*a1*fp+u2*z1*k11*n32/(n31+n32)*(-a1)*fp-mu_LLN(a1,b1,c1)*(a1*z3)*fpp
L2h1 = v1*z1*k11*(-c2*fp+(-c3+d3*(-z3+1))*fpp)+z2*k12*(c2*fp+(c3+d3*z3)*fpp)+z2*k13*(c2*fp+(c3+d3*z3)*fpp)+(v2-v1)*(z3-z1)*k21*(-c2*fp+(-c3+d3*(-z3+1))*fpp)+(M-z2-z3)*k22*(c2*fp+(c3+d3*z3)*fpp)+(M-z2-z3)*k23*(c2*fp+(c3+d3*z3)*fpp)
L2h2 = z1*n11*(-a2*fp+(-a3+d1*(-z1+1)+d4*(-z1-z2+1)+d6*(-z1-z3+1))*fpp)+(z3-z1)*n12*(a2*fp+(a3+d1*z1+d4*(z1+z2)+d6*(z1+z3))*fpp)+z2*n21*(-b2*fp+(-b3+d2*(-z2+1)+d4*(-z1-z2+1)+d5*(-z2-z3+1))*fpp)+(M-z2-z3)*n22*(b2*fp+(b3+d2*z2+d4*(z1+z2)+d5*(z2+z3))*fpp)

# Drift in CLT
coeff_fp(z1,z2,z3) = (L0f+L1g+L2h1+L2h2).coefficient(fp)
a(z1) = ((coeff_fp(z1,0,0) - coeff_fp(0,0,0))/z1).full_simplify()
coeff_fp_z1 = a(0)
a(z2) = ((coeff_fp(0,z2,0) - coeff_fp(0,0,0))/z2).full_simplify()
coeff_fp_z2 = a(0)
a(z3) = ((coeff_fp(0,0,z3) - coeff_fp(0,0,0))/z3).full_simplify()
coeff_fp_z3 = a(0)
sol_b = solve([coeff_fp_z1==0, coeff_fp_z2==0, coeff_fp_z3==0], a2, b2,
c2, solution_dict=True)
a2(a1,b1,c1) =  (sol_b[0][a2]).full_simplify()
b2(a1,b1,c1) =  (sol_b[0][b2]).full_simplify()
c2(a1,b1,c1) =  (sol_b[0][c2]).full_simplify()
a2_sol = a2(a1_sol, b1_sol, c1_sol)
b2_sol = b2(a1_sol, b1_sol, c1_sol)
c2_sol = c2(a1_sol, b1_sol, c1_sol)
# print("a2 = ", a2_sol)
# print("b2 = ", b2_sol, "\n")
# print("c2 = ", c2_sol, "\n")
mu_CLT(a1, b1, c1, a2, b2, c2) = coeff_fp(0,0,0).full_simplify()
mu_CLT_sol = mu_CLT(a1_sol, b1_sol, c1_sol, a2_sol, b2_sol, c2_sol)
print("\nmu in CLT: ", mu_CLT(a1_sol, b1_sol, c1_sol, a2_sol, b2_sol,
c2_sol).factor())	
mu_CLT_fromk(k11,k12,k13,k21,k22,k23,n11,n12,n21,n22,n31,n32) = mu_CLT(a1_sol,
b1_sol, c1_sol, a2_sol, b2_sol, c2_sol).full_simplify()
print("...as N->infty: ",
limit(mu_CLT_fromk(k11,k12,k13,k21,k22,k23,N*n11,N*n12,N*n21,N*n22,N*n31,N*n32),
N=oo).factor())

# Diffusion in CLT
coeff_fpp(z1,z2,z3) = (L0f+L1g+L2h1+L2h2).coefficient(fpp)
a(z1) = ((coeff_fpp(z1,0,0) - coeff_fpp(0,0,0))/z1).full_simplify()
coeff_fpp_z1 = a(0)
coeff_fpp_z1z1 = ((coeff_fpp(z1,0,0) - z1*a(0) -
coeff_fpp(0,0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_fpp(0,z2,0) - coeff_fpp(0,0,0))/z2).full_simplify()
coeff_fpp_z2 = a(0)
coeff_fpp_z2z2 = ((coeff_fpp(0,z2,0) - z2*a(0) -
coeff_fpp(0,0,0))/z2/z2).full_simplify()
a(z3) = ((coeff_fpp(0,0,z3) - coeff_fpp(0,0,0))/z3).full_simplify()
coeff_fpp_z3 = a(0)
coeff_fpp_z3z3 = ((coeff_fpp(0,0,z3) - z3*a(0) -
coeff_fpp(0,0,0))/z3/z3).full_simplify()
coeff_fpp_z1z3 = ((coeff_fpp(z1,0,z3) - z1*coeff_fpp_z1 -
z3*coeff_fpp_z3 - z1^2*coeff_fpp_z1z1 - z3^2*coeff_fpp_z3z3 -
coeff_fpp(0,0,0))/z1/z3).full_simplify()
coeff_fpp_z1z2 = ((coeff_fpp(z1,z2,0) - z1*coeff_fpp_z1 -
z2*coeff_fpp_z2 - z1^2*coeff_fpp_z1z1 - z2^2*coeff_fpp_z2z2 -
coeff_fpp(0,0,0))/z1/z2).full_simplify()
coeff_fpp_z2z3 = ((coeff_fpp(0,z2,z3) - z2*coeff_fpp_z2 -
z3*coeff_fpp_z3 - z2^2*coeff_fpp_z2z2 - z3^2*coeff_fpp_z3z3 -
coeff_fpp(0,0,0))/z2/z3).full_simplify()

sol_c = solve([coeff_fpp_z1 == 0, coeff_fpp_z2 == 0, coeff_fpp_z3 == 0,
coeff_fpp_z1z1 == 0, coeff_fpp_z2z2 == 0, coeff_fpp_z3z3 == 0,
coeff_fpp_z1z3 == 0, coeff_fpp_z1z2 == 0, coeff_fpp_z2z3 == 0], a3, b3,
c3, d1, d2, d3, d4, d5, d6, solution_dict=True)
a3(a1,a2,b1,b2,c1,c2) =  (sol_c[0][a3]).full_simplify()
b3(a1,a2,b1,b2,c1,c2) =  (sol_c[0][b3]).full_simplify()
c3(a1,a2,b1,b2,c1,c2) =  (sol_c[0][c3]).full_simplify()
d1(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d1]).full_simplify()
d2(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d2]).full_simplify()
d3(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d3]).full_simplify()
d4(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d4]).full_simplify()
d5(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d5]).full_simplify()
d6(a1,a2,b1,b2,c1,c2) =  (sol_c[0][d6]).full_simplify()
a3_sol=a3(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
b3_sol=b3(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
c3_sol=c3(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d1_sol=d1(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d2_sol=d2(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d3_sol=d3(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d4_sol=d4(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d5_sol=d5(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)
d6_sol=d6(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol, c2_sol)

sigma2_CLT(a1, a2, b1, b2, c1, c2, a3, b3, c3, d1, d2, d3, d4, d5, d6) = coeff_fpp(0,0,0).full_simplify()
sigma2_CLT_sol = sigma2_CLT(a1_sol, a2_sol, b1_sol, b2_sol, c1_sol,
c2_sol, a3_sol, b3_sol, c3_sol, d1_sol, d2_sol, d3_sol, d4_sol, d5_sol,
d6_sol).full_simplify().factor()
print("\n\nsigma2 in CLT: ", sigma2_CLT_sol)
sigma2_CLT_fromk(k11,k12,k13,k21,k22,k23,n11,n12,n21,n22) = sigma2_CLT_sol
print("...as N->infty: ",
limit(sigma2_CLT_fromk(k11,k12,k13,k21,k22,k23,N*n11,N*n12,N*n21,N*n22),
N=oo).factor())


print("Section 3.8: Signalling proteins with 2 compartments, B fast")

# variables
N,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26,d27,d28,d29,d30,k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42,v2,v4,u2,u4,z1,z2,z3,z4,D2f,D4f,D22f,D44f,D24f = var('N,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26,d27,d28,d29,d30,k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42,v2,v4,u2,u4,z1,z2,z3,z4,D2f,D4f,D22f,D44f,D24f')
k11=k13=k14=k15=k16=k17=k18=k21=k22=k23=k24=k25=k26=k27=k28=n12=n22=n32=n11=n21=n31=n41=n42=1
k12=2
# LLN
G0f = (z1*k13-v2*n22/(n22+n21)*k15+(z2-z1)*k23-v2*n21/(n22+n21)*k25)*D2f+(z3*v2*n22/(n21+n22)*k17-v4*n42/(n41+n42)*k18+v2*n21/(n22+n21)*(z4-z3)*k27-v4*n41/(n41+n42)*k28)*D4f
g1(z2,z4)=(a1*z2+a2*z4)*D2f+(a3*z2+a4*z4)*D4f
G1g = k11*(g1(z2+1,z4)-g1(z2,z4))+z1*k12*(g1(z2-1,z4)-g1(z2,z4))+z1*k14*(g1(z2,z4+1)-g1(z2,z4))+z3*k16*(g1(z2,z4-1)-g1(z2,z4))+k21*(g1(z2+1,z4)-g1(z2,z4))+(z2-z1)*k22*(g1(z2-1,z4)-g1(z2,z4))+(z2-z1)*k24*(g1(z2,z4+1)-g1(z2,z4))+(z4-z3)*k26*(g1(z2,z4-1)-g1(z2,z4))
g2(z1,z3)=(a5*z1+a6*z3)*D2f+(a7*z1+a8*z3)*D4f
G2g = z1*n11*(g2(z1-1,z3)-g2(z1,z3))+(z2-z1)*n12*(g2(z1+1,z3)-g2(z1,z3))+z3*n31*(g2(z1,z3-1)-g2(z1,z3))+(z4-z3)*n32*(g2(z1,z3+1)-g2(z1,z3))

coeff_LLN_D2f(z1,z2,z3,z4) = (G0f+G1g+G2g).coefficient(D2f)
a(z1) = ((coeff_LLN_D2f(z1,0,0,0) - coeff_LLN_D2f(0,0,0,0))/z1).full_simplify()
coeff_LLN_D2f_z1 = a(0)
a(z2) = ((coeff_LLN_D2f(0,z2,0,0) - coeff_LLN_D2f(0,0,0,0))/z2).full_simplify()
coeff_LLN_D2f_z2 = a(0)
a(z3) = ((coeff_LLN_D2f(0,0,z3,0) - coeff_LLN_D2f(0,0,0,0))/z3).full_simplify()
coeff_LLN_D2f_z3 = a(0)
a(z2) = ((coeff_LLN_D2f(0,0,0,z4) - coeff_LLN_D2f(0,0,0,0))/z4).full_simplify()
coeff_LLN_D2f_z4 = a(0)
coeff_LLN_D4f(z1,z2,z3,z4) = (G0f+G1g+G2g).coefficient(D4f)
a(z1) = ((coeff_LLN_D4f(z1,0,0,0) - coeff_LLN_D4f(0,0,0,0))/z1).full_simplify()
coeff_LLN_D4f_z1 = a(0)
a(z2) = ((coeff_LLN_D4f(0,z2,0,0) - coeff_LLN_D4f(0,0,0,0))/z2).full_simplify()
coeff_LLN_D4f_z2 = a(0)
a(z3) = ((coeff_LLN_D4f(0,0,z3,0) - coeff_LLN_D4f(0,0,0,0))/z3).full_simplify()
coeff_LLN_D4f_z3 = a(0)
a(z2) = ((coeff_LLN_D4f(0,0,0,z4) - coeff_LLN_D4f(0,0,0,0))/z4).full_simplify()
coeff_LLN_D4f_z4 = a(0)


sol_ab = solve([coeff_LLN_D2f_z1==0, coeff_LLN_D2f_z2==0, coeff_LLN_D2f_z3==0, coeff_LLN_D2f_z4==0, coeff_LLN_D4f_z1==0, coeff_LLN_D4f_z2==0, coeff_LLN_D4f_z3==0, coeff_LLN_D4f_z4==0], a1, a2, a3, a4, a5, a6, a7, a8, solution_dict=True)
a1_sol =  (sol_ab[0][a1]).full_simplify()
a2_sol =  (sol_ab[0][a2]).full_simplify()
a3_sol =  (sol_ab[0][a3]).full_simplify()
a4_sol =  (sol_ab[0][a4]).full_simplify()
a5_sol =  (sol_ab[0][a5]).full_simplify()
a6_sol =  (sol_ab[0][a6]).full_simplify()
a7_sol =  (sol_ab[0][a7]).full_simplify()
a8_sol =  (sol_ab[0][a8]).full_simplify()


mu_LLN_2(a1, a2, a3, a4, a5, a6, a7, a8) = coeff_LLN_D2f(0,0).full_simplify()
mu_LLN_2_sol = mu_LLN_2(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol)
mu_LLN_4(a1, a2, a3, a4, a5, a6, a7, a8) = coeff_LLN_D4f(0,0).full_simplify()
mu_LLN_4_sol = mu_LLN_4(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol)
mu_LLN_2_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = mu_LLN_2_sol
mu_LLN_4_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = mu_LLN_4_sol

print("LLN:")
print("coeff of D2f as N->infty: ", limit(mu_LLN_2_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo).full_simplify())
print("coeff of D4f as N->infty: ", limit(mu_LLN_4_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo).full_simplify())

# CLT
k11=k13=k14=k15=k16=k17=k18=k21=k22=k23=k24=k25=k26=k27=k28=n12=n22=n32=n11=n21=n31=n41=n42=1
k12=2
L0f = 1/2*(z1*k13+v2*n22/(n21+n22)*k15+(z2-z1)*k23+v2*n21/(n21+n22)*k25)*D22f+1/2*(v2*n22/(n22+n21)*k17*z3+v4*n42/(n42+n41)*k18+v2*n21/(n22+n21)*k27*(z4-z3)+v4*n41/(n42+n41)*k28)*D44f+(-u2*n22/(n22+n21)*k15-u2*n21/(n22+n21)*k25)*D2f+(u2*n22/(n21+n22)*k17*z3-u4*n42/(n41+n42)*k18+u2*n21/(n21+n22)*k27*(z4-z3)-u4*n41/(n41+n42)*k28)*D4f
D2g(z2,z4) = (a1*z2+a2*z4)*D22f + (a3*z2+a4*z4)*D24f
D4g(z2,z4) = (a1*z2+a2*z4)*D24f + (a3*z2+a4*z4)*D44f
L1g = (z1*k13-v2*n22/(n21+n22)*k15+(z2-z1)*k23-v2*n21/(n22+n21)*k25)*D2g(z2,z4)+(v2*n22/(n21+n22)*k17*z3-v4*n42/(n41+n42)*k18+v2*n21/(n21+n22)*k27*(z4-z3)-v4*n41/(n41+n42)*k28)*D4g(z2,z4) - (-v2+4/3)*D2g(z2,z4)-(2/3*v2-v4)*D4g(z2,z4)

h1(z2,z4)  = (b1*z2+b2*z4)*D2f + (b3*z2+b4*z4)*D4f + (c1*z2+c2*z4+d1*z2*(z2-1)/2 + d2*z4*(z4-1)/2 + d3*(z2+z4)*(z2+z4-1)/2) * D22f + (c3*z2+c4*z4+d4*z2*(z2-1)/2 + d5*z4*(z4-1)/2 + d6*(z2+z4)*(z2+z4-1)/2) * D44f + (c5*z2+c6*z4+d7*z2*(z2-1)/2 + d8*z4*(z4-1)/2 + d9*(z2+z4)*(z2+z4-1)/2) * D24f

L2h1 = k11*(h1(z2+1,z4)-h1(z2,z4))+z1*k12*(h1(z2-1,z4)-h1(z2,z4))+z1*k14*(h1(z2,z4+1)-h1(z2,z4))+z3*k16*(h1(z2,z4-1)-h1(z2,z4))+k21*(h1(z2+1,z4)-h1(z2,z4))+(z2-z1)*k22*(h1(z2-1,z4)-h1(z2,z4))+(z2-z1)*k24*(h1(z2,z4+1)-h1(z2,z4))+(z4-z3)*k26*(h1(z2,z4-1)-h1(z2,z4))

h2(z1,z3,z2,z4)  = (b5*z1+b6*z3)*D2f + (b7*z1+b8*z3)*D4f + (c7*z1+c8*z3+d10*z1*(z1-1)/2 + d11*z3*(z3-1)/2+d12*(z1+z2)*(z1+z2-1)/2+d13*(z1+z3)*(z1+z3-1)/2+d14*(z1+z4)*(z1+z4-1)/2+d15*(z2+z3)*(z2+z3-1)/2+d16*(z3+z4)*(z3+z4-1)/2) * D22f + (c9*z1+c10*z3+d17*z1*(z1-1)/2 + d18*z3*(z3-1)/2+d19*(z1+z2)*(z1+z2-1)/2+d20*(z1+z3)*(z1+z3-1)/2+d21*(z1+z4)*(z1+z4-1)/2+d22*(z2+z3)*(z2+z3-1)/2+d23*(z3+z4)*(z3+z4-1)/2) * D44f + (c11*z1+c12*z3+d24*z1*(z1-1)/2 + d25*z3*(z3-1)/2+d26*(z1+z2)*(z1+z2-1)/2+d27*(z1+z3)*(z1+z3-1)/2+d28*(z1+z4)*(z1+z4-1)/2+d29*(z2+z3)*(z2+z3-1)/2+d30*(z3+z4)*(z3+z4-1)/2) * D24f

L2h2 = z1*n11*(h2(z1-1,z3,z2,z4)-h2(z1,z3,z2,z4))+(z2-z1)*n12*(h2(z1+1,z3,z2,z4)-h2(z1,z3,z2,z4))+z3*n31*(h2(z1,z3-1,z2,z4)-h2(z1,z3,z2,z4))+(z4-z3)*n32*(h2(z1,z3+1,z2,z4)-h2(z1,z3,z2,z4))

# Drift in CLT
coeff_CLT_D2f(z1,z2,z3,z4) = (L0f+L1g+L2h1+L2h2).coefficient(D2f)
coeff_CLT_D4f(z1,z2,z3,z4) = (L0f+L1g+L2h1+L2h2).coefficient(D4f)
a(z1) = ((coeff_CLT_D2f(z1,0,0,0) - coeff_CLT_D2f(0,0,0,0))/z1).full_simplify()
coeff_CLT_D2f_z1 = a(0)
a(z2) = ((coeff_CLT_D2f(0,z2,0,0) - coeff_CLT_D2f(0,0,0,0))/z2).full_simplify()
coeff_CLT_D2f_z2 = a(0)
a(z3) = ((coeff_CLT_D2f(0,0,z3,0) - coeff_CLT_D2f(0,0,0,0))/z3).full_simplify()
coeff_CLT_D2f_z3 = a(0)
a(z4) = ((coeff_CLT_D2f(0,0,0,z4) - coeff_CLT_D2f(0,0,0,0))/z4).full_simplify()
coeff_CLT_D2f_z4 = a(0)
a(z1) = ((coeff_CLT_D4f(z1,0,0,0) - coeff_CLT_D4f(0,0,0,0))/z1).full_simplify()
coeff_CLT_D4f_z1 = a(0)
a(z2) = ((coeff_CLT_D4f(0,z2,0,0) - coeff_CLT_D4f(0,0,0,0))/z2).full_simplify()
coeff_CLT_D4f_z2 = a(0)
a(z3) = ((coeff_CLT_D4f(0,0,z3,0) - coeff_CLT_D4f(0,0,0,0))/z3).full_simplify()
coeff_CLT_D4f_z3 = a(0)
a(z4) = ((coeff_CLT_D4f(0,0,0,z4) - coeff_CLT_D4f(0,0,0,0))/z4).full_simplify()
coeff_CLT_D4f_z4 = a(0)

sol_b = solve([coeff_CLT_D2f_z1==0, coeff_CLT_D2f_z2==0, coeff_CLT_D2f_z3==0, coeff_CLT_D2f_z4==0, coeff_CLT_D4f_z1==0, coeff_CLT_D4f_z2==0, coeff_CLT_D4f_z3==0, coeff_CLT_D4f_z4==0], b1, b2, b3, b4, b5, b6, b7, b8, solution_dict=True)
b1(a1,a2,a3,a4,a5,a6,a7,a8) =  (sol_b[0][b1]).full_simplify()
b2(a1,a2,a3,a4,a5,a6,a7,a8) =  (sol_b[0][b2]).full_simplify()
b3(a1,a2,a3,a4,a5,a6,a7,a8) =  (sol_b[0][b3]).full_simplify()
b4(a1,a2,a3,a4,a5,a6,a7,a8) =  (sol_b[0][b4]).full_simplify()
b5(a1,a2,a3,a4,a5,a6,a7,a8) =  (sol_b[0][b5]).full_simplify()
b6(a1,a2,a3,a4,a5,a6,a7,a8) =  (sol_b[0][b6]).full_simplify()
b7(a1,a2,a3,a4,a5,a6,a7,a8) =  (sol_b[0][b7]).full_simplify()
b8(a1,a2,a3,a4,a5,a6,a7,a8) =  (sol_b[0][b8]).full_simplify()
b1_sol = b1(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol).full_simplify()
b2_sol = b2(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol).full_simplify()
b3_sol = b3(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol).full_simplify()
b4_sol = b4(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol).full_simplify()
b5_sol = b5(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol).full_simplify()
b6_sol = b6(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol).full_simplify()
b7_sol = b7(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol).full_simplify()
b8_sol = b8(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol).full_simplify()


mu_CLT_D2f(a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8) = coeff_CLT_D2f(0,0).full_simplify()
mu_CLT_D2f_sol = mu_CLT_D2f(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol)
mu_CLT_D4f(a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8) = coeff_CLT_D4f(0,0).full_simplify()
mu_CLT_D4f_sol = mu_CLT_D4f(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol)
mu_CLT_D2f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = mu_CLT_D2f_sol
mu_CLT_D4f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = mu_CLT_D4f_sol

print("CLT:")
print("coeff of D1f as N->infty: ", limit(mu_CLT_D2f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))
print("coeff of D3f as N->infty: ", limit(mu_CLT_D4f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))


# Diffusiont in CLT
coeff_CLT_D22f(z1,z2,z3,z4) = (L0f+L1g+L2h1+L2h2).coefficient(D22f).full_simplify()
coeff_CLT_D44f(z1,z2,z3,z4) = (L0f+L1g+L2h1+L2h2).coefficient(D44f).full_simplify()
coeff_CLT_D24f(z1,z2,z3,z4) = (L0f+L1g+L2h1+L2h2).coefficient(D24f).full_simplify()


a(z1) = ((coeff_CLT_D22f(z1,0,0,0) - coeff_CLT_D22f(0,0,0,0))/z1).full_simplify()
coeff_CLT_D22f_z1 = a(0).full_simplify()
coeff_CLT_D22f_z1z1 = ((coeff_CLT_D22f(z1,0,0,0) - z1*a(0) - coeff_CLT_D22f(0,0,0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D22f(0,z2,0,0) - coeff_CLT_D22f(0,0,0,0))/z2).full_simplify()
coeff_CLT_D22f_z2 = a(0).full_simplify()
coeff_CLT_D22f_z2z2 = ((coeff_CLT_D22f(0,z2,0,0) - z2*a(0) - coeff_CLT_D22f(0,0,0,0))/z2/z2).full_simplify()
coeff_CLT_D22f_z1z2 = ((coeff_CLT_D22f(z1,z2,0,0) - z1*coeff_CLT_D22f_z1 - z2*coeff_CLT_D22f_z2 - z1^2*coeff_CLT_D22f_z1z1 - z2^2*coeff_CLT_D22f_z2z2 - coeff_CLT_D22f(0,0,0,0))/z1/z2).full_simplify()
a(z3) = ((coeff_CLT_D22f(0,0,z3,0) - coeff_CLT_D22f(0,0,0,0))/z3).full_simplify()
coeff_CLT_D22f_z3 = a(0).full_simplify()
coeff_CLT_D22f_z3z3 = ((coeff_CLT_D22f(0,0,z3,0) - z3*a(0) - coeff_CLT_D22f(0,0,0,0))/z3/z3).full_simplify()
coeff_CLT_D22f_z1z3 = ((coeff_CLT_D22f(z1,0,z3,0) - z1*coeff_CLT_D22f_z1 - z3*coeff_CLT_D22f_z3 - z1^2*coeff_CLT_D22f_z1z1 - z3^2*coeff_CLT_D22f_z3z3 - coeff_CLT_D22f(0,0,0,0))/z1/z3).full_simplify()
coeff_CLT_D22f_z2z3 = ((coeff_CLT_D22f(0,z2,z3,0) - z2*coeff_CLT_D22f_z2 - z3*coeff_CLT_D22f_z3 - z2^2*coeff_CLT_D22f_z2z2 - z3^2*coeff_CLT_D22f_z3z3 - coeff_CLT_D22f(0,0,0,0))/z2/z3).full_simplify()
a(z4) = ((coeff_CLT_D22f(0,0,0,z4) - coeff_CLT_D22f(0,0,0,0))/z4).full_simplify()
coeff_CLT_D22f_z4 = a(0).full_simplify()
coeff_CLT_D22f_z4z4 = ((coeff_CLT_D22f(0,0,0,z4) - z4*a(0) - coeff_CLT_D22f(0,0,0,0))/z4/z4).full_simplify()
coeff_CLT_D22f_z1z4 = ((coeff_CLT_D22f(z1,0,0,z4) - z1*coeff_CLT_D22f_z1 - z4*coeff_CLT_D22f_z4 - z1^2*coeff_CLT_D22f_z1z1 - z4^2*coeff_CLT_D22f_z4z4 - coeff_CLT_D22f(0,0,0,0))/z1/z4).full_simplify()
coeff_CLT_D22f_z2z4 = ((coeff_CLT_D22f(0,z2,0,z4) - z2*coeff_CLT_D22f_z2 - z4*coeff_CLT_D22f_z4 - z2^2*coeff_CLT_D22f_z2z2 - z4^2*coeff_CLT_D22f_z4z4 - coeff_CLT_D22f(0,0,0,0))/z2/z4).full_simplify()
coeff_CLT_D22f_z3z4 = ((coeff_CLT_D22f(0,0,z3,z4) - z3*coeff_CLT_D22f_z3 - z4*coeff_CLT_D22f_z4 - z3^2*coeff_CLT_D22f_z3z3 - z4^2*coeff_CLT_D22f_z4z4 - coeff_CLT_D22f(0,0,0,0))/z3/z4).full_simplify()

a(z1) = ((coeff_CLT_D44f(z1,0,0,0) - coeff_CLT_D44f(0,0,0,0))/z1).full_simplify()
coeff_CLT_D44f_z1 = a(0).full_simplify()
coeff_CLT_D44f_z1z1 = ((coeff_CLT_D44f(z1,0,0,0) - z1*a(0) - coeff_CLT_D44f(0,0,0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D44f(0,z2,0,0) - coeff_CLT_D44f(0,0,0,0))/z2).full_simplify()
coeff_CLT_D44f_z2 = a(0).full_simplify()
coeff_CLT_D44f_z2z2 = ((coeff_CLT_D44f(0,z2,0,0) - z2*a(0) - coeff_CLT_D44f(0,0,0,0))/z2/z2).full_simplify()
coeff_CLT_D44f_z1z2 = ((coeff_CLT_D44f(z1,z2,0,0) - z1*coeff_CLT_D44f_z1 - z2*coeff_CLT_D44f_z2 - z1^2*coeff_CLT_D44f_z1z1 - z2^2*coeff_CLT_D44f_z2z2 - coeff_CLT_D44f(0,0,0,0))/z1/z2).full_simplify()
a(z3) = ((coeff_CLT_D44f(0,0,z3,0) - coeff_CLT_D44f(0,0,0,0))/z3).full_simplify()
coeff_CLT_D44f_z3 = a(0).full_simplify()
coeff_CLT_D44f_z3z3 = ((coeff_CLT_D44f(0,0,z3,0) - z3*a(0) - coeff_CLT_D44f(0,0,0,0))/z3/z3).full_simplify()
coeff_CLT_D44f_z1z3 = ((coeff_CLT_D44f(z1,0,z3,0) - z1*coeff_CLT_D44f_z1 - z3*coeff_CLT_D44f_z3 - z1^2*coeff_CLT_D44f_z1z1 - z3^2*coeff_CLT_D44f_z3z3 - coeff_CLT_D44f(0,0,0,0))/z1/z3).full_simplify()
coeff_CLT_D44f_z2z3 = ((coeff_CLT_D44f(0,z2,z3,0) - z2*coeff_CLT_D44f_z2 - z3*coeff_CLT_D44f_z3 - z2^2*coeff_CLT_D44f_z2z2 - z3^2*coeff_CLT_D44f_z3z3 - coeff_CLT_D44f(0,0,0,0))/z2/z3).full_simplify()
a(z4) = ((coeff_CLT_D44f(0,0,0,z4) - coeff_CLT_D44f(0,0,0,0))/z4).full_simplify()
coeff_CLT_D44f_z4 = a(0).full_simplify()
coeff_CLT_D44f_z4z4 = ((coeff_CLT_D44f(0,0,0,z4) - z4*a(0) - coeff_CLT_D44f(0,0,0,0))/z4/z4).full_simplify()
coeff_CLT_D44f_z1z4 = ((coeff_CLT_D44f(z1,0,0,z4) - z1*coeff_CLT_D44f_z1 - z4*coeff_CLT_D44f_z4 - z1^2*coeff_CLT_D44f_z1z1 - z4^2*coeff_CLT_D44f_z4z4 - coeff_CLT_D44f(0,0,0,0))/z1/z4).full_simplify()
coeff_CLT_D44f_z2z4 = ((coeff_CLT_D44f(0,z2,0,z4) - z2*coeff_CLT_D44f_z2 - z4*coeff_CLT_D44f_z4 - z2^2*coeff_CLT_D44f_z2z2 - z4^2*coeff_CLT_D44f_z4z4 - coeff_CLT_D44f(0,0,0,0))/z2/z4).full_simplify()
coeff_CLT_D44f_z3z4 = ((coeff_CLT_D44f(0,0,z3,z4) - z3*coeff_CLT_D44f_z3 - z4*coeff_CLT_D44f_z4 - z3^2*coeff_CLT_D44f_z3z3 - z4^2*coeff_CLT_D44f_z4z4 - coeff_CLT_D44f(0,0,0,0))/z3/z4).full_simplify()

a(z1) = ((coeff_CLT_D24f(z1,0,0,0) - coeff_CLT_D24f(0,0,0,0))/z1).full_simplify()
coeff_CLT_D24f_z1 = a(0).full_simplify()
coeff_CLT_D24f_z1z1 = ((coeff_CLT_D24f(z1,0,0,0) - z1*a(0) - coeff_CLT_D24f(0,0,0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D24f(0,z2,0,0) - coeff_CLT_D24f(0,0,0,0))/z2).full_simplify()
coeff_CLT_D24f_z2 = a(0).full_simplify()
coeff_CLT_D24f_z2z2 = ((coeff_CLT_D24f(0,z2,0,0) - z2*a(0) - coeff_CLT_D24f(0,0,0,0))/z2/z2).full_simplify()
coeff_CLT_D24f_z1z2 = ((coeff_CLT_D24f(z1,z2,0,0) - z1*coeff_CLT_D24f_z1 - z2*coeff_CLT_D24f_z2 - z1^2*coeff_CLT_D24f_z1z1 - z2^2*coeff_CLT_D24f_z2z2 - coeff_CLT_D24f(0,0,0,0))/z1/z2).full_simplify()
a(z3) = ((coeff_CLT_D24f(0,0,z3,0) - coeff_CLT_D24f(0,0,0,0))/z3).full_simplify()
coeff_CLT_D24f_z3 = a(0).full_simplify()
coeff_CLT_D24f_z3z3 = ((coeff_CLT_D24f(0,0,z3,0) - z3*a(0) - coeff_CLT_D24f(0,0,0,0))/z3/z3).full_simplify()
coeff_CLT_D24f_z1z3 = ((coeff_CLT_D24f(z1,0,z3,0) - z1*coeff_CLT_D24f_z1 - z3*coeff_CLT_D24f_z3 - z1^2*coeff_CLT_D24f_z1z1 - z3^2*coeff_CLT_D24f_z3z3 - coeff_CLT_D24f(0,0,0,0))/z1/z3).full_simplify()
coeff_CLT_D24f_z2z3 = ((coeff_CLT_D24f(0,z2,z3,0) - z2*coeff_CLT_D24f_z2 - z3*coeff_CLT_D24f_z3 - z2^2*coeff_CLT_D24f_z2z2 - z3^2*coeff_CLT_D24f_z3z3 - coeff_CLT_D24f(0,0,0,0))/z2/z3).full_simplify()
a(z4) = ((coeff_CLT_D24f(0,0,0,z4) - coeff_CLT_D24f(0,0,0,0))/z4).full_simplify()
coeff_CLT_D24f_z4 = a(0).full_simplify()
coeff_CLT_D24f_z4z4 = ((coeff_CLT_D24f(0,0,0,z4) - z4*a(0) - coeff_CLT_D24f(0,0,0,0))/z4/z4).full_simplify()
coeff_CLT_D24f_z1z4 = ((coeff_CLT_D24f(z1,0,0,z4) - z1*coeff_CLT_D24f_z1 - z4*coeff_CLT_D24f_z4 - z1^2*coeff_CLT_D24f_z1z1 - z4^2*coeff_CLT_D24f_z4z4 - coeff_CLT_D24f(0,0,0,0))/z1/z4).full_simplify()
coeff_CLT_D24f_z2z4 = ((coeff_CLT_D24f(0,z2,0,z4) - z2*coeff_CLT_D24f_z2 - z4*coeff_CLT_D24f_z4 - z2^2*coeff_CLT_D24f_z2z2 - z4^2*coeff_CLT_D24f_z4z4 - coeff_CLT_D24f(0,0,0,0))/z2/z4).full_simplify()
coeff_CLT_D24f_z3z4 = ((coeff_CLT_D24f(0,0,z3,z4) - z3*coeff_CLT_D24f_z3 - z4*coeff_CLT_D24f_z4 - z3^2*coeff_CLT_D24f_z3z3 - z4^2*coeff_CLT_D24f_z4z4 - coeff_CLT_D24f(0,0,0,0))/z3/z4).full_simplify()



sol_cA = solve([coeff_CLT_D22f_z1==0, coeff_CLT_D22f_z2==0, coeff_CLT_D22f_z3==0, coeff_CLT_D22f_z4==0, coeff_CLT_D22f_z1z1==0, coeff_CLT_D22f_z2z2==0, coeff_CLT_D22f_z3z3==0, coeff_CLT_D22f_z4z4==0, coeff_CLT_D22f_z1z2==0, coeff_CLT_D22f_z1z3==0, coeff_CLT_D22f_z2z3==0, coeff_CLT_D22f_z1z4==0, coeff_CLT_D22f_z2z4==0, coeff_CLT_D22f_z3z4==0], c1, c2, c7, c8, d1, d2, d3, d10, d11, d12, d13, d14, d15, d16, solution_dict=True)
c1(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][c1]).full_simplify()
c2(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][c2]).full_simplify()
c7(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][c7]).full_simplify()
c8(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][c8]).full_simplify()
d1(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][d1]).full_simplify()
d2(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][d2]).full_simplify()
d3(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][d3]).full_simplify()
d10(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][d10]).full_simplify()
d11(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][d11]).full_simplify()
d12(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][d12]).full_simplify()
d13(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][d13]).full_simplify()
d14(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][d14]).full_simplify()
d15(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][d15]).full_simplify()
d16(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cA[0][d16]).full_simplify()
c1_sol = c1(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
c2_sol = c2(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
c7_sol = c7(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
c8_sol = c8(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d1_sol = d1(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d2_sol = d2(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d3_sol = d3(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d10_sol = d10(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d11_sol = d11(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d12_sol = d12(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d13_sol = d13(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d14_sol = d14(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d15_sol = d15(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d16_sol = d16(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()



sol_cB = solve([coeff_CLT_D44f_z1==0, coeff_CLT_D44f_z2==0, coeff_CLT_D44f_z3==0, coeff_CLT_D44f_z4==0, coeff_CLT_D44f_z1z1==0, coeff_CLT_D44f_z2z2==0, coeff_CLT_D44f_z3z3==0, coeff_CLT_D44f_z4z4==0, coeff_CLT_D44f_z1z2==0, coeff_CLT_D44f_z1z3==0, coeff_CLT_D44f_z2z3==0, coeff_CLT_D44f_z1z4==0, coeff_CLT_D44f_z2z4==0, coeff_CLT_D44f_z3z4==0], c3, c4, c9, c10, d4, d5, d6, d17, d18, d19, d20, d21, d22, d23, solution_dict=True)
c3(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][c3]).full_simplify()
c4(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][c4]).full_simplify()
c9(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][c9]).full_simplify()
c10(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][c10]).full_simplify()
d4(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][d4]).full_simplify()
d5(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][d5]).full_simplify()
d6(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][d6]).full_simplify()
d17(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][d17]).full_simplify()
d18(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][d18]).full_simplify()
d19(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][d19]).full_simplify()
d20(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][d20]).full_simplify()
d21(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][d21]).full_simplify()
d22(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][d22]).full_simplify()
d23(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cB[0][d23]).full_simplify()
c3_sol = c3(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
c4_sol = c4(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
c9_sol = c9(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
c10_sol = c10(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d4_sol = d4(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d5_sol = d5(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d6_sol = d6(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d17_sol = d17(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d18_sol = d18(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d19_sol = d19(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d20_sol = d20(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d21_sol = d21(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d22_sol = d22(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d23_sol = d23(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()


sol_cC = solve([coeff_CLT_D24f_z1==0, coeff_CLT_D24f_z2==0, coeff_CLT_D24f_z3==0, coeff_CLT_D24f_z4==0, coeff_CLT_D24f_z1z1==0, coeff_CLT_D24f_z2z2==0, coeff_CLT_D24f_z3z3==0, coeff_CLT_D24f_z4z4==0, coeff_CLT_D24f_z1z2==0, coeff_CLT_D24f_z1z3==0, coeff_CLT_D24f_z2z3==0, coeff_CLT_D24f_z1z4==0, coeff_CLT_D24f_z2z4==0, coeff_CLT_D24f_z3z4==0], c5, c6, c11, c12, d7, d8, d9, d24, d25, d26, d27, d28, d29, d30, solution_dict=True)
c5(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][c5]).full_simplify()
c6(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][c6]).full_simplify()
c11(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][c11]).full_simplify()
c12(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][c12]).full_simplify()
d7(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][d7]).full_simplify()
d8(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][d8]).full_simplify()
d9(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][d9]).full_simplify()
d24(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][d24]).full_simplify()
d25(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][d25]).full_simplify()
d26(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][d26]).full_simplify()
d27(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][d27]).full_simplify()
d28(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][d28]).full_simplify()
d29(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][d29]).full_simplify()
d30(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8) =  (sol_cC[0][d30]).full_simplify()
c5_sol = c5(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
c6_sol = c6(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
c11_sol = c11(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
c12_sol = c12(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d7_sol = d7(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d8_sol = d8(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d9_sol = d9(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d24_sol = d24(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d25_sol = d25(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d26_sol = d26(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d27_sol = d27(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d28_sol = d28(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d29_sol = d29(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()
d30_sol = d30(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, a7_sol, a8_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol, b7_sol, b8_sol).full_simplify()

sigma2_CLT_D22f(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,c1,c2,c7,c8,d1,d2,d3,d10,d11,d12,d13,d14,d15,d16,c3,c4,c9,c10,d4,d5,d6,d17,d18,d19,d20,d21,d22,d23,c5,c6,c11,c12,d7,d8,d9,d24,d25,d26,d27,d28,d29,d30) = coeff_CLT_D22f(0,0).full_simplify()
sigma2_CLT_D22f_sol = sigma2_CLT_D22f(a1_sol,a2_sol,a3_sol,a4_sol,a5_sol,a6_sol,a7_sol,a8_sol,b1_sol,b2_sol,b3_sol,b4_sol,b5_sol,b6_sol,b7_sol,b8_sol,c1_sol,c2_sol,c7_sol,c8_sol,d1_sol,d2_sol,d3_sol,d10_sol,d11_sol,d12_sol,d13_sol,d14_sol,d15_sol,d16_sol,c3_sol,c4_sol,c9_sol,c10_sol,d4_sol,d5_sol,d6_sol,d17_sol,d18_sol,d19_sol,d20_sol,d21_sol,d22_sol,d23_sol,c5_sol,c6_sol,c11_sol,c12_sol,d7_sol,d8_sol,d9_sol,d24_sol,d25_sol,d26_sol,d27_sol,d28_sol,d29_sol,d30_sol).full_simplify()
sigma2_CLT_D44f(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,c1,c2,c7,c8,d1,d2,d3,d10,d11,d12,d13,d14,d15,d16,c3,c4,c9,c10,d4,d5,d6,d17,d18,d19,d20,d21,d22,d23,c5,c6,c11,c12,d7,d8,d9,d24,d25,d26,d27,d28,d29,d30) = coeff_CLT_D44f(0,0).full_simplify()
sigma2_CLT_D44f_sol = sigma2_CLT_D44f(a1_sol,a2_sol,a3_sol,a4_sol,a5_sol,a6_sol,a7_sol,a8_sol,b1_sol,b2_sol,b3_sol,b4_sol,b5_sol,b6_sol,b7_sol,b8_sol,c1_sol,c2_sol,c7_sol,c8_sol,d1_sol,d2_sol,d3_sol,d10_sol,d11_sol,d12_sol,d13_sol,d14_sol,d15_sol,d16_sol,c3_sol,c4_sol,c9_sol,c10_sol,d4_sol,d5_sol,d6_sol,d17_sol,d18_sol,d19_sol,d20_sol,d21_sol,d22_sol,d23_sol,c5_sol,c6_sol,c11_sol,c12_sol,d7_sol,d8_sol,d9_sol,d24_sol,d25_sol,d26_sol,d27_sol,d28_sol,d29_sol,d30_sol).full_simplify()
sigma2_CLT_D24f(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,c1,c2,c7,c8,d1,d2,d3,d10,d11,d12,d13,d14,d15,d16,c3,c4,c9,c10,d4,d5,d6,d17,d18,d19,d20,d21,d22,d23,c5,c6,c11,c12,d7,d8,d9,d24,d25,d26,d27,d28,d29,d30) = coeff_CLT_D24f(0,0).full_simplify()
sigma2_CLT_D24f_sol = sigma2_CLT_D24f(a1_sol,a2_sol,a3_sol,a4_sol,a5_sol,a6_sol,a7_sol,a8_sol,b1_sol,b2_sol,b3_sol,b4_sol,b5_sol,b6_sol,b7_sol,b8_sol,c1_sol,c2_sol,c7_sol,c8_sol,d1_sol,d2_sol,d3_sol,d10_sol,d11_sol,d12_sol,d13_sol,d14_sol,d15_sol,d16_sol,c3_sol,c4_sol,c9_sol,c10_sol,d4_sol,d5_sol,d6_sol,d17_sol,d18_sol,d19_sol,d20_sol,d21_sol,d22_sol,d23_sol,c5_sol,c6_sol,c11_sol,c12_sol,d7_sol,d8_sol,d9_sol,d24_sol,d25_sol,d26_sol,d27_sol,d28_sol,d29_sol,d30_sol).full_simplify()


sigma2_CLT_D22f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = sigma2_CLT_D22f_sol
sigma2_CLT_D44f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = sigma2_CLT_D44f_sol
sigma2_CLT_D24f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = sigma2_CLT_D24f_sol

print("coeff of D22f as N->infty: ", limit(sigma2_CLT_D22f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))
print("coeff of D44f as N->infty: ", limit(sigma2_CLT_D44f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))
print("coeff of D24f as N->infty: ", limit(sigma2_CLT_D24f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))



print("Section 3.8: Signalling proteins with 2 compartments, B slow")

# variables
N,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26,d27,d28,d29,d30,k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42,v2,v4,v6,u2,u4,u6,z1,z2,z3,z4,D2f,D4f,D22f,D44f,D24f,D6f,D66f,D26f,D46f = var('N,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26,d27,d28,d29,d30,k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42,v2,v4,v6,u2,u4,u6,z1,z2,z3,z4,D2f,D4f,D22f,D44f,D24f,D6f,D26f,D46f,D66f')
k11=k13=k14=k15=k16=k17=k18=k21=k22=k23=k24=k25=k26=k27=k28=n12=n22=n32=n11=n21=n31=n41=n42=1
k12=2
# LLN
G0f = (z1*k13-v2*n22/(n22+n21)*k15+(z2-z1)*k23-v2*n21/(n22+n21)*k25)*D2f+(v2*n22/(n21+n22)*v6*n32/(n31+n32)*k17-v4*n42/(n41+n42)*k18+v2*n21/(n22+n21)*v6*n31/(n31+n32)*k27-v4*n41/(n41+n42)*k28)*D4f+(z1*k14-v6*n32/(n32+n31)*k16+(z2-z1)*k24-v6*n31/(n32+n31)*k26)*D6f
g1(z2)=a1*z2*D2f+a2*z2*D4f+a3*z2*D6f
G1g = k11*(g1(z2+1)-g1(z2))+z1*k12*(g1(z2-1)-g1(z2))+k21*(g1(z2+1)-g1(z2))+(z2-z1)*k22*(g1(z2-1)-g1(z2))
g2(z1)=a4*z1*D2f+a5*z1*D4f+a6*z1*D6f
G2g = z1*n11*(g2(z1-1)-g2(z1))+(z2-z1)*n12*(g2(z1+1)-g2(z1))

coeff_LLN_D2f(z1,z2) = (G0f+G1g+G2g).coefficient(D2f)
a(z1) = ((coeff_LLN_D2f(z1,0) - coeff_LLN_D2f(0,0))/z1).full_simplify()
coeff_LLN_D2f_z1 = a(0)
a(z2) = ((coeff_LLN_D2f(0,z2) - coeff_LLN_D2f(0,0))/z2).full_simplify()
coeff_LLN_D2f_z2 = a(0)
coeff_LLN_D4f(z1,z2) = (G0f+G1g+G2g).coefficient(D4f)
a(z1) = ((coeff_LLN_D4f(z1,0) - coeff_LLN_D4f(0,0))/z1).full_simplify()
coeff_LLN_D4f_z1 = a(0)
a(z2) = ((coeff_LLN_D4f(0,z2) - coeff_LLN_D4f(0,0))/z2).full_simplify()
coeff_LLN_D4f_z2 = a(0)
coeff_LLN_D6f(z1,z2) = (G0f+G1g+G2g).coefficient(D6f)
a(z1) = ((coeff_LLN_D6f(z1,0) - coeff_LLN_D6f(0,0))/z1).full_simplify()
coeff_LLN_D6f_z1 = a(0)
a(z2) = ((coeff_LLN_D6f(0,z2) - coeff_LLN_D6f(0,0))/z2).full_simplify()
coeff_LLN_D6f_z2 = a(0)


sol_ab = solve([coeff_LLN_D2f_z1==0, coeff_LLN_D2f_z2==0, coeff_LLN_D4f_z1==0, coeff_LLN_D4f_z2==0, coeff_LLN_D6f_z1==0, coeff_LLN_D6f_z2==0], a1, a2, a3, a4, a5, a6, solution_dict=True)
a1_sol =  (sol_ab[0][a1]).full_simplify()
a2_sol =  (sol_ab[0][a2]).full_simplify()
a3_sol =  (sol_ab[0][a3]).full_simplify()
a4_sol =  (sol_ab[0][a4]).full_simplify()
a5_sol =  (sol_ab[0][a5]).full_simplify()
a6_sol =  (sol_ab[0][a6]).full_simplify()


mu_LLN_2(a1, a2, a3, a4, a5, a6) = coeff_LLN_D2f(0,0).full_simplify()
mu_LLN_2_sol = mu_LLN_2(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol)
mu_LLN_4(a1, a2, a3, a4, a5, a6) = coeff_LLN_D4f(0,0).full_simplify()
mu_LLN_4_sol = mu_LLN_4(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol)
mu_LLN_6(a1, a2, a3, a4, a5, a6) = coeff_LLN_D6f(0,0).full_simplify()
mu_LLN_6_sol = mu_LLN_6(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol)
mu_LLN_2_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = mu_LLN_2_sol
mu_LLN_4_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = mu_LLN_4_sol
mu_LLN_6_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = mu_LLN_6_sol

print("LLN:")
print("coeff of D2f as N->infty: ", limit(mu_LLN_2_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo).full_simplify())
print("coeff of D4f as N->infty: ", limit(mu_LLN_4_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo).full_simplify())
print("coeff of D6f as N->infty: ", limit(mu_LLN_6_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo).full_simplify())

# CLT
k11=k13=k14=k15=k16=k17=k18=k21=k22=k23=k24=k25=k26=k27=k28=n12=n22=n32=n11=n21=n31=n41=n42=1
k12=2
L0f = 1/2*(z1*k13+v2*n22/(n21+n22)*k15+(z2-z1)*k23+v2*n21/(n21+n22)*k25)*D22f+1/2*(z1*k14+v6*n32/(n31+n32)*k16+(z2-z1)*k24+v6*n31/(n31+n32)*k26)*D66f+1/2*(v2*n22/(n22+n21)*v6*n32/(n31+n32)*k17+v4*n42/(n42+n41)*k18+v2*n21/(n22+n21)*v6*n31/(n31+n32)*k27+v4*n41/(n42+n41)*k28)*D44f+(-u2*n22/(n22+n21)*k15-u2*n21/(n22+n21)*k25)*D2f+(-u6*n32/(n32+n31)*k16-u6*n31/(n32+n31)*k26)*D6f+(u2*n22/(n21+n22)*u6*n32/(n32+n31)*k17-u4*n42/(n41+n42)*k18+u2*n21/(n21+n22)*u6*n31/(n32+n31)*k27-u4*n41/(n41+n42)*k28)*D4f
D2g(z2) = a1*z2*D22f + a2*z2*D24f+ a3*z2*D26f
D4g(z2) = a1*z2*D24f + a2*z2*D44f+ a3*z2*D46f
D6g(z2) = a1*z2*D26f + a2*z2*D46f+ a3*z2*D66f
L1g = (z1*k13-v2*n22/(n21+n22)*k15+(z2-z1)*k23-v2*n21/(n22+n21)*k25)*D2g(z2)+(z1*k14-v6*n32/(n31+n32)*k16+(z2-z1)*k24-v6*n31/(n32+n31)*k26)*D6g(z2)+(v2*n22/(n21+n22)*v6*n32/(n32+n31)*k17-v4*n42/(n41+n42)*k18+v2*n21/(n21+n22)*v6*n31/(n32+n31)*k27-v4*n41/(n41+n42)*k28)*D4g(z2) - (-v2+4/3)*D2g(z2)-(1/2*v2*v6-v4)*D4g(z2)-(-v6+4/3)*D6g(z2)

h1(z2)  = b1*z2*D2f + b2*z2*D4f+b3*z2*D6f + (c1*z2+d1*z2*(z2-1)/2) * D22f + (c2*z2+d2*z2*(z2-1)/2) * D44f+(c3*z2+d3*z2*(z2-1)/2) * D66f+(c4*z2+d4*z2*(z2-1)/2) * D24f+(c5*z2+d5*z2*(z2-1)/2) * D26f+(c6*z2+d6*z2*(z2-1)/2) * D46f

L2h1 = k11*(h1(z2+1)-h1(z2))+z1*k12*(h1(z2-1)-h1(z2))+k21*(h1(z2+1)-h1(z2))+(z2-z1)*k22*(h1(z2-1)-h1(z2))

h2(z1,z2)  = b4*z1*D2f + b5*z1*D4f+b6*z1*D6f + (c7*z1+d7*z1*(z1-1)/2 +d8*(z1+z2)*(z1+z2-1)/2) * D22f + (c8*z1+d9*z1*(z1-1)/2 +d10*(z1+z2)*(z1+z2-1)/2) * D44f +(c9*z1+d11*z1*(z1-1)/2 +d12*(z1+z2)*(z1+z2-1)/2) * D66f +(c10*z1+d13*z1*(z1-1)/2 +d14*(z1+z2)*(z1+z2-1)/2) * D24f +(c11*z1+d15*z1*(z1-1)/2 +d16*(z1+z2)*(z1+z2-1)/2) * D26f+ (c12*z1+d17*z1*(z1-1)/2 +d18*(z1+z2)*(z1+z2-1)/2) * D46f

L2h2 = z1*n11*(h2(z1-1,z2)-h2(z1,z2))+(z2-z1)*n12*(h2(z1+1,z2)-h2(z1,z2))

# Drift in CLT
coeff_CLT_D2f(z1,z2) = (L0f+L1g+L2h1+L2h2).coefficient(D2f)
coeff_CLT_D4f(z1,z2) = (L0f+L1g+L2h1+L2h2).coefficient(D4f)
coeff_CLT_D6f(z1,z2) = (L0f+L1g+L2h1+L2h2).coefficient(D6f)
a(z1) = ((coeff_CLT_D2f(z1,0) - coeff_CLT_D2f(0,0))/z1).full_simplify()
coeff_CLT_D2f_z1 = a(0)
a(z2) = ((coeff_CLT_D2f(0,z2) - coeff_CLT_D2f(0,0))/z2).full_simplify()
coeff_CLT_D2f_z2 = a(0)
a(z1) = ((coeff_CLT_D4f(z1,0) - coeff_CLT_D4f(0,0))/z1).full_simplify()
coeff_CLT_D4f_z1 = a(0)
a(z2) = ((coeff_CLT_D4f(0,z2) - coeff_CLT_D4f(0,0))/z2).full_simplify()
coeff_CLT_D4f_z2 = a(0)
a(z1) = ((coeff_CLT_D6f(z1,0) - coeff_CLT_D6f(0,0))/z1).full_simplify()
coeff_CLT_D6f_z1 = a(0)
a(z2) = ((coeff_CLT_D6f(0,z2) - coeff_CLT_D6f(0,0))/z2).full_simplify()
coeff_CLT_D6f_z2 = a(0)


sol_b = solve([coeff_CLT_D2f_z1==0, coeff_CLT_D2f_z2==0, coeff_CLT_D4f_z1==0, coeff_CLT_D4f_z2==0, coeff_CLT_D6f_z1==0, coeff_CLT_D6f_z2==0], b1, b2, b3, b4, b5, b6, solution_dict=True)
b1(a1,a2,a3,a4,a5,a6) =  (sol_b[0][b1]).full_simplify()
b2(a1,a2,a3,a4,a5,a6) =  (sol_b[0][b2]).full_simplify()
b3(a1,a2,a3,a4,a5,a6) =  (sol_b[0][b3]).full_simplify()
b4(a1,a2,a3,a4,a5,a6) =  (sol_b[0][b4]).full_simplify()
b5(a1,a2,a3,a4,a5,a6) =  (sol_b[0][b5]).full_simplify()
b6(a1,a2,a3,a4,a5,a6) =  (sol_b[0][b6]).full_simplify()
b1_sol = b1(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol).full_simplify()
b2_sol = b2(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol).full_simplify()
b3_sol = b3(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol).full_simplify()
b4_sol = b4(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol).full_simplify()
b5_sol = b5(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol).full_simplify()
b6_sol = b6(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol).full_simplify()



mu_CLT_D2f(a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6) = coeff_CLT_D2f(0,0).full_simplify()
mu_CLT_D2f_sol = mu_CLT_D2f(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol)
mu_CLT_D4f(a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6) = coeff_CLT_D4f(0,0).full_simplify()
mu_CLT_D4f_sol = mu_CLT_D4f(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol)
mu_CLT_D6f(a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6) = coeff_CLT_D6f(0,0).full_simplify()
mu_CLT_D6f_sol = mu_CLT_D6f(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol)
mu_CLT_D2f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = mu_CLT_D2f_sol
mu_CLT_D4f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = mu_CLT_D4f_sol
mu_CLT_D6f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = mu_CLT_D6f_sol

print("CLT:")
print("coeff of D2f as N->infty: ", limit(mu_CLT_D2f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))
print("coeff of D4f as N->infty: ", limit(mu_CLT_D4f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))
print("coeff of D6f as N->infty: ", limit(mu_CLT_D6f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))



# Diffusiont in CLT
coeff_CLT_D22f(z1,z2) = (L0f+L1g+L2h1+L2h2).coefficient(D22f).full_simplify()
coeff_CLT_D44f(z1,z2) = (L0f+L1g+L2h1+L2h2).coefficient(D44f).full_simplify()
coeff_CLT_D66f(z1,z2) = (L0f+L1g+L2h1+L2h2).coefficient(D66f).full_simplify()
coeff_CLT_D24f(z1,z2) = (L0f+L1g+L2h1+L2h2).coefficient(D24f).full_simplify()
coeff_CLT_D26f(z1,z2) = (L0f+L1g+L2h1+L2h2).coefficient(D26f).full_simplify()
coeff_CLT_D46f(z1,z2) = (L0f+L1g+L2h1+L2h2).coefficient(D46f).full_simplify()

a(z1) = ((coeff_CLT_D22f(z1,0) - coeff_CLT_D22f(0,0))/z1).full_simplify()
coeff_CLT_D22f_z1 = a(0).full_simplify()
coeff_CLT_D22f_z1z1 = ((coeff_CLT_D22f(z1,0) - z1*a(0) - coeff_CLT_D22f(0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D22f(0,z2) - coeff_CLT_D22f(0,0))/z2).full_simplify()
coeff_CLT_D22f_z2 = a(0).full_simplify()
coeff_CLT_D22f_z2z2 = ((coeff_CLT_D22f(0,z2) - z2*a(0) - coeff_CLT_D22f(0,0))/z2/z2).full_simplify()
coeff_CLT_D22f_z1z2 = ((coeff_CLT_D22f(z1,z2) - z1*coeff_CLT_D22f_z1 - z2*coeff_CLT_D22f_z2 - z1^2*coeff_CLT_D22f_z1z1 - z2^2*coeff_CLT_D22f_z2z2 - coeff_CLT_D22f(0,0))/z1/z2).full_simplify()


a(z1) = ((coeff_CLT_D44f(z1,0) - coeff_CLT_D44f(0,0))/z1).full_simplify()
coeff_CLT_D44f_z1 = a(0).full_simplify()
coeff_CLT_D44f_z1z1 = ((coeff_CLT_D44f(z1,0) - z1*a(0) - coeff_CLT_D44f(0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D44f(0,z2) - coeff_CLT_D44f(0,0))/z2).full_simplify()
coeff_CLT_D44f_z2 = a(0).full_simplify()
coeff_CLT_D44f_z2z2 = ((coeff_CLT_D44f(0,z2) - z2*a(0) - coeff_CLT_D44f(0,0))/z2/z2).full_simplify()
coeff_CLT_D44f_z1z2 = ((coeff_CLT_D44f(z1,z2) - z1*coeff_CLT_D44f_z1 - z2*coeff_CLT_D44f_z2 - z1^2*coeff_CLT_D44f_z1z1 - z2^2*coeff_CLT_D44f_z2z2 - coeff_CLT_D44f(0,0))/z1/z2).full_simplify()


a(z1) = ((coeff_CLT_D66f(z1,0) - coeff_CLT_D66f(0,0))/z1).full_simplify()
coeff_CLT_D66f_z1 = a(0).full_simplify()
coeff_CLT_D66f_z1z1 = ((coeff_CLT_D66f(z1,0) - z1*a(0) - coeff_CLT_D66f(0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D66f(0,z2) - coeff_CLT_D66f(0,0))/z2).full_simplify()
coeff_CLT_D66f_z2 = a(0).full_simplify()
coeff_CLT_D66f_z2z2 = ((coeff_CLT_D66f(0,z2) - z2*a(0) - coeff_CLT_D66f(0,0))/z2/z2).full_simplify()
coeff_CLT_D66f_z1z2 = ((coeff_CLT_D66f(z1,z2) - z1*coeff_CLT_D66f_z1 - z2*coeff_CLT_D66f_z2 - z1^2*coeff_CLT_D66f_z1z1 - z2^2*coeff_CLT_D66f_z2z2 - coeff_CLT_D66f(0,0))/z1/z2).full_simplify()


a(z1) = ((coeff_CLT_D24f(z1,0) - coeff_CLT_D24f(0,0))/z1).full_simplify()
coeff_CLT_D24f_z1 = a(0).full_simplify()
coeff_CLT_D24f_z1z1 = ((coeff_CLT_D24f(z1,0) - z1*a(0) - coeff_CLT_D24f(0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D24f(0,z2) - coeff_CLT_D24f(0,0))/z2).full_simplify()
coeff_CLT_D24f_z2 = a(0).full_simplify()
coeff_CLT_D24f_z2z2 = ((coeff_CLT_D24f(0,z2) - z2*a(0) - coeff_CLT_D24f(0,0))/z2/z2).full_simplify()
coeff_CLT_D24f_z1z2 = ((coeff_CLT_D24f(z1,z2) - z1*coeff_CLT_D24f_z1 - z2*coeff_CLT_D24f_z2 - z1^2*coeff_CLT_D24f_z1z1 - z2^2*coeff_CLT_D24f_z2z2 - coeff_CLT_D24f(0,0))/z1/z2).full_simplify()


a(z1) = ((coeff_CLT_D26f(z1,0) - coeff_CLT_D26f(0,0))/z1).full_simplify()
coeff_CLT_D26f_z1 = a(0).full_simplify()
coeff_CLT_D26f_z1z1 = ((coeff_CLT_D26f(z1,0) - z1*a(0) - coeff_CLT_D26f(0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D26f(0,z2) - coeff_CLT_D26f(0,0))/z2).full_simplify()
coeff_CLT_D26f_z2 = a(0).full_simplify()
coeff_CLT_D26f_z2z2 = ((coeff_CLT_D26f(0,z2) - z2*a(0) - coeff_CLT_D26f(0,0))/z2/z2).full_simplify()
coeff_CLT_D26f_z1z2 = ((coeff_CLT_D26f(z1,z2) - z1*coeff_CLT_D26f_z1 - z2*coeff_CLT_D26f_z2 - z1^2*coeff_CLT_D26f_z1z1 - z2^2*coeff_CLT_D26f_z2z2 - coeff_CLT_D26f(0,0))/z1/z2).full_simplify()


a(z1) = ((coeff_CLT_D46f(z1,0) - coeff_CLT_D46f(0,0))/z1).full_simplify()
coeff_CLT_D46f_z1 = a(0).full_simplify()
coeff_CLT_D46f_z1z1 = ((coeff_CLT_D46f(z1,0) - z1*a(0) - coeff_CLT_D46f(0,0))/z1/z1).full_simplify()
a(z2) = ((coeff_CLT_D46f(0,z2) - coeff_CLT_D46f(0,0))/z2).full_simplify()
coeff_CLT_D46f_z2 = a(0).full_simplify()
coeff_CLT_D46f_z2z2 = ((coeff_CLT_D46f(0,z2) - z2*a(0) - coeff_CLT_D46f(0,0))/z2/z2).full_simplify()
coeff_CLT_D46f_z1z2 = ((coeff_CLT_D46f(z1,z2) - z1*coeff_CLT_D46f_z1 - z2*coeff_CLT_D46f_z2 - z1^2*coeff_CLT_D46f_z1z1 - z2^2*coeff_CLT_D46f_z2z2 - coeff_CLT_D46f(0,0))/z1/z2).full_simplify()


sol_cA = solve([coeff_CLT_D22f_z1==0, coeff_CLT_D22f_z2==0, coeff_CLT_D22f_z1z1==0, coeff_CLT_D22f_z2z2==0, coeff_CLT_D22f_z1z2==0], c1, c7, d1, d7, d8, solution_dict=True)
c1(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cA[0][c1]).full_simplify()
c7(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cA[0][c7]).full_simplify()
d1(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cA[0][d1]).full_simplify()
d7(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cA[0][d7]).full_simplify()
d8(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cA[0][d8]).full_simplify()
c1_sol = c1(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
c7_sol = c7(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d1_sol = d1(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d7_sol = d7(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d8_sol = d8(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()



sol_cB = solve([coeff_CLT_D44f_z1==0, coeff_CLT_D44f_z2==0, coeff_CLT_D44f_z1z1==0, coeff_CLT_D44f_z2z2==0, coeff_CLT_D44f_z1z2==0], c2, c8, d2, d9, d10, solution_dict=True)
c2(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cB[0][c2]).full_simplify()
c8(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cB[0][c8]).full_simplify()
d2(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cB[0][d2]).full_simplify()
d9(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cB[0][d9]).full_simplify()
d10(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cB[0][d10]).full_simplify()
c2_sol = c2(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
c8_sol = c8(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d2_sol = d2(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d9_sol = d9(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d10_sol = d10(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()


sol_cC = solve([coeff_CLT_D66f_z1==0, coeff_CLT_D66f_z2==0, coeff_CLT_D66f_z1z1==0, coeff_CLT_D66f_z2z2==0, coeff_CLT_D66f_z1z2==0], c3, c9, d3, d11, d12, solution_dict=True)
c3(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cC[0][c3]).full_simplify()
c9(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cC[0][c9]).full_simplify()
d3(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cC[0][d3]).full_simplify()
d11(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cC[0][d11]).full_simplify()
d12(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cC[0][d12]).full_simplify()
c3_sol = c3(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
c9_sol = c9(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d3_sol = d3(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d11_sol = d11(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d12_sol = d12(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()

sol_cD = solve([coeff_CLT_D24f_z1==0, coeff_CLT_D24f_z2==0, coeff_CLT_D24f_z1z1==0, coeff_CLT_D24f_z2z2==0, coeff_CLT_D24f_z1z2==0], c4, c10, d4, d13, d14, solution_dict=True)
c4(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cD[0][c4]).full_simplify()
c10(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cD[0][c10]).full_simplify()
d4(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cD[0][d4]).full_simplify()
d13(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cD[0][d13]).full_simplify()
d14(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cD[0][d14]).full_simplify()
c4_sol = c4(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
c10_sol = c10(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d4_sol = d4(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d13_sol = d13(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d14_sol = d14(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()

sol_cE = solve([coeff_CLT_D26f_z1==0, coeff_CLT_D26f_z2==0, coeff_CLT_D26f_z1z1==0, coeff_CLT_D26f_z2z2==0, coeff_CLT_D26f_z1z2==0], c5, c11, d5, d15, d16, solution_dict=True)
c5(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cE[0][c5]).full_simplify()
c11(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cE[0][c11]).full_simplify()
d5(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cE[0][d5]).full_simplify()
d15(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cE[0][d15]).full_simplify()
d16(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cE[0][d16]).full_simplify()
c5_sol = c5(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
c11_sol = c11(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d5_sol = d5(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d15_sol = d15(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d16_sol = d16(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()

sol_cF = solve([coeff_CLT_D46f_z1==0, coeff_CLT_D46f_z2==0, coeff_CLT_D46f_z1z1==0, coeff_CLT_D46f_z2z2==0, coeff_CLT_D46f_z1z2==0], c6, c12, d6, d17, d18, solution_dict=True)
c6(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cF[0][c6]).full_simplify()
c12(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cF[0][c12]).full_simplify()
d6(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cF[0][d6]).full_simplify()
d17(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cF[0][d17]).full_simplify()
d18(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6) =  (sol_cF[0][d18]).full_simplify()
c6_sol = c6(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
c12_sol = c12(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d6_sol = d6(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d17_sol = d17(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()
d18_sol = d18(a1_sol, a2_sol, a3_sol, a4_sol, a5_sol, a6_sol, b1_sol, b2_sol, b3_sol, b4_sol, b5_sol, b6_sol).full_simplify()



sigma2_CLT_D22f(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,c1,c7,d1,d7,d8,c2,c8,d2,d9,d10,c3,c9,d3,d11,d12,c4,c10,d4,d13,d14,c5,c11,d5,d15,d16,c6,c12,d17,d18) = coeff_CLT_D22f(0,0).full_simplify()
sigma2_CLT_D22f_sol = sigma2_CLT_D22f(a1_sol,a2_sol,a3_sol,a4_sol,a5_sol,a6_sol,b1_sol,b2_sol,b3_sol,b4_sol,b5_sol,b6_sol,c1_sol,c7_sol,d1_sol,d7_sol,d8_sol,c2_sol,c8_sol,d2_sol,d9_sol,d10_sol,c3_sol,c9_sol,d3_sol,d11_sol,d12_sol,c4_sol,c10_sol,d4_sol,d13_sol,d14_sol,c5_sol,c11_sol,d5_sol,d15_sol,d16_sol,c6_sol,c12_sol,d17_sol,d18_sol).full_simplify()
sigma2_CLT_D44f(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,c1,c7,d1,d7,d8,c2,c8,d2,d9,d10,c3,c9,d3,d11,d12,c4,c10,d4,d13,d14,c5,c11,d5,d15,d16,c6,c12,d17,d18) = coeff_CLT_D44f(0,0).full_simplify()
sigma2_CLT_D44f_sol = sigma2_CLT_D44f(a1_sol,a2_sol,a3_sol,a4_sol,a5_sol,a6_sol,b1_sol,b2_sol,b3_sol,b4_sol,b5_sol,b6_sol,c1_sol,c7_sol,d1_sol,d7_sol,d8_sol,c2_sol,c8_sol,d2_sol,d9_sol,d10_sol,c3_sol,c9_sol,d3_sol,d11_sol,d12_sol,c4_sol,c10_sol,d4_sol,d13_sol,d14_sol,c5_sol,c11_sol,d5_sol,d15_sol,d16_sol,c6_sol,c12_sol,d17_sol,d18_sol).full_simplify()
sigma2_CLT_D66f(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,c1,c7,d1,d7,d8,c2,c8,d2,d9,d10,c3,c9,d3,d11,d12,c4,c10,d4,d13,d14,c5,c11,d5,d15,d16,c6,c12,d17,d18) = coeff_CLT_D66f(0,0).full_simplify()
sigma2_CLT_D66f_sol = sigma2_CLT_D66f(a1_sol,a2_sol,a3_sol,a4_sol,a5_sol,a6_sol,b1_sol,b2_sol,b3_sol,b4_sol,b5_sol,b6_sol,c1_sol,c7_sol,d1_sol,d7_sol,d8_sol,c2_sol,c8_sol,d2_sol,d9_sol,d10_sol,c3_sol,c9_sol,d3_sol,d11_sol,d12_sol,c4_sol,c10_sol,d4_sol,d13_sol,d14_sol,c5_sol,c11_sol,d5_sol,d15_sol,d16_sol,c6_sol,c12_sol,d17_sol,d18_sol).full_simplify()
sigma2_CLT_D24f(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,c1,c7,d1,d7,d8,c2,c8,d2,d9,d10,c3,c9,d3,d11,d12,c4,c10,d4,d13,d14,c5,c11,d5,d15,d16,c6,c12,d17,d18) = coeff_CLT_D24f(0,0).full_simplify()
sigma2_CLT_D24f_sol = sigma2_CLT_D24f(a1_sol,a2_sol,a3_sol,a4_sol,a5_sol,a6_sol,b1_sol,b2_sol,b3_sol,b4_sol,b5_sol,b6_sol,c1_sol,c7_sol,d1_sol,d7_sol,d8_sol,c2_sol,c8_sol,d2_sol,d9_sol,d10_sol,c3_sol,c9_sol,d3_sol,d11_sol,d12_sol,c4_sol,c10_sol,d4_sol,d13_sol,d14_sol,c5_sol,c11_sol,d5_sol,d15_sol,d16_sol,c6_sol,c12_sol,d17_sol,d18_sol).full_simplify()
sigma2_CLT_D26f(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,c1,c7,d1,d7,d8,c2,c8,d2,d9,d10,c3,c9,d3,d11,d12,c4,c10,d4,d13,d14,c5,c11,d5,d15,d16,c6,c12,d17,d18) = coeff_CLT_D26f(0,0).full_simplify()
sigma2_CLT_D26f_sol = sigma2_CLT_D26f(a1_sol,a2_sol,a3_sol,a4_sol,a5_sol,a6_sol,b1_sol,b2_sol,b3_sol,b4_sol,b5_sol,b6_sol,c1_sol,c7_sol,d1_sol,d7_sol,d8_sol,c2_sol,c8_sol,d2_sol,d9_sol,d10_sol,c3_sol,c9_sol,d3_sol,d11_sol,d12_sol,c4_sol,c10_sol,d4_sol,d13_sol,d14_sol,c5_sol,c11_sol,d5_sol,d15_sol,d16_sol,c6_sol,c12_sol,d17_sol,d18_sol).full_simplify()
sigma2_CLT_D46f(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,c1,c7,d1,d7,d8,c2,c8,d2,d9,d10,c3,c9,d3,d11,d12,c4,c10,d4,d13,d14,c5,c11,d5,d15,d16,c6,c12,d17,d18) = coeff_CLT_D46f(0,0).full_simplify()
sigma2_CLT_D46f_sol = sigma2_CLT_D46f(a1_sol,a2_sol,a3_sol,a4_sol,a5_sol,a6_sol,b1_sol,b2_sol,b3_sol,b4_sol,b5_sol,b6_sol,c1_sol,c7_sol,d1_sol,d7_sol,d8_sol,c2_sol,c8_sol,d2_sol,d9_sol,d10_sol,c3_sol,c9_sol,d3_sol,d11_sol,d12_sol,c4_sol,c10_sol,d4_sol,d13_sol,d14_sol,c5_sol,c11_sol,d5_sol,d15_sol,d16_sol,c6_sol,c12_sol,d17_sol,d18_sol).full_simplify()




sigma2_CLT_D22f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = sigma2_CLT_D22f_sol
sigma2_CLT_D44f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = sigma2_CLT_D44f_sol
sigma2_CLT_D66f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = sigma2_CLT_D66f_sol
sigma2_CLT_D24f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = sigma2_CLT_D24f_sol
sigma2_CLT_D26f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = sigma2_CLT_D26f_sol
sigma2_CLT_D46f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,n11,n12,n21,n22,n31,n32,n41,n42) = sigma2_CLT_D46f_sol

print("coeff of D22f as N->infty: ", limit(sigma2_CLT_D22f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))
print("coeff of D44f as N->infty: ", limit(sigma2_CLT_D44f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))
print("coeff of D66f as N->infty: ", limit(sigma2_CLT_D66f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))
print("coeff of D24f as N->infty: ", limit(sigma2_CLT_D24f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))
print("coeff of D26f as N->infty: ", limit(sigma2_CLT_D26f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))
print("coeff of D46f as N->infty: ", limit(sigma2_CLT_D46f_sol_fromk(k11,k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,N^2*n11,N^2*n12,N^2*n21,N^2*n22,N^2*n31,N^2*n32,N^2*n41,N^2*n42), N=oo))


