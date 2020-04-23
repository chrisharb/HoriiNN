using Plots
n = 250
γ = LinRange(0,π/2, n)
rng = LinRange(-2,1,n)
lc =  0.1#10 .^rng#1
θ = π/4#LinRange(0,π/2,nx)
σ1 = -10
KI = zeros(n)
KII = zeros(n)
μ = 0.3
#Kout  = straightcrack(lc[1],γ,θ,σ1)
for i = 1:n
    (Kout,s)  = straightcrack(lc,γ[i],θ,σ1,μ)
    KII[i] = real(Kout[1])
    KI[i] = imag(Kout[1])
end
plot(γ/π,KII./((abs(σ1))*sqrt(π/lc)))
# plot(real(s),imag(s))
