using Plots
nx = 200
γ = LinRange(0,π/2, nx)
rng = LinRange(-2,1,nx)
lc =  0.1#10 .^rng#1
θ = π/4#LinRange(0,π/2,nx)
σ1 = -100
KI = zeros(nx)
KII = zeros(nx)
μ = -1
#Kout  = straightcrack(lc[1],γ,θ,σ1)
for i = 1:nx
    Kout  = straightcrack(lc,γ[i],θ,σ1,μ)
    KII[i] = real(Kout[1])
    KI[i] = imag(Kout[1])
end

plot!(γ/π,KI./maximum(abs.(KI)))
