using Plots
pyplot()
include("a11.jl")
n = 250
γ = 0.24π # LinRange(0,π/2, n)
rng = LinRange(-2,0,n)
lc =  LinRange(0.01,3,n)
θ = π/4#LinRange(0,π/2,nx)
σ1 = -100
KI = zeros(n)
KII = zeros(n)
μ = 0.3
#Kout  = straightcrack(lc[1],γ,θ,σ1)
for i = 1:n
    (Kout,s,w)  = straightcrack(lc[i],γ,θ,σ1,μ)
    KII[i] = real(Kout[1])
    KI[i] = imag(Kout[1])
end
plot(lc,KI./(abs.(σ1*sqrt.(π./lc))))
# plot(real(s),imag(s))
