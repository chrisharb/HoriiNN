using PyPlot
include("HNN_funs.jl")
tot = 1
L = 1#[0.01,0.1,0.5,1]#range(1,10,length=tot)#[0.01,0.1,0.5,1]
θ = range(π/4,length=tot)
γ = π/5#range(π/20,π/2,length=tot)
KI = zeros(tot)
KII = zeros(tot)
χ1 = zeros(ComplexF64,tot)
τc = 0
μ = 0.3
β = (1.0+1.0im*μ)/2
n = 10
λ = -0
# Due to a misunderstanding all indexing and notations are similar to Hills textbook
w = W(n)
#fig, (ax1,ax2) = subplots(ncols=2,nrows=1)
I = Array{Float64,1}(undef,n)
I1 = Array{Float64,1}(undef,n)
for i = 1:n
    I[i] = 2*sin((i*π/(2n+1))*(2n-1))/sin((i*π)/(2n+1))/(2n+1)
    I1[i] = cot((2i-1)/(2n+1)*π/2)*sin((2i-1)/(2n+1)*n*π)
end

h=1
for J = 1:length(L)
    l = L[J]
    for k = 1:tot
        s = l*(η(n).+1)/2 #
        r = l*(ξ(n).+1)/2
        Z = z.(1.0, s, θ[k]) # Complex mapping of collocation points. Performs correctly
        Z0 = z0.(1.0, r, θ[k]) # Complex mapping of integration points. Performs correctly
        E = zeros(ComplexF64,n,n)
        L1₀ = Array{ComplexF64,2}(undef,n,n)
        L2₀ = Array{ComplexF64,2}(undef,n,n)
        S₀ = Array{ComplexF64,1}(undef,n)
        for i = 1:n
            for j = 1:n
                E[i,j] = w[j]/(l*(s[i]+1)/2-l*(r[j]+1)/2)
                L1₀[i,j] = w[j]*L1(Z[i], Z0[j], θ[k], β, U, F, G, Fd, Gd)
                L2₀[i,j] = w[j]*L2(Z[i], Z0[j], θ[k], β, U, F, G, Fd, Gd)
            end
            S₀[i] = S(Z[i], θ[k])
        end

        M = 2*exp(2im*θ[k])
        σx = 0.5*(1+λ)+0.5*(1-λ)*cos(2γ[h])
        σy = 0.5*(1+λ)-0.5*(1-λ)*cos(2γ[h])
        σθ = 0.5*(1+λ)-0.5*(1-λ)*cos(2θ[k])
        τxy = 0.5*(1-λ)*sin(2γ[h])
        τθ = 0.5*(1-λ)*sin(2θ[k])
        ϕ∞₀ = (σy+σx)/4
        ψ∞₀ = (σy+σx)/2*exp(-2im*γ[h])
        ϕ∞R = 0.5*(τxy-μ*σy+τc)*1im.*(1 .-Z ./(Z.^2 .-1).^0.5)
        ψ∞R = 0.5*(τxy-μ*σy+τc)*1im.*(-Z./(Z.^2 .-1).^(3/2) .+2Z./(Z.^2 .-1).^0.5 .-2)
        ϕ∞ = ϕ∞₀ .+ϕ∞R
        ψ∞ = ψ∞₀ .+ψ∞R
        ϕ2∞ = 0.5im*(τxy-μ*σy+τc).*(Z.^2 ./(Z.^2 .-1).^1.5 .-1 ./(Z.^2 .-1).^0.5)
        ψ1∞ = ϕ∞₀+conj(ϕ∞₀)+exp(2im*θ[k])*(ψ∞₀)
        #σf = ϕ∞.+conj(ϕ∞).+exp(2im*θ[k]).*(conj(Z).*ϕ2∞.+ψ∞) # Alternative calculation for σf, produces identical results
        σf = (τxy-μ*σy+τc)*S₀.+ψ1∞
        σ1 = -real(σf)
        σ2 = -imag(σf)
        ax1.plot(s,σ1)
        ax2.plot(s,σ2)
        A1 = (M.*E +L1₀ +L2₀)
        A2 = (-conj(M)*E +real(L1₀) -conj(L2₀))
        B1 = A1\σ1
        B2 = A2\σ2
        χ = real(B1)+imag(B1)+1im*(imag(B2)+real(B2))
        χ1[k] = 2/(2n+1)*sum(conj(χ).*I)
        # ax1.plot(s,real(χ))
        # ax2.plot(s,imag(χ))
        KI[k] = real((2π)^(1.5)*M*χ1[k])/sqrt(l)
        KII[k] = imag((2π)^(1.5)*M*χ1[k])/sqrt(l)
    end
    # ax1.plot(θ/π,KI,string(l/2+0.25))
    # # ax1.set_ylim(-0.2,0.6)
    # # ax1.set_xlim(0,0.5)
    # ax2.plot(θ/π,KII,string(l/2+0.25))
    # ax2.set_ylim(-0.2,0.6)
    # ax1.set_ylim(0,1)
end
