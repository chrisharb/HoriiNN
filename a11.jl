using Plots

function L1(r, z, z0, θ, β, c)
    L1 = β*U((F(z, z0, c)+ F(z, conj(z0), c)), z0, c, θ) +
        conj(β)*U( (z0-z)*G(z, z0, c), z0, c, θ) -
        1/(z+z0) + exp(2im*θ)*(conj(z)+conj(z0))/(z+z0)^2
end

function L2(r, z, z0, θ, β, c)
    L2 = -conj(β)*U((F( z, z0, c) + F(z, conj(z0), c)), z0, c, θ) -
        β*U( ((conj(z0)-z0)*G(z,conj(z0), c)), z0, c, θ) -
        1/(conj(z)+z0) + exp(2im*θ)*(1/(z+z0))
end

function F(z, z0, c)
    F = (1-(z/z0)*(sqrt(complex(z0^2-c^2))/sqrt(complex(z^2-c^2))))*(z0/(z^2-z0^2))
end

function G(z, z0, c)
    G = (1-(z/z0)*(sqrt(complex(z0^2-c^2))/sqrt(complex(z^2-c^2)))) * (2z0^2/(z^2-z0^2)^2)+
         (1-((z*z0)/sqrt(complex(z0^2-c^2))))/(z^2-z0^2)
end

function U(z, z0, c, θ)
    U = -F(z, z0, c) + F(conj(z), z0, c) +
    exp(2im*θ)*(2*F(z, z0, c)+(z-conj(z)*F(z, z0, c)))
end

function S(z, z0, θ, c)
    S = 0.5im*(conj(z)/sqrt(complex(z^2-c^2)) -
        (z/sqrt(complex(z^2-c^2))) +
        exp(2im*θ)*(((conj(z)-z)*c^2)/((z^2-z0^2)^(3/2))+
        (2*z)/sqrt(complex(z^2-c^2))-2))
end

function M1(r, s, θ)
    h = (g(s)-g(r)+1im*(f(s)-f(r)))
    M1 = 1/h-exp(2im*θ)*(conj(h)/h^2)
end

function M2(r, s, θ)
    h = (g(s)-g(r)+1im*(f(s)-f(r)))
    M2 = 1/conj(h)+2*exp(2im*θ)
end

pyplot()

function CurvedCrack(l)
μ = 0.5 #friction along shear crack
#l =  #wing crack length
c = 0.1 #shear crack length
r = 10 #length along wing crack = l for straight crack
n = 1000 #Steps
β = (1+μ*1im)/2 #complex function of friction

g(r) = (r/10)^(2/3)
Δg(r) = (2/3)*(r/10)^-(1/3)
f(r) = (6r)^(1/3)
Δf(r) = 2r^-(2/3)
s = LinRange(r/n,r,n)
z = c .+ g.(s) .+ 1im*f.(s)
z0 = c + g(r) + 1im*f(r)
θ = atan.(Δf.(s)./Δg.(s))
#plot(real(z),imag(z))
#plot!(imag(z0),real(z0))

Z1 = π/2*((8/l*exp.(θ*1im)) .+
    L1.(r, z, z0, θ, β, c) .+
    L2.(r, z, z0, θ, β, c) .+
    L1.(0, z, z0, θ, β, c) .+
    L2.(0, z, z0, θ, β, c))

Z2 = π/2*(4/l*exp.(θ*1im))-
    L1.(r, z, z0, θ, β, c) .+
    L2.(r, z, z0, θ, β, c)

Ss = S.(z, z0, θ, c)
η = zeros(n,2)
for i = 1:n
    A = [real(Z1[i]),imag(Z2[i])] ; [imag(Z1[i]),-real(Z2[i])]
    x = [(-3-μ*5)*real(Ss[i])+5, (-3-μ*5)*imag(Ss[i])-10]
    η[i , :] .= A\x
end
α = 1 ./(r.*sqrt.(complex(l .-s))).*(η[:,1].*(2*s./l .-1).+1im*η[:,2].*s/l)
plot(real(α[1:end-100]),-imag(α[1:end-100]))
end

function straightcrack(lc,γ,θ,σ1,μ)
    #μ = 0 #friction along shear crack
    r = 1
    c = r/lc #shear crack length

    #r = 0.2 #length along wing crack = l for straight crack
    l = r;
    n = 1000 #Steps
    β = (1+μ*1im)/2 #complex function of friction
    s = Array{Complex{Float64}}(undef,n)
    X = Array{Complex{Float64}}(undef,n)
    for N = 1:n
        s[N] = cos(π*(2N-1)/2n)
        X[N] = cos(π*N/n)
    end
    w = π/n
    #θ = π/6
    z = c.+s.*exp(1im*θ)

    dz = Array{Complex{Float64}}(undef,n)
    dz[2:end] .= z[2:end] .-z[1:end-1]
    dz[1] = z[1]
    z0 = c+r*exp(1im*θ)

    Z1 = π/2*((8/l*exp.(θ*1im)) .-
        L1.(r, z, z0, θ, β, c) .-
        L2.(r, z, z0, θ, β, c) .+
        L1.(0, z, z0, θ, β, c) .+
        L2.(0, z, z0, θ, β, c))

    Z2 = π/2*(4/l*exp.(θ*1im)).+
        L1.(r, z, z0, θ, β, c) .+
        L2.(r, z, z0, θ, β, c)

    Ss = S.(z, z0, θ, c)
    η = zeros(n,2)
    #γ = π/8
    #σ1 = -1e6
    σ2 = 0
    σx∞ = 0.5*(σ1+σ2)+0.5*(σ1-σ2)*cos(2γ)
    σy∞ = 0.5*(σ1+σ2)+0.5*(σ1-σ2)*cos(2γ)
    σθ∞ = 0
    τr∞ = 0
    τ∞ = 0.5*(σ1-σ2)*sin(2γ)
    τc = 0
    σa1 = τ∞-μ*σy∞+τc

    for i = 1:n
        A = [real(Z1[i]) imag(Z2[i]) ; imag(Z1[i]) -real(Z2[i])]
        x = [(τ∞ - μ*σy∞ + τc)*real(Ss[i]) + σθ∞, (τ∞ - μ*σy∞ + τc)*imag(Ss[i]) + τr∞]
        η[i , :] .= A\x
    end

    α = 1 ./sqrt.(complex(s.*(l .-s))).*((2*s./l .-1).+1im*1 .*s./l)

    # # plot(real(α),imag(α))
    #
    K = 0
    for i = 1:n
        K = K+π/n*((2π)^3/2*exp(1im*θ)*(η[i,1]-1im*η[i,2]))*(X[i]-s[i])
    end

    KI = real(K)
    KII = imag(K)
    return K
end
