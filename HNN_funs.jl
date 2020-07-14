function L1(z, z0, θ, β, U, F, G, Fdash, Gdash)::ComplexF64 #Checked 29/06/20
    F1(z,z0) = F(z,z0)+F(z,conj(z0))
    F1d(z,z0) = Fdash(z,z0)+Fdash(z,conj(z0))
    F2(z,z0) = (z0-conj(z0))*G(z,z0)
    F2d(z,z0) = (z0-conj(z0))*Gdash(z,z0)
    # U = -Fin(z,z0) + Fin(conj(z),z0) + exp(2im*θ)*(2*Fin(z,z0)+(z-conj(z))*FinD)
    # L1 = -F1(z,z0) + F1(conj(z),z0) + exp(2im*θ)*(2*F1(z,z0)+(z-conj(z))*F1D(z,z0))
    L1 = β*U(F1, F1d, z, z0, θ)+
        conj(β)*U(F2, F2d, z, z0, θ)-
        1/(z+z0)+
        exp(2im*θ)*((conj(z)+conj(z0))/(z+z0)^2)
    # η = -β*(F(z,z0)+F(z,conj(z0)))+(conj(z0)-z0)*conj(β)*G(z,z0)
    # L1 = η + conj(η) + 2z0/(z^2-z0^2) +
    #     exp(2im*θ)*((4z*z0)/(z^2-z0^2)+(2*(z^2+z0^2))/(z^2-z0^2)^2+conj(η)-η)
end

function L2(z::ComplexF64, z0::ComplexF64, θ, β::ComplexF64, U, F, G, Fdash, Gdash)::ComplexF64 #Checked 29/06/20
    F1(z,z0) = F(z,z0)+F(z,conj(z0))
    F1d(z,z0) = Fdash(z,z0)+Fdash(z,conj(z0))
    F2(z,z0) = (conj(z0)-z0)*G(z,z0)
    F2d(z,z0) = (conj(z0)-z0)*Gdash(z,conj(z0))
    L2 = -conj(β)*U(F1, F1d, z, z0, θ) -
        β*U(F2, F2d, z, z0, θ)-
        1/(conj(z)+conj(z0))+
        exp(2im*θ)*1/(z+z0)
    # η = conj(β)*(F(z,z0)+F(z,conj(z0)))+(conj(z0)-z0)*β*G(z,conj(z0))
    # L1 = η + conj(η) + 2*conj(z0)/(conj(z)^2-conj(z0)^2) +
    #     exp(2im*θ)*(2z0/(z^2-z0^2)+conj(η)-η)
end


function F(z::ComplexF64, z0::ComplexF64)::ComplexF64 #Checked 14/07/20
    F = (1-(z/z0)*(z0^2-1)^0.5/(z^2-1)^0.5)*(z0/(z^2-z0^2))
end

function G(z::ComplexF64, z0::ComplexF64)::ComplexF64 #Checked 14/07/20
    G = (1-(z/z0)*(z0^2-1)^0.5/(z^2-1)^0.5)*(2z0^2/(z^2-z0^2)^2)+
          (1-((z*z0)/((z0^2-1)^0.5*(z^2-1)^0.5)))/(z^2-z0^2)
    #G = (2*z^2-z*z0+2z0^2)/(2*(z^2-z0^2)^2) # Based on Mathematica derivation
end

function U(Fin, FinD, z::ComplexF64, z0::ComplexF64, θ)::ComplexF64 #Checked 14/07/20
    # Fin is a function of z and z0 passed to this function
    # FinD is the derivative of the passed function WRT z
    U = -Fin(z,z0) + Fin(conj(z),z0) + exp(2im*θ)*(2*Fin(z,z0)+(conj(z)-z)*FinD(z,z0))
end

function S(z::ComplexF64, θ)::ComplexF64 #Checked 29/06/20
    S = 0.5im*((conj(z)/(conj(z)^2-1)^0.5) -
            (z/(z^2-1)^0.5) +
            exp(2im*θ)*((conj(z)-z)/((z^2-1)^(3/2))+
            (2*z)/(z^2-1)^0.5-2))
    #S = 0.5im*(z/(z^2-1)^0.5 + (1-z)/(z^2-1)^1.5 +
        # conj(z)/(conj(z)^2-1)+
        # exp(2im*θ)*(2z/(z^2-1)^0.5-z/(z^2-1)^1.5-2)-1)
end

# Function of the crack geometry
function M1(r, s, θ)
    h = (g(s)-g(r)+1im*(f(s)-f(r)))
    M1 = 1/h-exp(2im*θ)*(conj(h)/h^2)
end

# Function of the crack geometry
function M2(r, s, θ)
    h = (g(s)-g(r)+1im*(f(s)-f(r)))
    M2 = 1/conj(h)+2*exp(2im*θ)
end

# Function to map collocation points to the complex plane
function z(c, s, θ)::ComplexF64
    z = c .+s*exp(complex(1im*θ))
end

# Function to map chebyshev nodes to the complex plane
function z0(c, r, θ)::ComplexF64
    z = c .+r*exp(complex(1im*θ))
end

# Chebyshev nodes
function ξ(n)
    ξ = zeros(n)
    for j = 1:n
        #ξ[j] = cos(π*(2j-1)/(2n))
        ξ[j] = cos(π*(2j-1)/(2n+1)) #Type 2 nodes
    end
    return ξ
end

# Chebyshev collocation points
function η(n)
    η = zeros(n)
    for i = 1:n
        #x[i] = cos(π*i/n)
        η[i] = cos(π*(2i)/(2n+1)) #Type 2 nodes
    end
    return η
end

function W(n)
    W = π*(1 .+ξ(n))./(n+0.5)
end
# Derivative of F wrt z, derived using Mathematica
function Fd(z::ComplexF64, z0::ComplexF64)::ComplexF64
    #Fd = (z0^2-1)^0.5/((z^2-1)^1.5*(z^2-z0^2))-
    #(2z*(z0*(z^2-1)^0.5)-z*(z0^2-1))/((z^2-1)^0.5*(z^2-z0^2)^2)
    Fd = (8*z*z0-16*z^3*z0+8*z^5*z0-z0^2+z0^4-3*z^4*(z0^2-1)+z^2*(z0^4-1))/
        (4*(z^2-1)^2*(z^2-z0^2)^2) #Mathematica solution d/dz F(z,z0)
end

# Derivative of G wrt z, derived using Mathematica
function Gd(z::ComplexF64,z0::ComplexF64)::ComplexF64
    # Gd = (-4*z^3*(1-z^2)^2+z^2*(3+8*z^2-5*z^2-24*z^4)*z0+
    # 4*z*(1-z^2)^2*(-3+z^2)*z0^2+(1-5*z^2+2*(16+5)*z^4)*z0^3+
    # 12*z*(1-z^2)^2*z0^4-(2*(4+1)-(-8+1)*z^2+
    # 5*z^4)*z0^5+(1+z^2)*z0^7)/(2*(1-z^2)^2*(1-z0^2)*(z^2-z0^2)^3) #Mathematica solution based on G in HNN
    Gd = (-4*z^3+3*z^2*z0-12*z*z0^2+z0^3)/(2*(z^2-z0^2)^3) #Mathematica solution derived from F in HNN
    # a1 = (2*z0^2*((z^2*(z0^2-1))/(z0*(z^2-1)^1.5)-(z0^2-1)^0.5/(z0*(z^2-1)^0.5)))/(z^2-z0^2)^2
    # a2 = (8z*z0^2*(1-(z*(z0^2-1)^1/2)/(z0*(z^2-1)^1/2)))/(z^2-z0^2)^3
    # a3 = ((z^2*z0)/((z^2-1)^3/2*(z0^2-1)^1/2)-z0/((z^2-1)^0.5*(z0^2-1)^0.5))/(z^2-z0^2)
    # a4 = (2z*(1-(z*z0)/((z^2-1)*(z0^2-1))))/(z^2-z0^2)^2
    # Gd = a1 - a2 + a3 - a4
end
