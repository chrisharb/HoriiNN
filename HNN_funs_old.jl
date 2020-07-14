function L1(r, z, z0, θ, β, U, F, G) #Checked 29/06/20
    L1 = β*U((F(z, z0)+ F(z, conj(z0))), z0, θ) +
        conj(β)*U((z0-conj(z0))*G(z, z0), z0, θ) -
        1/(z+z0) + exp(2im*θ)*(conj(z)+conj(z0))/(z+z0)^2
end

function L2(r, z, z0, θ, β, U, F, G) #Checked 29/06/20
    L2 = -conj(β)*U((F(z, z0) + F(z, conj(z0))), z0, θ) -
        β*U((conj(z0)-z0)*G(z, conj(z0)), z0, θ) -
        1/(conj(z)+conj(z0)) + exp(2im*θ)*(1/(z+z0))
end

function F(z, z0)
    F = (1-(z/z0)*sqrt(complex(z0^2-1))/sqrt(complex(z^2-1)))*(z0/(z^2-z0^2))
end

function G(z, z0) #Checked 29/06/20
    G = (1-(z/z0)*(sqrt(complex(z0^2-1))/sqrt(complex(z^2-1)))) * (2z0^2/(z^2-z0^2)^2)+
         (1-(z*z0)/(sqrt(complex(z0^2-1)))*sqrt(complex(z^2-1)))/(z^2-z0^2)
end

function U(z, z0, θ) #Checked 29/06/20
    U = -F(z, z0) + F(conj(z), z0) +
    exp(2im*θ)*(2*F(z, z0)+(z-conj(z))*F(z, z0))
end

function S(z, θ) #Checked 29/06/20
    S = 0.5im*(conj(z)/sqrt(complex(z^2-1)) -
        (z/sqrt(complex(z^2-1))) +
        exp(2im*θ)*(((conj(z)-z))/((z^2-1)^(3/2))+
        (2*z)/sqrt(complex(z^2-1))-2))
end

function M1(r, s, θ)
    h = (g(s)-g(r)+1im*(f(s)-f(r)))
    M1 = 1/h-exp(2im*θ)*(conj(h)/h^2)
end

function M2(r, s, θ)
    h = (g(s)-g(r)+1im*(f(s)-f(r)))
    M2 = 1/conj(h)+2*exp(2im*θ)
end

function z(c, s, θ)
    z = c .+s*exp(complex(1im*θ))
end

function z0(c, r, θ)
    z = c .+r*exp(complex(1im*θ))
end

function ξ(n)
    ξ = zeros(n)
    for j = 1:n
        ξ[j] = cos(π*(j-1/2)/n)
    end
    return ξ
end

function x(n)
    x = zeros(n-1)
    for i = 1:n-1
        x[i] = cos(π*i/n)
    end
    return x
end
