using PyPlot
f(x) = 1/(1+25x^2)
k = 5
xp = ones(k)
for i = 1:k
    xp[i] = cos(π*(2*i-1)/(2*(k+1)))
end
yp = f.(xp)
x = LinRange(-1,1,m)
w = ones(k,m)
for i = 1:k
    for j = 1:k
        if i != j
            w[i,:] = w[i,:].*(x.-xp[j])./(xp[i]-xp[j])
        end
    end
end
y = zeros(m)
for i = 1:k
    for j = 1:m
        y[j] += w[i,j]*yp[i]
    end
end
plot(x,f.(x))
plot(x,y)
