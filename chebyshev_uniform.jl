using PyPlot
n = 100
s = zeros(n)
xp = [-1,0,1]
yp = [0,1,0]
w = π/n
for i = 1:n
    s[i] = cos(π*(2i-1)/2n)
end
W = zeros(n,3)
for i = 1:n
    for j = 1:3
        if i != j
            W[i,j] = w/(xp[j]-s[i])
        end
    end
end
y = zeros(n)
for i = 1:n
    for j = 1:3
        y[i] += W[i,j]
    end
end
