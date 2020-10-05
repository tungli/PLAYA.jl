print("GSA....")
f1 = function (p)
  prob1 = remake(prob;p=p)
  sol = solve(prob1,Rodas5();
              saveat=range(tspan[2]/10., stop=tspan[2], length=100),
              reltol=1e-8,
              abstol=1e-8,
             )
  [mean(sol[6,:]), maximum(sol[22,:])]
end

bounds = [ [0.5*i, 1.5*i] for i in ps ]
m = gsa(f1, Morris(total_num_trajectory=1000, num_trajectory=150), bounds)
println("done")

# output most significant reactions
inds = sortperm(m.means[1,:])
for i in reverse(inds)[1:5]
    r =  rs.eqs[i]
    f = x-> join(map(string, x), " + ")
    println(f(r.substrates), " => ", f(r.products), ":  μ = $(m.means[i]),  σ = $(m.variances[i])")
end
