prob = ODEProblem(odesys, u0, tspan, ps)
println("done.")

print("Solving....")
sol = solve(prob, Rodas4P(), reltol=1e-8, abstol=1e-8)
println("done")
