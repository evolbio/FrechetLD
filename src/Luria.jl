module Luria
using Distributions, CSV, DataFrames, StatsBase, Integrals, Optimization, Optim, Plots
export read_salvador, opt_param, plot_cdf

function read_salvador(param="5_6_9")
	dir = "/Users/steve/text/Submit/Luria/Analysis/LuriaFit/salvadorOutput"
	run(`unzip $dir`);
	df = CSV.read("salvadorOutput/ld_"*param*".csv", DataFrame);
	rm("salvadorOutput", recursive=true)
	rm("__MACOSX", recursive=true)
	return df[:,2]
end

function opt(p, data, upperb)
	emp = ecdf(data)
	prob = IntegralProblem((x,pp) ->
					(cdf(Frechet(p[1],p[2]),x-p[3]) - emp(x))^2, p[3], upperb)
	solve(prob,HCubatureJL())
end

function callback(p, lossval)
	println("p = ", p, "; loss = ", lossval)
	return false
end

function opt_param(upperb)
	data = read_salvador();
	pp = [1.3, 2e4, 6e4]	# for 5_5_9
	pp = [1.3, 2e3, 6e3]	# for 5_6_9
	prob = OptimizationProblem((p,x) -> opt(p, data, upperb), pp, ones(3))
	solve(prob, NelderMead(), reltol=1e-3, callback=callback)
end

function plot_cdf(p_fr::Vector{Float64}, lb::Float64, ub::Float64)
	data = read_salvador();
	emp = ecdf(data)
	plt=plot(x->emp(x),lb,ub,xaxis=:log)
	plot!(x->cdf(Frechet(p_fr[1], p_fr[2]), x-p_fr[3]),lb,ub,xaxis=:log)
	display(plt)
	return plt
end

end # module
