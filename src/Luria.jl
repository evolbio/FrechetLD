module Luria
using Distributions, CSV, DataFrames, StatsBase, Integrals, Optimization, OptimizationOptimJL, 
		Plots
export read_salvador, opt_param, plot_cdf, plot_diff, plot_all

# file is e.g., "6_7_10"
function read_salvador(file)
	zipf = "/Users/steve/text/Submit/Luria/Analysis/LuriaFit/salvadorOutput"
	run(`unzip $zipf`);
	df = CSV.read("salvadorOutput/ld_"*file*".csv", DataFrame);
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

# file is param for read_salvador, e.g., "6_7_10"
function opt_param(file, N, u; tol=1e-3)
	data = read_salvador(file);
	pp = [1.32, 2.7*N*u, log(N*u) - 2.35]
	prob = OptimizationProblem((p,x) -> opt(p, data, N), pp, ones(3))
	solve(prob, NelderMead(), reltol=tol, callback=callback)
end


function plot_cdf(file, p_fr::Vector{Float64}, Nu::Float64, lb::Float64, ub::Float64)
	data = read_salvador(file);
	emp = ecdf(data)
	plt=plot(x->emp(x),lb,ub,xaxis=:log)
	plot!(x->cdf(Frechet(p_fr[1], p_fr[2]), x-p_fr[3]),lb,ub,xaxis=:log)
	e = MathConstants.e
	pp = [e/2, e*Nu, Nu*(log(Nu) - 2.35)]
	plot!(x->cdf(Frechet(pp[1], pp[2]), x-pp[3]),lb,ub,xaxis=:log)
	display(plt)
	return plt
end

function plot_cdf(file, Nu::Float64; lb::Float64=Nu, ub::Float64=Nu*100.,
			disp=true, xlabel=:none)
	data = read_salvador(file);
	emp = ecdf(data)
	nu = Int(Nu)
	plt=plot(x->emp(x),lb,ub,xaxis=:log, title="Nu = $nu", xlabel=xlabel,
				xguidefontsize=16,titlefontsize=14)
	e = MathConstants.e
	pp = [e/2, e*Nu, Nu*(log(Nu) - 2.35)]
	plot!(x->cdf(Frechet(pp[1], pp[2]), x-pp[3]),lb,ub,xaxis=:log,legend=:none)
	if disp display(plt) end
	return plt
end

function plot_all()
	files = ["7_10_10", "7_9_10", "6_8_10", "6_7_10", "6_6_10"]
	Nu_list = [1e0, 1e1, 1e2, 1e3, 1e4]
	n = length(files)
	plt = Array{Any}(nothing, n)
	for (f,Nu,i) in zip(files,Nu_list,1:n)
		plt[i] = plot_cdf(f,Nu;disp=false,
					xlabel=((i==3) ? "Number of mutants" : ""))
	end
	final = plot(plt...,layout=(1,5), size=(1300,300), bottommargins=50Plots.px)
	display(final)
	return final
end

function plot_diff(file, Nu::Float64; lb::Float64=Nu, ub::Float64=Nu*10.)
	data = read_salvador(file);
	emp = ecdf(data)
	e = MathConstants.e
	pp = [e/2, e*Nu, Nu*(log(Nu) - 2.35)]
	println("\nloc = ", pp[3], "; scale = ", sqrt(pp[2]))
	plt=plot(x->emp(x)-cdf(Frechet(pp[1], pp[2]), x-pp[3]),lb,ub)
	display(plt)
	return plt
end

end # module
