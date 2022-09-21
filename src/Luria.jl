module Luria
using Distributions, CSV, DataFrames, StatsBase, Integrals, Optimization, OptimizationOptimJL, 
		Plots
export read_salvador, opt_param, plot_cdf, plot_diff, plot_all

# To use, obtain file salvadorOutput.zip for Zenodo at

# Reset location of salvadorOutput.zip file as needed
# Unzip zipfile, reads file ls_FILE.csv, returns array with data, deletes unzipped files
# FILE is e.g., "6_7_10" in which 6 => 10^6 samples, 7 => 10^-7 mutation (u), and
#		10 => 10^10 final population (N)
# In most cases do not call this directly.
function read_salvador(file)
	zipf = "/Users/steve/text/Submit/Luria/Analysis/LuriaFit/salvadorOutput"
	run(`unzip $zipf`);
	df = CSV.read("salvadorOutput/ld_"*file*".csv", DataFrame);
	rm("salvadorOutput", recursive=true)
	rm("__MACOSX", recursive=true)
	return df[:,2]
end

# Loss function, could use maximum difference instead of summed squared deviations
# or other loss function. Not called directly.
function opt(p, data, upperb)
	emp = ecdf(data)
	prob = IntegralProblem((x,pp) ->
					(cdf(Frechet(p[1],p[2]),x-p[3]) - emp(x))^2, p[3], upperb)
	solve(prob,HCubatureJL())
end

# Some printing to show progress of numerical optimization of parameters
function callback(p, lossval)
	println("p = ", p, "; loss = ", lossval)
	return false
end

# Call this function to find an optimized fit of parameters to simulated
# distribution data. N and u must match the values used for the simulation.
# file is param for read_salvador, e.g., "6_7_10", for which N=10^10 and u = 10^-7
# returns numerically optimized fit for parameters of frechet
function opt_param(file, N, u; tol=1e-3)
	data = read_salvador(file);
	pp = [1.32, 2.7*N*u, log(N*u) - 2.35]
	prob = OptimizationProblem((p,x) -> opt(p, data, N), pp, ones(3))
	solve(prob, NelderMead(), reltol=tol, callback=callback)
end

# plots the empirical CDF of the simulated data vs the optimized Frechet fit
# p_fr is optimized frechet parameter vector for shape, scale, and location
# pass in as for example [1.36, 271.0, 34.8]
# Nu = N*u, lb,ub are lower and upper bound for plotting, e.g., 1e4, 1e6 for Nu=1e4
function plot_cdf(file, p_fr::Vector{Float64}, Nu::Float64, lb::Float64, ub::Float64)
	data = read_salvador(file);
	emp = ecdf(data)
	plt=plot(x->emp(x),lb,ub,xaxis=:log)
	plot!(x->cdf(Frechet(p_fr[1], p_fr[2]), x-p_fr[3]),lb,ub,xaxis=:log)
	e = MathConstants.e
	# using 1+alpha = 1+e/2 does not change visual fit
	# pp = [e/2, e*Nu, Nu*(log(Nu) - 2.35)] 	# original estimate 2.35
	pp = [e/2, e*Nu, Nu*(log(Nu) - (1+e/2))]	# new guess, (1+e/2) approx 2.359
	pp = [e/2, e*Nu, Nu*log(Nu*e^(- (1+e/2)))]	# rewritten in more elegant form
	plot!(x->cdf(Frechet(pp[1], pp[2]), x-pp[3]),lb,ub,xaxis=:log)
	display(plt)
	return plt
end

# plots empirical CDF vs Frechet using single parameter Nu to calculate the theoretical
# values for Frechet parameters given in the article 
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

# plots for a range of Nu values, used to produce figure in article
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

# plots difference between empirical CDF and Frechet approximation, shows mismatch
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
