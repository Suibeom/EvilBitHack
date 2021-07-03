### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 60521604-db74-11eb-184f-c31686bbe37b
using Plots

# ╔═╡ 5da37e73-963d-4a3d-86b0-68c0369f079d
#This notebook is intended to be an interactive companion to Demean's video on the fast inverse square root hack https://www.youtube.com/watch?v=p8u_k2LIZyo

# ╔═╡ adc69f34-e2fb-4c61-a736-2775551e5b75
#We're going to set up a particular trick first; log(x + 1) ~ x + μ
#The title display the average error between the two over [0,1]

# ╔═╡ 4b279d65-7cbf-4ea0-95f6-5f991a8902a5
@bind μ HTML("<input type='range' value='0.0' max='0.1' step='0.0001'")

# ╔═╡ daea3956-d6ac-4bc8-9607-bdbfbc666823
begin
	x = 0:0.01:1
	y = [(x .+ μ) log.(2,x .+ 1)]
	e = sum(abs.((x .+ μ)-log.(2,x .+ 1)))/length(x)
	
	plot(x,y, leg=false, title=string(e))
end

# ╔═╡ 0da371ad-fd8a-4c82-adee-33219cc6c290
#Once we've picked μ, we can put together the 'magic number' we need to line up the bitwise representations of y and 1/√(y)

# ╔═╡ 7ca61f82-43c5-4c23-8842-bbfdee65683e
magic_number = round(Int32, (3/2 * 2^23 * ( 127 - μ)))

# ╔═╡ b340938f-d2a1-4627-9f20-0e3c0a042f14
string(magic_number, base=16)

# ╔═╡ d6824e83-6fa7-47a0-8a7b-95b4a5292b0b
#First, lets see how close this bit-shifting trickery gets us without the Newton's Method step. This implementation excludes it, and we'll graph it against Julia's.

# ╔═╡ b2cac6c8-a8d6-4e74-ba87-24805ee496f2
function fast_inv_sq(x::Float32)
	float_bits = "0"*bitstring(x)[1:end-1]
	fuzzy_log = parse(Int32,float_bits,base=2)
	log_invs_sq = magic_number - fuzzy_log
	float_bits = bitstring(log_invs_sq)
	print(float_bits)
	#	Recasting the bits to a float is trivial in C, but Julia
	#	absolutely hates this!
	exp = parse(Int32, float_bits[2:9],base=2)
	mant = parse(Int32, float_bits[10:end], base=2)
	return 2.0^(exp-127)*(1 + mant/(2.0^23))
	
end

# ╔═╡ 0d2e4834-0322-47e4-b464-4a0bd72f28a0
#Not bad! Depending on the circumstances, I'd wonder if this would be good enough on its own. Gamers have certainly tolerated worst imprecision-based visual artifacts...

# ╔═╡ dc696731-29cf-4687-95ce-408a9b296fd1
begin
	f_x = Float32.(1.0:0.1:100.0)
	f_y = [fast_inv_sq.(f_x) f_x.^(-1/2)]
	e2 = sum(abs.(fast_inv_sq.(f_x)-f_x.^(-1/2)))/length(f_x)
	
	plot(f_x,f_y,title=string(e2))
end

# ╔═╡ d2ad81a6-b159-428c-bf95-e7637f5d63aa
#Now we add the Newton's Method step. This makes a hell of a difference.

# ╔═╡ bf2a4a66-4fb4-490c-bf0b-a11ce1066a58
function fast_inv_sq_v2(x::Float32)
	x2 = 0.5 * x
	float_bits = "0"*bitstring(x)[1:end-1]
	fuzzy_log = parse(Int32,float_bits,base=2)
	log_invs_sq = magic_number - fuzzy_log
	float_bits = bitstring(log_invs_sq)
	print(float_bits)
	exp = parse(Int32, float_bits[2:9],base=2)
	mant = parse(Int32, float_bits[10:end], base=2)
	y = 2.0^(exp-127)*(1 + mant/(2.0^23))
	#	Now for the single step of Newton's method
	y2 = y*(1.5 - x2 * y * y)
	return y2
	
end

# ╔═╡ 26cfea0b-afe9-49f7-bc50-138b88fcfaaf
begin
	ff_x = Float32.(1.0:0.1:100.0)
	ff_y = [fast_inv_sq_v2.(ff_x) ff_x.^(-1/2)]
	e3 = sum(abs.(fast_inv_sq_v2.(ff_x)-ff_x.^(-1/2)))/length(ff_x)
	plot(ff_x,ff_y, title=string(e3))
end

# ╔═╡ Cell order:
# ╠═60521604-db74-11eb-184f-c31686bbe37b
# ╠═5da37e73-963d-4a3d-86b0-68c0369f079d
# ╠═adc69f34-e2fb-4c61-a736-2775551e5b75
# ╠═4b279d65-7cbf-4ea0-95f6-5f991a8902a5
# ╠═daea3956-d6ac-4bc8-9607-bdbfbc666823
# ╠═0da371ad-fd8a-4c82-adee-33219cc6c290
# ╠═7ca61f82-43c5-4c23-8842-bbfdee65683e
# ╠═b340938f-d2a1-4627-9f20-0e3c0a042f14
# ╠═d6824e83-6fa7-47a0-8a7b-95b4a5292b0b
# ╠═b2cac6c8-a8d6-4e74-ba87-24805ee496f2
# ╠═0d2e4834-0322-47e4-b464-4a0bd72f28a0
# ╠═dc696731-29cf-4687-95ce-408a9b296fd1
# ╠═d2ad81a6-b159-428c-bf95-e7637f5d63aa
# ╠═bf2a4a66-4fb4-490c-bf0b-a11ce1066a58
# ╠═26cfea0b-afe9-49f7-bc50-138b88fcfaaf
