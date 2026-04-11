using Random
using Printf

function monte_carlo_pi(total_points::UInt64)

    rng = Random.default_rng()  # default RNG

    inside_circle = UInt64(0)

    start_time = time()  # start timer

    for i in 1:total_points
        x = rand(rng)
        y = rand(rng)

        if x*x + y*y <= 1.0
            inside_circle += 1
        end
    end

    elapsed = time() - start_time  # elapsed time

    pi_est = 4.0 * inside_circle / total_points

    error_percent = ((pi_est - π) / π) * 100

    @printf("Estimated value of Pi: %.18f\n", pi_est)
    @printf("Error in estimated value of Pi: %.5f %%\n", error_percent)
    @printf("Execution time: %.6f seconds\n", elapsed)
end


# ================================
# Main (command-line handling)
# ================================
if length(ARGS) < 1
    println("Usage: julia monte_pi.jl <number_of_points>")
    exit(1)
end

try
    N = parse(UInt64, ARGS[1])
    monte_carlo_pi(N)
catch
    println("Error: Input must be a positive integer.")
end