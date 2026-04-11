using Base.Threads
using Random
using Printf

function monte_carlo_pi_parallel(N::UInt64)

    inside = Threads.Atomic{UInt64}(0)

    start = time()

    Threads.@threads for i in 1:N
        x = rand()
        y = rand()

        if x*x + y*y <= 1.0
            atomic_add!(inside, UInt64(1))
        end
    end

    elapsed = time() - start

    pi_est = 4.0 * inside[] / N
    err = ((pi_est - π) / π) * 100

    @printf("Estimated value of Pi: %.18f\n", pi_est)
    @printf("Error in estimated value of Pi: %.5f %%\n", err)
    @printf("Execution time: %.6f seconds\n", elapsed)
end


# ================================
# Main (command-line handling)
# ================================
if length(ARGS) < 1
    println("Usage: julia monte_pi.jl <number_of_points>")
    exit(1)
end


N = parse(UInt64, ARGS[1])
# print(N)
monte_carlo_pi_parallel(N)
