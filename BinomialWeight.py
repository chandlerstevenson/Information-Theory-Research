import math

def binomial(N, w, epsilon):
    binomial_result = math.comb(N, w) * (epsilon ** w) * ((1 - epsilon) ** (N - w))
    return binomial_result

delta = 0.5 # User input: "What % of the distribution should be covered"  
            # For effective decoding, delta >= .5
N_dim = 128 # User Input: Length of the codeword
w_init = 0 # Initializer: should stay zero
eps = 0.01 # User Input: Probability of error values < .1 are reasonable
binomial_sum = binomial(N_dim, w_init, eps)

while binomial_sum < (1-delta) and w_init < N_dim:
    w_init += 1
    binomial_sum = binomial(N_dim, w_init, eps) + binomial_sum
    # print(binomial_sum)

print(f"N = {N_dim}, epsilon = {eps}, optimal w = {w_init}, 1-delta = {1-delta}")
