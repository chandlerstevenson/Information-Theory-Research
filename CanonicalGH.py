import numpy as np 

np.set_printoptions(threshold=np.inf)
def create_canonical_matrices(n, k, density):
    if n <= k:
        raise ValueError("n must be greater than k")
    if density < 0 or density > 1:
        raise ValueError("density must be between 0 and 1")

    # Create the canonical generator matrix G
    P_continuous = np.random.rand(k, n-k)
    P = (P_continuous < density).astype(int)
    G = np.hstack((np.eye(k, dtype=int), P))

    # Create the parity check matrix H
    H = np.hstack((P.T, np.eye(n-k, dtype=int)))

    # Verify GH is equal to the zero matrix
    product = np.matmul(G, H.T) % 2
    assert np.array_equal(product, np.zeros((k, n-k), dtype=int)), "GH is not equal to the zero matrix"

    return G, H

n = 6 #User input
k = 3 #User input
density = 0.6 #User input
G, H = create_canonical_matrices(n, k, density)

print("Generator matrix G:")
print(G)

print("\nParity check matrix H:")
print(H) 
