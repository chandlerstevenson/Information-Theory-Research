import time
import numpy as np


#Helps display entire matrix 
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

n_dim = 32 #User Input
k_dim = 26 #User Input
density = .5 #User Input 

G, H = create_canonical_matrices(n_dim, k_dim, density)

#NOTE: BERSimulator.py operates on this loop for each probability of error. These results  
#      are only for a single trial. Naturally, with higher probabilites of error will come 
#      a need to look deeper into the noise sequences => P increases 

binary_list = [int(x) for x in bin(n_dim)[2:].zfill(n_dim)]

P = 3 #Use Input (rough estimation of partition of noise sequences:  
       #For example, P = 3 implies the result is on average, found in the upper 3rd 
       #of noise sequences  

num_messages = 1000 #User Input (number of messages in BERSimulator.py)
# Time the matrix-vector multiplication
start_time = time.time()
for i in range(num_messages* int((2**n_dim)/P)): 
    y = np.dot(H, binary_list) % 2
end_time = time.time()

# Calculate the total time
total_time = end_time - start_time
print("Total time:", total_time, "seconds")
