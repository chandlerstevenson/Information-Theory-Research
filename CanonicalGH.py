import random
import numpy as np

def generator_matrix(K, N, thresh):
    Q = N - K

    G = np.zeros((K, N), dtype=int)

    # Identity portion
    for n in range(N-K, N):
        for k in range(K):
            if k == n - (N-K):
                G[k, n] = 1

    # Parity portion
    for k in range(K):
        j = random.randint(N-K, N-1)
        G[k, j] = 1

    for n in range(N-K):
        for k in range(K):
            # Skip the presets
            if G[k, n] != 1:
                if random.random() < thresh:
                    G[k, n] = 1

    return G

def parity_check_matrix(K, N, thresh):
    Q = N - K

    G = generator_matrix(K, N, thresh)
    H = np.zeros((Q, N), dtype=int)

    for p in range(Q):
        for n in range(Q):
            if p == n:
                H[p, n] = 1

    for p in range(Q):
        for n in range(N-K, N):
            H[p, n] = G[n-(N-K), p]

    return H

# Example usage
K = 3
N = 6
thresh = 0.5

G = generator_matrix(K, N, thresh)
H = parity_check_matrix(K, N, thresh)

GH = np.mod(np.matmul(G, H.T), 2)

print("G:\n", G)
print("H:\n", H)
print("GH:\n", GH)


def encode_messages(messages, Gen):
    """
    Encode an array of n bit messages by multiplying each one by a generator matrix G.

    Args:
        messages: An array of n bit messages.
        G: A generator matrix.

    Returns:
        An array of codewords produced by encoding each message with the generator matrix G.
    """
    codewords = np.dot(Gen, messages) % 2
    return codewords
 
codeword = np.dot(G.T, [1, 1, 1])  
print(np.mod(codeword, 2))
check = np.mod(np.dot(H, codeword), 2)
print(check)
