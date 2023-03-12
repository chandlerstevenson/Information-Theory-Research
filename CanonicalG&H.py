# # import numpy as np
# # import random
# # import sys

# # def main(args):
# #     if len(args) != 4:
# #         print("usage: lbcgen.py K N threshold")
# #         return

# #     K = int(args[1])
# #     N = int(args[2])
# #     thresh = float(args[3])
# #     Q = N - K

# #     print(f"{K} X {N} generator")
# #     print(f"Code Rate: {K/N}")
# #     print(f"Parity: {Q}")
# #     print(f"Parity Check Threshold: {thresh}")

# #     if Q < 0:
# #         answer = input(f"{Q} parity bits\nYou sure you want to continue? (y/n): ")
# #         if answer != 'y':
# #             print("exiting")
# #             return

# #     # GENERATOR MATRIX
# #     # Identity portion
# #     G = np.zeros((K, N))
# #     for n in range(N - K, N):
# #         for k in range(K):
# #             if k == (n - (N - K)):
# #                 G[k][n] = 1
# #     # Parity Portion 
# #     # Each of the K rows has to have at LEAST one 1, otherwise the corresponding bit isn't checked
# #     # So preload
# #     for k in range(K):
# #         j = random.randint(N - K, N - 1)
# #         G[k][j] = 1
# #     for n in range(N - K):
# #         for k in range(K):
# #             # Skip the presets
# #             if G[k][n] != 1:
# #                 if random.random() < thresh:
# #                     G[k][n] = 1
# #     print("Generator matrix G:")
# #     print(G)

# #     # PARITY CHECK MATRIX
# #     H = np.zeros((N - K, N))
# #     for p in range(N - K):
# #         for n in range(N - K):
# #             if p == n:
# #                 H[p][n] = 1
# #     for p in range(N - K):
# #         for n in range(N - K, N):
# #             H[p][n] = G[n - (N - K)][p]
# #     print("Parity check matrix H:")
# #     print(H)

# #     # Verify GH = 0
# #     GH = np.dot(G, H.T) % 2
# #     if np.array_equal(GH, np.zeros((K, N - K))):
# #         print("GH = 0, code is linear")
# #     else:
# #         print("GH is not zero, code is not linear")

# # if __name__ == '__main__':
# #     main(sys.argv)

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
