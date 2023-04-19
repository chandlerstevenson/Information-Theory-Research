
# GRAND

# K. R. Duffy, J. Li, and M. Medard, "Capacity-achieving guessing random 
# additive noise decoding," IEEE Trans. Inf. Theory, vol. 65, no. 7, pp. 
#  4023–4040, 2019.

import numpy as np  
import math 
import scipy as sp  
import random  
import matplotlib as plt
import itertools 
from itertools import product 
import time # Miles ad
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
mat_begin = time.time()
# Here, we create the Generator and Parity check Matrices:  
np.set_printoptions(threshold=np.inf)
def create_canonical_matrices(n, k, density=0.6):
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
mat_end = time.time() 

mat_time = mat_end - mat_begin
print(f"make matrix: {mat_time}")
n_dim = 6
k_dim = 3
print(f"starting {n_dim}x{k_dim} RLC simulation")
density = 0.5

G, H = create_canonical_matrices(n_dim, k_dim, density)



def generate_error_sequences(length_codeword):
    # generate all possible binary sequences of a given size
    sequences = []
    for i in range(2 ** length_codeword):
        seq = [int(bin_num) for bin_num in bin(i)[2:].zfill(length_codeword)]
        sequences.append(seq)
    return sequences  


def calculate_error_probability(error_vector, prob_error, n_dim):
    num_errors = np.sum(error_vector)
    return (1 - prob_error) ** (n_dim - num_errors) * prob_error ** num_errors

def rank_errors(prob_error, n_dim):
    error_array = np.array(list(product([0, 1], repeat=n_dim)))
    error_array = error_array[1:]

    error_list = []
    for error in error_array:
        error_prob = calculate_error_probability(error, prob_error, n_dim)
        error_list.append((error, error_prob))
    return error_list



def choose_error_vector(error_list):
    # calculate total probability of all error vectors
    total_prob = sum(prob for _, prob in error_list)
    # choose a random number between 0 and the total probability
    rand_num = random.uniform(0, total_prob)

    # iterate through error vectors until cumuslative probability exceeds rand_num
    cum_prob = 0
    for error, prob in error_list:
        cum_prob += prob
        if cum_prob >= rand_num:
            return np.array([error])  # return a 2D array instead of 1D
print("generating noise...")


#Define a random binary list of integers to serve as the message, m 
def random_binary_list(length):
    # Generate a list of random binary integers of given length
    binary_list = [random.randint(0, 1) for _ in range(length)]
    return binary_list
    
message = np.array(random_binary_list(k_dim)) 

notcodeword = np.mod(np.dot(G.T,message), 2) 

check = np.mod(np.dot(H,notcodeword), 2) 





def find_error(received_signal, ranked_errors):
    for error, _ in ranked_errors:
        syndrome = np.dot(H, (received_signal - error) % 2) % 2
        if np.array_equal(syndrome, np.zeros(n_dim - k_dim)):
            decoded_bits = (received_signal % 2).tolist()
            return decoded_bits

    return []



print("finding error...") 


def encode_messages(messages, generator_matrix):
    return np.dot(messages, generator_matrix) % 2



def ber(bits, decoded_bits):
    num_errors = np.sum(bits != decoded_bits)
    return num_errors / len(bits)



def introduce_error(codeword, prob_error):
    error = np.random.choice([0, 1], size=codeword.shape, p=[1 - prob_error, prob_error])
    return (codeword + error) % 2
print("error defined")

num_messages = 30000



Eb_N0_db = np.arange(0, 7.5, 1)
# ...




print("decoding...")

def payload_ber(bits, decoded_bits, k_dim): 
    decoded_payload = decoded_bits[:, :k_dim]
    num_errors = np.sum(bits.ravel() != decoded_payload.ravel()) #Error is here 
    return num_errors / (bits.shape[0] * k_dim)
print("payload defined")

def generate_ranked_errors_for_dbs(Eb_N0_db):
    ranked_errors = {}
    for snr_db in Eb_N0_db:
        prob_error = 0.5 * math.erfc(math.sqrt(10 ** (snr_db / 10)))
        ranked_errors[prob_error] = rank_errors(prob_error, n_dim)
    return ranked_errors
start = time.time()
ranked_errors_for_dbs = generate_ranked_errors_for_dbs(Eb_N0_db) 
stop = time.time() 
result = stop - start
print(f"ranked db error time: {result}")
# Define an empty list to store the decoded bits
block_error_values = []
print("block error vals defined")
decoded_bits = [] 
ber_values = []
payload_ber_values = []

# Precompute the ranked_errors for each probability of error
start = time.time()
ranked_errors_for_dbs = {}
for snr_db in Eb_N0_db:
    prob_error = 0.5 * math.erfc(math.sqrt(10 ** (snr_db / 10)))
    ranked_errors_for_dbs[prob_error] = rank_errors(prob_error, n_dim)
end = time.time()

print('first loop took:',end - start)

for snr_db in Eb_N0_db:
    start = time.time()
    # Clear the decoded_bits list at the beginning of each iteration
    decoded_bits_list = []

    bits = np.random.randint(0, 2, num_messages * k_dim) #these are messages
    messages = bits.reshape(num_messages, k_dim)
    # Calculate the probability of error for the current SNR value
    prob_error = 0.5 * math.erfc(math.sqrt(10 ** (snr_db / 10)))

    codewords = encode_messages(messages, G)
    end1 = time.time()
    print('encoding takes:',end1-start)

    flat_signal = codewords.ravel()

    received_signal = introduce_error(flat_signal, prob_error)
    received_signal = received_signal.reshape(num_messages, n_dim)

    decoding_begin = time.time()
    for codeword in received_signal:
        decoded_bits = find_error(codeword, ranked_errors_for_dbs[prob_error])
        decoded_bits_list.append(decoded_bits) 
    decoding_end = time.time() 
    print("decoding time: ", decoding_end-decoding_begin)

    decoded_bits_array = np.array(decoded_bits_list)
    ber_val = ber(flat_signal, decoded_bits_array.ravel())
    ber_values.append(ber_val)

    payload_begin = time.time()
    payload_ber_val = payload_ber(messages, decoded_bits_array, k_dim)
    payload_ber_values.append(payload_ber_val) 
    payload_end = time.time() 

    print("payload time:", payload_end-payload_begin)







# Calculate probability of error for each Eb/N0 value
prob_error_values = [0.5 * math.erfc(math.sqrt(10 ** (snr_db / 10))) for snr_db in Eb_N0_db]

# Plot Eb/N0 vs. BER
plt.plot(Eb_N0_db, ber_values, marker='o', linestyle='-', label='Codeword BER')

# Plot Eb/N0 vs. probability of error
plt.plot(Eb_N0_db, prob_error_values, marker='s', linestyle='--', color='r', label='Probability of Error')

plt.plot(Eb_N0_db, payload_ber_values, marker='o', linestyle=':', color='g', label='Payload BER')

# plt.plot(Eb_N0_db, block_error_values, marker='o', linestyle=':', color='b', label='Block Error')

# Set y-axis scale to logarithmic
plt.yscale('log')

# Add axis labels, title, and legend
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('BER / Probability of Error')
plt.title(f"Eb/N0 vs. BER, BLER, and Probability of Error ({n_dim}x{k_dim}, n = {num_messages})")
plt.legend()

# Show plot
plt.show()

