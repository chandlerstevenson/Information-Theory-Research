# GRAND

# K. R. Duffy, J. Li, and M. Medard, "Capacity-achieving guessing random 
# additive noise decoding," IEEE Trans. Inf. Theory, vol. 65, no. 7, pp. 
#  4023–4040, 2019.

from base64 import decode
import numpy as np
import itertools 
import numpy as np  
import math 
import scipy as sp  
import random  
import itertools 
from itertools import product 
import matplotlib.pyplot as plt 
import time 

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

n_dim = 128 #user input
k_dim = 116 #user input
density = 0.5 #user input

G, H = create_canonical_matrices(n_dim, k_dim, density)


def increment_error_locations(locations, err_loc_vec):
    max_index = len(locations) - 1

    # Increment to the next error vector
    while True: 
        for ii in range(len(err_loc_vec) -1, -1, -1):
            if err_loc_vec[ii] == locations[-1]:
                continue
            ind_in_mask = locations.index(err_loc_vec[ii])
            err_loc_vec[ii] = locations[ind_in_mask + 1]
            for jj in range(ii + 1, len(err_loc_vec)):
                new_ind = ind_in_mask + jj - ii
                if new_ind <= max_index:
                    err_loc_vec[jj] = locations[new_ind]
                else:
                    break
            break
        else:
            return None

        # Check if there are any repeated indices in the current combination
        if len(set(err_loc_vec)) == len(err_loc_vec):
            return err_loc_vec



# length = 128
# weight = 3
# locations = list(range(1, length + 1))          TESTING AREA
# err_loc_vec = []

    
def generate_sequence(length, location_vector):
    binary_sequence = [0] * (length)

    ones_placed = 0
    for loc in location_vector: 
        if ones_placed < length:
            binary_sequence[loc - 1] = 1
            ones_placed += 1
        else:
            break

    return binary_sequence


import numpy as np

def check_decoding(nois_codeword, H):
    syndrome = np.dot(H, nois_codeword) % 2
    return np.all(syndrome == 0)


def calculate_error_probability(error_vector, prob_error, n_dim):
    num_errors = np.sum(error_vector)
    return (1 - prob_error) ** (n_dim - num_errors) * prob_error ** num_errors


def encode_messages(messages, generator_matrix):
    return np.dot(messages, generator_matrix) % 2

def ber(codewords, decoded_codewords_list):
    num_errors = 0
    num_compared_bits = 0
    bler_ber = 0 
    for codeword, decoded_codeword in zip(codewords, decoded_codewords_list):
        decod_array = np.array(decoded_codeword)
        if decod_array.size != 0:  # Check if the decoded_codeword is not [] (i.e., decoding was successful)
            num_errors += np.sum(codeword != decoded_codeword)
            num_compared_bits += n_dim 
        else: 
            bler_ber += 1  #increment block error val by 1 

    return (num_errors / num_compared_bits), bler_ber




def introduce_error(codeword, prob_error):
    error = np.random.choice([0, 1], size=codeword.shape, p=[1 - prob_error, prob_error])
    return (codeword + error) % 2



def decoder(noisy_codeword, length, max_weight, H):
    # Check if the given noisy codeword already satisfies Hy = 0
    initial_syndrome = np.dot(H, noisy_codeword) % 2
    if np.array_equal(initial_syndrome, np.zeros(n_dim - k_dim)):
        return np.array(noisy_codeword), 1

    locations = list(range(1, length + 1))
    num_tries = 0
    for weight in range(1, max_weight + 1):
        err_loc_vec = list(range(1, weight + 1))
        while True:
            err_loc_vec = increment_error_locations(locations, err_loc_vec)
            if err_loc_vec is None:
                break
            noise_sequence = generate_sequence(length, err_loc_vec)
            corrected_codeword = np.mod(np.add(noisy_codeword, noise_sequence), 2)
            num_tries += 1
            syndrome = np.dot(H, corrected_codeword) % 2
            if np.array_equal(syndrome, np.zeros(n_dim - k_dim)):
                decoded_bits = np.mod(corrected_codeword, 2).tolist()
                return np.array(decoded_bits), num_tries

    return np.array([]), num_tries





def payload_ber(messages, decoded_codewords_list, k_dim):
    num_compared_bits = 0
    block_errors = 0  
    total_bit_errors = 0 
    for message, decoded_bit in zip(messages, decoded_codewords_list):
        if len(decoded_bit) == 0:  
            block_errors += 1  
        else: 
            decoded_payload = decoded_bit[:k_dim] 
            bit_errors = np.sum(message != decoded_payload) 
            total_bit_errors += bit_errors 
            num_compared_bits += k_dim

    return total_bit_errors / num_compared_bits

def channel_capacity(prob_error):   
    p = prob_error  
    q = 1-p
    Entropy = -p * math.log2(p) - q * math.log2(q)
    chan_cap = 1 - Entropy 
    return chan_cap 


num_messages = 1000
import numpy as np
from scipy.stats import norm

Eb_N0_db = np.arange(1, 7.5, 1)
ber_values = []

bler_values = []
print("decoding...")

payload_ber_values = []
block_error_values = [] 
total_trial_tries = []  
channel_cap_values = []
max_weight = 8 #user input 

for snr_db in Eb_N0_db:
    # Clear the decoded_bits list at the beginning of each iteration
    decoded_bits = []
    # Initialize a counter for failed decoding attempts
    failed_decoding_attempts = 0
    bits = np.random.randint(0, 2, num_messages * k_dim)
    messages = bits.reshape(num_messages, k_dim) 
    # Calculate the probability of error for the current SNR value
    prob_error = 0.5 * math.erfc(math.sqrt(10 ** (snr_db / 10)))
    encoding_begin = time.time()
    codewords = encode_messages(messages, G) 
    encoding_end = time.time() 
    encoding_time = encoding_end - encoding_begin 
    print("encoding_time: ",encoding_time)

    flat_begin = time.time() 
    flat_signal = codewords.ravel()
    received_signal = introduce_error(flat_signal, prob_error)
    flat_end = time.time() 
    print("flatten time: ", flat_end-flat_begin)

    received_signal = received_signal.reshape(num_messages, n_dim)  
    # print(received_signal, "original")
    counter_vals = []
    for codeword in received_signal:
        decoded_codeword, num_tries = decoder(codeword, n_dim, max_weight, H)
        counter_vals.append(num_tries)
        if len(decoded_codeword) != 0 :
            decoded_bits.append(decoded_codeword)  # Append the decoded codeword as a list (or numpy array) 
            # counter_vals.append(counter_num_result)
        else:
            failed_decoding_attempts += 1

    avg_single_trial_tries = np.sum(counter_vals)/num_messages
    print(avg_single_trial_tries)
    total_trial_tries.append(avg_single_trial_tries) 
    # Calculate the block error rate for the current SNR value
    bler = failed_decoding_attempts / num_messages
    bler_values.append(bler)
    
    ber_val,_ = ber(codewords, decoded_bits)  # Pass decoded_bits directly without converting to numpy array
    ber_values.append(ber_val)
    
    payload_ber_val = payload_ber(messages, decoded_bits, k_dim) 
    payload_ber_values.append(payload_ber_val)  

    chan_cap_val = channel_capacity(prob_error) 
    channel_cap_values.append(chan_cap_val)

print(channel_cap_values)  
print(ber_values, "ber values")
print(payload_ber_values, "payload")
print(block_error_values, "block_error")
print(total_trial_tries, "total_tries")
print(channel_cap_values, "channel_cap_values")
print("code rate", k_dim/n_dim)
avg_num_trials_through = np.sum(total_trial_tries)/7 

first_trial = [0.0700234375, 0.0505625, 0.0305703125, 0.0146484375, 0.0045078125, 0.0013671875, 4.6875e-05] 
second_trial = [0.070078125, 0.05025, 0.0305703125, 0.01346875, 0.0045859375, 0.001046875, 9.375e-05]


# print(len(decoded_bits))

# print(ber_values) 
# print(bits)
# Calculate probability of error for each Eb/N0 value
prob_error_values = [0.5 * math.erfc(math.sqrt(10 ** (snr_db / 10))) for snr_db in Eb_N0_db]
# Turn on interactive mode
plt.ion()
plt.figure
plt.plot(Eb_N0_db, ber_values, marker='o', linestyle='-', label='BER')

plt.plot(Eb_N0_db, payload_ber_values, marker='o', linestyle=':', color='g', label='Payload BER')
# Plot Eb/N0 vs. probability of error
plt.plot(Eb_N0_db, prob_error_values, marker='s', linestyle='--', color='r', label='Probability of Error')

# Set y-axis scale to logarithmic
plt.yscale('log')

# Add axis labels, title, and legend
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('BER / Probability of Error')
plt.title(f"(N,k) = ({n_dim},{k_dim}) M = {num_messages}, Z_avg = {avg_num_trials_through: .2f}, W_max = {max_weight}")
plt.legend()
# Show plot
plt.show() 
plt.pause(0.001)  

plt.figure()
plt.plot(Eb_N0_db, total_trial_tries, marker='s', linestyle='--', color='m', label='Avg. Number of Queries') 

plt.xlabel('Eb/N0 (dB)')
plt.ylabel('Avg. Number of Queries per Codeword')
plt.title(f"(N,k) = ({n_dim},{k_dim}) M = {num_messages}, Z_avg = {avg_num_trials_through: .2f}, W_max = {max_weight}")
plt.legend() 
plt.show()
plt.pause(0.001)   


plt.figure()
plt.plot(Eb_N0_db, channel_cap_values, marker='s', linestyle='--', color='b', label='Channel Capacity') 

plt.xlabel('Eb/N0 (dB)')
plt.ylabel('Channel Capacity')
plt.title(f"(N,k) = ({n_dim},{k_dim}) M = {num_messages}, Z_avg = {avg_num_trials_through: .2f}, W_max = {max_weight}")
plt.legend() 
plt.show()
plt.pause(0.001)   

# Turn off interactive mode
plt.ioff()

# Keep the figures open until the user closes them
plt.show(block=True)
