import numpy as np  
import math 
import scipy as sp  
import random 

#------------------------------------------------------------------------------#
'''
Preliminary work:  
* Establish ranked order of noise sequences 

''' 
def rank_errors(channel_prob, error_size):
    # generate all possible error sequences of a given size 
    error_sequences = generate_error_sequences(error_size)

    # calculate the probability of each error sequence
    error_probs = {}
    for error_seq in error_sequences:
        error_prob = channel_prob ** sum(error_seq) * (1 - channel_prob) ** (error_size - sum(error_seq))
        error_probs[''.join(map(str, error_seq))] = error_prob

    # sort error probabilities in descending order
    sorted_errors = sorted(error_probs.items(), key=lambda x: x[1], reverse=True)

    # return a list of tuples, where each tuple contains an error sequence and its probability
    return [(error, prob) for error, prob in sorted_errors]


def generate_error_sequences(length_codeword):
    # generate all possible binary sequences of a given size
    sequences = []
    for i in range(2 ** length_codeword):
        seq = [int(bin_num) for bin_num in bin(i)[2:].zfill(length_codeword)]
        sequences.append(seq)
    return sequences  


#-------------------------------------------------------
# INPUT: rank_errors(a,b)
#   a: channel error probability  
#   b: length of codewords, N 
full_error_list = rank_errors(.48, 6)
for error1 in full_error_list:
    error1
genstr_lst = [t[0] for t in full_error_list] 
genarr_lst = [list(map(int, s)) for s in genstr_lst]
print(full_error_list) 
#
#--------------------------------------------------------

#------------------------------------------------------------
# CHOOSE AN ERROR VECTOR BASED ON PROBABILITY 
#
# Define an empty list to store the new tuples
new_error_seq = []
# Iterate through the original list of tuples
for string, float_value in full_error_list:
    # Split the string into a list of integers
    string_list = [int(x) for x in string]
    # Create a new tuple with the list and the float value
    new_tuple = (string_list, float_value)
    # Append the new tuple to the new list
    new_error_seq.append(new_tuple)

# print(new_error_seq)

def choose_error_vector(error_list):
    # calculate total probability of all error vectors
    total_prob = sum(prob for _, prob in error_list) 
    # print(total_prob)
    prob_alone = [t[1] for t in error_list]
    # choose a random number between 0 and the total probability
    rand_num = random.uniform(0, total_prob/1.5) 
    # print(prob_alone)
    # rand_num = np.random.choice(prob_alone, p=[prob/total_prob for _, prob in error_list])
    # print(rand_num)
    # iterate through error vectors until cumuslative probability exceeds rand_num
    cum_prob = 0
    for error, prob in error_list:
        cum_prob += prob
        if cum_prob >= rand_num:
            return error

# Call the choose_error_vector function to select an error vector
chosen_error = np.mod(choose_error_vector(new_error_seq), 2)
print(f"error: {chosen_error}")
# Print the selected error vector
# print(f"chosen error: {chosen_error}") 
#--------------------------------------------------------------------------------------
#CHOOSE A GENERATOR AND PARITY CHECK MATRIX (tbd for function):  

G = np.array(\
[
[1, 0, 0], 
[0, 1, 0], 
[0, 0, 1], 
[0, 1, 1],  
[1, 1, 0], 
[1, 0, 1]
] 
)


H = np.array(\
[
[0, 1, 1, 1, 0, 0], 
[1, 1, 0, 0, 1, 0],
[1, 0, 1, 0, 0, 1] 
] 
) 

#INPUT MESSAGE HERE
message = np.array([1,0,1])
notcodeword = np.mod(np.dot(G,message.T), 2) 
# print(codeword)
check = np.mod(np.dot(H,notcodeword), 2) 
# print(check)
#------------------------------------------------------------------------------------- 

def find_error(codewort, error_sequences, parity): 
    for i, error in enumerate(error_sequences):
        iplus = i+1 
        
        if np.array_equal(np.mod(np.dot(parity, codewort - error), 2), np.zeros(len(message))):
            noerror = np.mod((codewort - error), 2)
            return f"# of tries: {iplus} originial codeword: {noerror}"

    return f"Number of tries: {len(error_sequences)} originial codeword: {codewort}"

codeword = np.mod((notcodeword + chosen_error), 2)
print(f"original codeword {notcodeword}") 
print(f"codeword with error {codeword}") 
grand = find_error(codeword, genarr_lst, H) 

print(grand)
