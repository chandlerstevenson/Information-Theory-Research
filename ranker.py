import numpy as np 
import math 
def rank_errors(channel_prob, error_size=3):
    # here, we generate all possible error sequences of a given size 
    # it follows that because this is a bernoulli distribution, we should see 
    #    a sequence of 2**n where n is the number of bits in the sequence 
    # I: tuple -> (channel_prob, error_size) => (float, int) 
        #RI: Map integers in the sequence to strings to store in a dictionary  
        # RO: Sort said dictionary by thte the probability of each occurence.  
    # O: 
    error_sequences = generate_error_sequences(error_size)

    # calculate the probability of each error sequence
    error_probs = {}
    for error_seq in error_sequences:
        error_prob = channel_prob ** sum(error_seq) * (1 - channel_prob) ** (error_size - sum(error_seq))
        error_probs[''.join(map(str, error_seq))] = error_prob

    # sort error probabilities in descending order
    sorted_errors = sorted(error_probs.items(), key=lambda x: x[1], reverse=True)
    testkey = list(sorted_errors)[0:4] 
    
    # print the most likely error sequences and their probabilities
    for error, prob in sorted_errors:
        print(f"Error Sequence: {error}, Probability: {prob}")


def generate_error_sequences(length_codeword):
    # generate all possible binary sequences of a given size
    sequences = []
    for i in range(2 ** length_codeword):
        seq = [int(bin_num) for bin_num in bin(i)[2:].zfill(length_codeword)]
        sequences.append(seq)
    return sequences 


'''
entropy: measures the amount of 'surprise' in the set of given binary digits.  

'''
def entropy(sequence): 
    # count the number of occurrences of each symbol
    counts = {}
    for symbol in sequence:
        if symbol not in counts:
            counts[symbol] = 0
        counts[symbol] += 1

    # calculate the probability of each symbol
    total_count = len(sequence)
    probabilities = {symbol: count / total_count for symbol, count in counts.items()}

    # calculate the entropy rate
    entropy_rate = 0
    for symbol, probability in probabilities.items():
        entropy_rate = ((-1*probability * math.log2(probability)) - (1-probability)*math.log2(1- probability))

    return f"entropy: {entropy_rate}, probability: {probabilities}, length: {total_count}"

#Example  
ranked = rank_errors(.2, 20) 
print(ranked)

# Example 
sequence = np.array([0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
result = entropy(sequence)
print(result)