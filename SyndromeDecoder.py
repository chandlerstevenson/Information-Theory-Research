import numpy as np 
import random  
import sympy as sympy 
import itertools as listy   
import mpmath

generator_matrix = np.array([[1, 0, 0, 0, 1, 1],
                             [0, 1, 0, 1, 1, 0],
                             [0, 0, 1, 1, 0, 1]]) 

random_number = random.randint(0,1)

parity_check = np.array([[0, 1, 1, 1, 0, 0], 
                         [1, 1, 0, 0, 1, 0], 
                         [1, 0, 1, 0, 0, 1]])

# Let message be a 1x3 matrix
message = np.array([1, 0, 1]) 

# encode the message using the generator matrix 
# there's a risk that a number greater than 2 could appear  
# therefore, to conserve the intended values, mod 2 
c = np.mod(np.dot(generator_matrix.T,message), 2)



# Add an error to the codeword
codeword_with_error = np.mod(c + np.array([0, 0, 0, 1, 1, 1]), 2) 
print(f"This is the ouput of the error {codeword_with_error}")

# Decode the codeword using the parity-check matrix
syndrome = np.mod(np.dot(codeword_with_error.T, parity_check.T), 2) 
print(f"This is syndrome {syndrome}") 

# parity_checkmult = np.mod(np.dot(message, parity_check), 2)

# print(f"This should be zero {parity_checkmult}")

error_vector = np.zeros(len(parity_check[0])) #this makes an initial error vector with values of
 
for i in range(len(parity_check[0])):
    if np.array_equal(syndrome, parity_check[:,i]):
        error_vector[i] = 1
        break 




decoded_codeword = np.mod(codeword_with_error + error_vector, 2) 
final_decoded = decoded_codeword[0:3]


# Print the results
print("message:", message)
print("encoded codeword:", c)
print('Codeword with "noise":', codeword_with_error)
print("Decoded codeword:", final_decoded)

 