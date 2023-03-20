import numpy as np 

def fullmakepolygen(crc_poly): 

        def codeword_helper(crc_poly): 
            degree = len(crc_poly) - 1
            while degree >= 0 and crc_poly[degree] == 0:
                degree -= 1
            if degree < 0:
                return crc_poly + [0] * (len(crc_poly) + 1)
            else:
                return crc_poly + [0] * (degree) 
    
        made_codeword = codeword_helper(crc_poly)  #CODEWORD RESULT   

        def crc_polynomial_to_shifts(crc_poly):
            crc_poly = crc_poly[::-1]
            shifts = []
            for i in range(len(crc_poly)):
                if crc_poly[i] == 1:
                    shifts.append(len(crc_poly) - i - 1)
            return shifts[::-1]

        shift = crc_polynomial_to_shifts(made_codeword) 

        def generate_shift_matrix(shifts, input_list): 
            shifts = shifts
            shift_matrix = np.zeros((len(shifts), len(input_list)), dtype=int)
            for i, shift in enumerate(shifts):
                if shift >= len(input_list):
                    shift_matrix[i, :] = input_list
                else:
                    shift_matrix[i, shift:] = input_list[:-shift or None]
            return shift_matrix

        shift_matrix = generate_shift_matrix(shift, made_codeword)  
        
        return shift_matrix 

testing = fullmakepolygen([1,1,1,1])

print(testing)