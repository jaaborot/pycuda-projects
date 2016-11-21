###########################################################
# Author: Jeffrey A. Aborot                               #
# Institution: Algorithms & Complexity Laboratory         #
#              Department of Computer Science             #
#              University of the Philippines Diliman      #
###########################################################

from math import log, ceil, sqrt
import numpy as np
from numpy import int8, dot, kron, identity

 # Assumptions:
	#   1. A qRAM holds a copy of T and can be queried for any i-th element of T
	#   		by passing i as parameter value. qRAM puts its output register into
	#   		the state representing the i-th symbol of T. i.e. qRAM^T(i) = T[i]
	#   2. A qRAM holds data about the location of first occurrence in P of any
	#   		symbol in Sigma. i.e. qRAM^P(T[i]) = location of first occurrence of
	#   		symbol T[i] in P.
	#   3. Sigma = {a, c, t, g}
	#

#  Algo C Steps:
	#   1. Initialize.
	#   	1.1 index register with log N bits into state 0 [accomodates
	#   			integers 0,1,...,N-1]
	#   	1.2 symbol register with log |Sigma| bits into state 0 [accomodates
	#   			binary representation of symbols in Sigma]
	#   	1.3 location register with log |Sigma| + 1 bit into state 0
	#   			[accomodates integers -N,-(N-1),...,0,1,...,N-1]
	#   2. Put into a superposition the index register.
	#   3. Put the symbol register into the state |Beta(T[i])> where Beta(.) is
	#   			a function which maps a symbol in Sigma to an log Sigma-bit
	#   			binary number.
	#   4. Put the location register into the state |Psi(T[i])> where Psi(.) is
	#   			a function which maps a symbol in Sigma to its location of
	#   			first occurrence in P.
	#   5. Put the location register into the state |i-Psi(T[i])>.
	#   6. Measure the state of the location register and output the classical
	#   			value.

# Classes
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Kernels
def superpositioning_kernel():
    pass

def symbol_fetching_kernel():
    pass

def location_fectching_kernel():
    pass

def index_fetching_kernel():
    pass


# Methods
def superposition_matrix_operator():
    H = 1 / sqrt(2) * np.array([[1,1],[1,-1]])
    superpostion_operator_bit_count = _index_register_size - 1
    superposition_operator = H
    for bit in xrange(superpostion_operator_bit_count):
        superposition_operator = kron(superposition_operator, H)

    I = np.identity(2)
    identity_operator_bit_count = _symbol_register_size + _location_register_size - 1
    identity_operator = I
    for bit in xrange(identity_operator_bit_count):
        identity_operator = kron(identity_operator, I)

    operator = kron(superposition_operator, identity_operator)
    return operator

# The qRAM_text_matrix has two registers, input_register and output_register.
# input_register has size log(N) and holds the index i of the item to be fetched from the qRAM.
# output_register has size log(alphabet_size) and holds the binary code for symbol T[i].
# The dimension of the matrix is [2^{log(N)} x 2^{log(alphabet_size)}] x [2^{log(N)} x 2^{log(alphabet_size)}].
# The matrix is defined such that |i>|0> -> |i>|T[i]>.
def symbol_matrix_operator(text, alphabet):
    dimension = pow(2,_index_register_size) * pow(2,_symbol_register_size)
    matrix = np.identity(dimension)
    for index in xrange(N):
        symbol_decimal_code = qRAM_alphabet(text[index])  # decimal code of symbol in alphabet
        input_state = index * alphabet_size + 0
        output_state = index * alphabet_size + symbol_decimal_code
        matrix[input_state][input_state] = 0  # matrix[u][u] = 0
        matrix[output_state][output_state] = 0  # matrix[v][v] = 0
        matrix[input_state][output_state] = 1  # matrix[u][v] = 1
        matrix[output_state][input_state] = 1  # matrix[v][u] = 1

    I = np.identity(2)
    identity_operator_bit_count = _location_register_size - 1
    identity_operator = I
    for bit in xrange(identity_operator_bit_count):
        identity_operator = kron(identity_operator, I)

    operator = kron(matrix, identity_operator)
    return operator

# TODO: Verify correctness of method qRAM_pattern_matrix
def location_matrix_operator(pattern, alphabet):
    location_register_state_dimension = pow(2,_location_register_size)
    dimension = pow(2,_symbol_register_size) * location_register_state_dimension
    matrix = np.identity(dimension)
    for alphabet_index in range(alphabet_size):
        location_decimal_code = qRAM_location(alphabet[alphabet_index]) # location of first occurrence in P of a given symbol
        location_decimal_code = location_decimal_code + N # to compensate for negative values so that a location_decimal_code of -N will be mapped to the 0-th element of the pow(2, int(ceil(log(N) + 1))) values -N, -(N-1),...,-1,0,1,...,N-1
        input_state = alphabet_index * location_register_state_dimension + 0
        output_state = alphabet_index * location_register_state_dimension + location_decimal_code
        matrix[input_state][input_state] = 0  # matrix[u][u] = 0
        matrix[output_state][output_state] = 0  # matrix[v][v] = 0
        matrix[input_state][output_state] = 1  # matrix[u][v] = 1
        matrix[output_state][input_state] = 1  # matrix[v][u] = 1

    I = np.identity(2)
    identity_operator_bit_count = _index_register_size - 1
    identity_operator = I
    for bit in xrange(identity_operator_bit_count):
        identity_operator = kron(identity_operator, I)

    operator = kron(identity_operator, matrix)
    return operator

def subtraction_matrix_operator():
    total_register_size = _index_register_size + _symbol_register_size + _location_register_size
    dimension = pow(2,_index_register_size) * pow(2,_symbol_register_size) * pow(2,_location_register_size)
    matrix = np.zeros((dimension,dimension), dtype = float)
    print 'subtraction_matrix_operator: rows =', dimension, ', columns =', dimension
    for row in xrange(dimension):
        bin_row = format(row, '0'+str(total_register_size)+'b') # binary notation for input state |i>|T[i]>|G(T[i])>
        bin_index_symbol_component = bin_row[0 : _index_register_size + _symbol_register_size] # binary notation for index (|i>) and symbol (|T[i]>) register components of current row index
        bin_index_component = bin_row[0 : _index_register_size] # binary notation for index register (|i>) component of current row index
        bin_location_component = bin_row[_index_register_size + _symbol_register_size : total_register_size] # binary notation for location register (|G(T[i])>) component of current row index
        dec_index = int(bin_index_component, base=2) # decimal notation for index register (|i>) component of current row index
        dec_location = int(bin_location_component, base=2) # decimal notation for location register (|G(T[i])>) component of current row index
        dec_location = dec_index - dec_location # decimal notation for |i - G(T[i])>
        # bin_location_component = format(dec_location, '0'+str(_location_register_size)+'b') # binary notation for |i - G(T[i])>
        bin_location_component = np.binary_repr(dec_location, width=_location_register_size) # binary notation for |i - G(T[i])>
        bin_row = bin_index_symbol_component + bin_location_component # binary notation for output state |i>|T[i]>|i-G(T[i])>
        col = int(bin_row, base=2)
        print 'matrix[', row, '][', col, '] = 1.0'
        matrix[row][col] = 1.0
    return matrix


def qRAM_alphabet_matrix(alphabet):
    pass

def qRAM_alphabet(symbol):
    return alphabet.index(symbol)


def qRAM_text(text_index):
    return text[text_index]

def qRAM_location(symbol):
    try:
        return pattern.index(symbol)
    except ValueError:
        return -1

def qRAM_pattern(text_index):
    try:
        return pattern.index(text[text_index])
    except ValueError:
        return -1 # return index -1 if symbol is not in P


# This method returns this decimal number's binary number format.
def convert_decimal_to_binary(decimal_number, register_size):
    return np.binary_repr(decimal_number, width=register_size)


# Express the register state in binary format
def binarize(register):
    binary_register_state = [[0 for x in range(3)] for y in range(N)]
    for index in range(N):
        # index register
        binary_register_state[index][0] = convert_decimal_to_binary(register[index][0], int(ceil(log(N,2))))

        # symbol register
        binary_register_state[index][1] = convert_decimal_to_binary(register[index][1], int(ceil(log(alphabet_size,2))))

        # location register
        binary_register_state[index][2] = convert_decimal_to_binary(register[index][2], int(ceil(log(M,2)) + 1))
    return binary_register_state

# This method returns the vector representation of a binary number. The size of the
# vector is the 2^register_size.
def convert_decimal_to_vector(decimal_number, register_code, register_size):
    vector = [0 for x in range(pow(2,register_size))]

    # index register [unsigned integer, 0 to N-1, N vectors]
    # symbol register [unsigned integer, 0 to alphabet_size - 1, alphabet_size number of vectors]
    if register_code == _index_register or register_code == _symbol_register:
        vector[decimal_number] = 1

    # location register [unsigned integer, -(M-1) to 0 to +(N-1), N + M - 1 number of vectors]
    elif register_code == _location_register:
        vector[decimal_number + N] = 1

    return vector

# This method expresses the state of the registers as individual vectors per register
def vectorize(register):
    vector_register_state = [[np.empty([1,1]) for x in range(3)] for y in range(N)]
    for index in range(N):
        # index register
        vector_register_state[index][0] = convert_decimal_to_vector(register[index][0], _index_register, int(ceil(log(N,2))))

        # symbol register
        vector_register_state[index][1] = convert_decimal_to_vector(register[index][1], _symbol_register, int(ceil(log(alphabet_size,2))))

        # location register
        vector_register_state[index][2] = convert_decimal_to_vector(register[index][2], _location_register, int(ceil(log(N,2)) + 1))

    return vector_register_state

# This method expresses the state of the registers as a single vector only
def vectorize_compact(register):
    print 'vectorize_compact:'
    vector_register_state = [0 for x in xrange(N)]
    for index in xrange(N):
        # index register
        index_vector = np.array(convert_decimal_to_vector(register[index][0], _index_register, _index_register_size))[np.newaxis, :].T

        # symbol register
        symbol_vector = np.array(convert_decimal_to_vector(register[index][1], _symbol_register, _symbol_register_size))[np.newaxis, :].T

        # location register
        location_vector = np.array(convert_decimal_to_vector(register[index][2], _location_register, _location_register_size))[np.newaxis, :].T

        vector_register_state[index] = np.kron(index_vector, np.kron(symbol_vector, location_vector))
        print 'vector_register_state[', index, ']: rows =', vector_register_state[index].shape[0], ', columns = ', vector_register_state[index].shape[1]

    # TODO: verify correctness of transposition of vector
    print vector_register_state
    return np.array(vector_register_state)

print 'Booyah!'

# =================== User input ===================

# alphabet
alphabet = raw_input('Specify alphabet symbols in comma-separated format: ').split(',')
alphabet_size = len(alphabet)
print 'Alphabet: ', alphabet
print 'Alphabet size: ', len(alphabet)

# text
text = raw_input('Specify text: ')
N = len(text)
print 'T: ', text
print 'N: ', N

# pattern
pattern = raw_input('Specify pattern: ')
M = len(pattern)
print 'P: ', pattern
print 'M: ', len(pattern)

# constant variables
_index_register = 0
_symbol_register = 1
_location_register = 2
_index_register_size = int(ceil(log(N,2)))
_symbol_register_size = int(ceil(log(alphabet_size,2)))
_location_register_size = int(ceil(log(N,2))) + 1
print '_index_register_size =', _index_register_size
print '_symbol_register_size =', _symbol_register_size
print '_location_register_size =', _location_register_size

# =================== Prepare qRAMs ===================
# qRAM_text: index in T -> symbol
# index = int(raw_input('Enter index of symbol in T: '))
# symbol = qRAM_text(index)
# print 'Symbol in index', index, 'in T:', symbol

# qRAM_pattern: symbol -> location of first occurrence of symbol in P
# location = qRAM_pattern(symbol)
# print 'Location of first occurrence of symbol', symbol, 'in P:', location

#  =================== Define registers  ===================
# 1. =================== Initialize ===================
# 1.1 Initialize index register into state |0>
# 1.2 Initialize symbol register into state |0>
# 1.3 Initialize symbol register into state |0>
register = [[0 for x in xrange(3)] for y in xrange(N)]

print '\n'
print 'Step 1 register state ============================'
print 'decimal:', register
print 'binary:', binarize(register)
print 'vector:', vectorize(register)
zero_state = vectorize_compact(register)
zero_state_rows = zero_state.shape[0]
zero_state_cols = zero_state.shape[1]
print 'zero_state (rows = ', zero_state_rows,', columns, ', zero_state_cols,'):'
print 'Register state in compact vector notation:'
for row in xrange(zero_state_rows):
    x_dim = zero_state[row].shape[0]
    y_dim = zero_state[row].shape[1]
    print 'zero_state[',row,']: rows =',x_dim,'columns =',y_dim
    print zero_state[row]

# 2. Put into a superposition the index register.
superposition_operator = superposition_matrix_operator()
superposition_operator_rows = superposition_operator.shape[0]
superposition_operator_cols = superposition_operator.shape[1]
print 'Superposition operator (rows = ', superposition_operator_rows, ', columns = ', superposition_operator_cols,'):'
# print '['
# for row in xrange(superposition_operator_rows):
#     print '[',
#     for col in xrange(superposition_operator_cols):
#         if superposition_operator[row][col]:
#             print bcolors.FAIL + str(superposition_operator[row][col]) + bcolors.ENDC,
#         else:
#             print superposition_operator[row][col],
#     print '],'
# print ']'
superpositioned_state = []
for vector in zero_state:
    superpositioned_state.append(dot(superposition_operator, vector))
superpositioned_state_rows = np.array(superpositioned_state).shape[0]
superpositioned_state_cols = np.array(superpositioned_state).shape[1]
for index in xrange(superpositioned_state_rows):
    print 'superpositioned_state[',index,'] (rows = ', np.array(superpositioned_state[index]).shape[0],', columns = ', np.array(superpositioned_state[index]).shape[1],'):'
    print superpositioned_state[index]
# print '['
# for row in xrange(superpositioned_state_rows):
#     print '[',
#     for col in xrange(superpositioned_state_cols):
#         if superpositioned_state[row][col]:
#             print bcolors.FAIL + str(superpositioned_state[row][col]) + bcolors.ENDC,
#         else:
#             print superpositioned_state[row][col],
#     print '],'
# print ']'
#
# for i in range(N):
#     register[i][0] = i
# print '\n'
# print 'Step 2 register state ============================'
# print 'decimal:', register
# print 'binary:', binarize(register)
# print 'vector:', vectorize(register)
# vectorized_compact_register = vectorize_compact(register)
# vectorized_compact_register_rows = vectorized_compact_register.shape[0]
# vectorized_compact_register_cols = vectorized_compact_register.shape[1]
# print 'Register state in compact vector notation:'
# print '['
# for row in xrange(vectorized_compact_register_rows):
#     print '[',
#     for col in xrange(vectorized_compact_register_cols):
#         if vectorized_compact_register[row][col]:
#             print bcolors.FAIL + str(vectorized_compact_register[row][col]) + bcolors.ENDC,
#         else:
#             print vectorized_compact_register[row][col],
#     print '],'
# print ']'

# 3. Put the symbol register into the state |Beta(T[i])> where Beta(.) is
#   	a function which maps a symbol in Sigma to an log Sigma-bit
#       binary number.
symbol_operator = symbol_matrix_operator(text, alphabet)
symbol_operator_rows = symbol_operator.shape[0]
symbol_operator_cols = symbol_operator.shape[1]
print 'symbol_operator: rows = ', symbol_operator.shape[0], ', columns = ', symbol_operator.shape[1]
symbol_state = []
for vector in superpositioned_state:
    symbol_state.append(dot(symbol_operator, vector))
symbol_state_rows = np.array(symbol_state).shape[0]
symbol_state_cols = np.array(symbol_state).shape[1]
for index in xrange(symbol_state_rows):
    print 'symbol_state[',index,']: rows =', np.array(symbol_state[index]).shape[0], ',columns =', np.array(symbol_state[index]).shape[1]
    print symbol_state[index]


# for i in range(N):
#     register[i][1] = qRAM_alphabet(qRAM_text(i))
# print '\n'
# print 'Step 3 register state ============================'
# print 'decimal:', register
# print 'binary:', binarize(register)
# print 'vector:', vectorize(register)
# vectorized_compact_register = vectorize_compact(register)
# vectorized_compact_register_rows = vectorized_compact_register.shape[0]
# vectorized_compact_register_cols = vectorized_compact_register.shape[1]
# print 'Register state in compact vector notation:'
# print '['
# for row in xrange(vectorized_compact_register_rows):
#     print '[',
#     for col in xrange(vectorized_compact_register_cols):
#         if vectorized_compact_register[row][col]:
#             print bcolors.FAIL + str(vectorized_compact_register[row][col]) + bcolors.ENDC,
#         else:
#             print vectorized_compact_register[row][col],
#     print '],'
# print ']'

# 4. Put the location register into the state |Psi(T[i])> where Psi(.) is
#       a function which maps a symbol in Sigma to its location of
#       first occurrence in P.
location_operator = location_matrix_operator(pattern, alphabet)
location_operator_rows = location_operator.shape[0]
location_operator_cols = location_operator.shape[1]
print 'location_operator: rows = ', location_operator.shape[0], ', columns = ', location_operator.shape[1]
location_state = []
for vector in symbol_state:
    location_state.append(dot(location_operator, vector))
location_state_rows = np.array(location_state).shape[0]
location_state_cols = np.array(location_state).shape[1]
for index in xrange(location_state_rows):
    print 'location_state[',index,']: rows =', np.array(location_state[index]).shape[0], ',columns =', np.array(location_state[index]).shape[1]
    print location_state[index]

# for i in range(N):
#     register[i][2] = qRAM_pattern(register[i][1])
# print '\n'
# print 'Step 4 register state ============================'
# print 'decimal:', register
# print 'binary:', binarize(register)
# print 'vector:', vectorize(register)
# vectorized_compact_register = vectorize_compact(register)
# vectorized_compact_register_rows = vectorized_compact_register.shape[0]
# vectorized_compact_register_cols = vectorized_compact_register.shape[1]
# print 'Register state in compact vector notation:'
# print '['
# for row in xrange(vectorized_compact_register_rows):
#     print '[',
#     for col in xrange(vectorized_compact_register_cols):
#         if vectorized_compact_register[row][col]:
#             print bcolors.FAIL + str(vectorized_compact_register[row][col]) + bcolors.ENDC,
#         else:
#             print vectorized_compact_register[row][col],
#     print '],'
# print ']'

# 5. Put the location register into the state |i-Psi(T[i])>.
subtraction_operator = subtraction_matrix_operator()
subtraction_operator_rows = subtraction_operator.shape[0]
subtraction_operator_cols = subtraction_operator.shape[1]
print 'subtraction_operator: rows =', subtraction_operator_rows, ', columns =', subtraction_operator_cols
subtraction_state = []
for vector in location_state:
    subtraction_state.append(dot(subtraction_operator, vector))
subtraction_state_rows = np.array(subtraction_state).shape[0]
subtraction_state_cols = np.array(subtraction_state).shape[1]
for index in xrange(subtraction_state_rows):
    print 'subtraction_state[',index,']: rows =', np.array(subtraction_state[index]).shape[0], ',columns =', np.array(subtraction_state[index]).shape[1]
    print subtraction_state[index]

# for i in range(N):
#     register[i][2] = i - register[i][2]
# print '\n'
# print 'Step 5 register state ============================'
# print 'decimal:', register
# print 'binary:', binarize(register)
# print 'vector:', vectorize(register)
# vectorized_compact_register = vectorize_compact(register)
# vectorized_compact_register_rows = vectorized_compact_register.shape[0]
# vectorized_compact_register_cols = vectorized_compact_register.shape[1]
# print 'Register state in compact vector notation:'
# print '['
# for row in xrange(vectorized_compact_register_rows):
#     print '[',
#     for col in xrange(vectorized_compact_register_cols):
#         if vectorized_compact_register[row][col]:
#             print bcolors.FAIL + str(vectorized_compact_register[row][col]) + bcolors.ENDC,
#         else:
#             print vectorized_compact_register[row][col],
#     print '],'
# print ']'

# print 'qRAM_text_matrix:\n'
# qRAM_text_matrix_operator = qRAM_text_matrix(text,alphabet)
# print 'qRAM_text_matrix_operator.shape =', qRAM_text_matrix_operator.shape
# qRAM_text_matrix_operator_rows = qRAM_text_matrix_operator.shape[0]
# qRAM_text_matrix_operator_cols = qRAM_text_matrix_operator.shape[1]
# for row in xrange(qRAM_text_matrix_operator_rows):
#     for col in xrange(qRAM_text_matrix_operator_cols):
#         if int(qRAM_text_matrix_operator[row][col]):
#             print bcolors.WARNING + str(int(qRAM_text_matrix_operator[row][col])) + bcolors.ENDC,
#         else:
#             print int(qRAM_text_matrix_operator[row][col]),
#     print '\n'


# print 'qRAM_pattern_matrix:\n'
# qRAM_pattern_matrix_operator = qRAM_pattern_matrix(pattern,alphabet)
# print 'qRAM_pattern_matrix_operator.shape = ', qRAM_pattern_matrix_operator.shape
# qRAM_pattern_matrix_operator_rows = qRAM_pattern_matrix_operator.shape[0]
# qRAM_pattern_matrix_operator_cols = qRAM_pattern_matrix_operator.shape[1]
# for row in xrange(qRAM_pattern_matrix_operator_rows):
#     for col in xrange(qRAM_pattern_matrix_operator_cols):
#         if int(qRAM_pattern_matrix_operator[row][col]):
#             print bcolors.WARNING + str(int(qRAM_pattern_matrix_operator[row][col])) + bcolors.ENDC,
#         else:
#             print int(qRAM_pattern_matrix_operator[row][col]),
#     print '\n'

# 6. Measure the state of the location register and output the classical
#       value. (pseudo-randomly select among the indices of T)

# TODO: Merge the remote branch main-method to remote branch master.
# TODO: Pull updates from remote branch master to local branch master.
