###########################################################
# Author: Jeffrey A. Aborot                               #
# Institution: Algorithms & Complexity Laboratory         #
#              Department of Computer Science             #
#              University of the Philippines Diliman      #
###########################################################

from math import log, ceil
import numpy as np
from numpy import int8

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

# Methods
def qRAM_alphabet(symbol):
    return alphabet.index(symbol)

def qRAM_text(text_index):
    return text[text_index]

def qRAM_pattern(symbol_index):
    try:
        return pattern.index(text[symbol_index])
    except ValueError:
        return -1 # return index -1 if symbol is not in P

# This method returns this decimal number's binary number format.
def convert_decimal_to_binary(decimal_number, register_size):
    # bit_list = []
    # if decimal_number < 0:
    #     # get absolute value of decimal number
    #     decimal_number = abs(decimal_number)
    #     # get binary representation of decimal number
    #     for x in xrange(0, register_size):
    #         i = register_size - x - 1
    #         if decimal_number >= pow(2, i):
    #             bit_list.append(1)
    #             decimal_number = decimal_number - pow(2, i)
    #         else:
    #             bit_list.append(0)
    #     bit_list.reverse()
    #     # invert binary number
    #     bit_list
    #     # add 1 to binary number
    #     # return resulting binary number
    # else:
    #     # get binary representation of decimal number
    #     # return resulting binary number
    return np.binary_repr(decimal_number, width=8)

# get the index of a symbol in alphabet in binary format
def convert_symbol_to_binary(symbol, register_size):
    alphabet_index = alphabet.index(symbol)
    return convert_decimal_to_binary(alphabet_index, register_size)


# Express the register state in binary format
def binarize(register):
    binary_register_state = [[0 for x in range(3)] for y in range(N)]
    for index in range(N):
        # index register
        binary_register_state[index][0] = convert_decimal_to_binary(register[index][0], int(ceil(log(N))))

        # symbol register
        binary_register_state[index][1] = convert_decimal_to_binary(register[index][1], int(ceil(log(alphabet_size))))

        # location register
        binary_register_state[index][2] = convert_decimal_to_binary(register[index][2], int(ceil(log(M)) + 1))
    return binary_register_state




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
register = [[0 for x in range(3)] for y in range(N)]
print '\n'
print 'Step 1 register state ============================'
print 'decimal:', register
print 'binary:', binarize(register)

# 2. Put into a superposition the index register.
for i in range(N):
    register[i][0] = i
print '\n'
print 'Step 2 register state ============================'
print 'decimal:', register
print 'binary:', binarize(register)

# 3. Put the symbol register into the state |Beta(T[i])> where Beta(.) is
#   	a function which maps a symbol in Sigma to an log Sigma-bit
#       binary number.
for i in range(N):
    register[i][1] = qRAM_alphabet(qRAM_text(i))
print '\n'
print 'Step 3 register state ============================'
print 'decimal:', register
print 'binary:', binarize(register)

# 4. Put the location register into the state |Psi(T[i])> where Psi(.) is
#       a function which maps a symbol in Sigma to its location of
#       first occurrence in P.
for i in range(N):
    register[i][2] = qRAM_pattern(register[i][1])
print '\n'
print 'Step 4 register state ============================'
print 'decimal:', register
print 'binary:', binarize(register)

# 5. Put the location register into the state |i-Psi(T[i])>.
for i in range(N):
    register[i][2] = i - register[i][2]
print '\n'
print 'Step 5 register state ============================'
print 'decimal:', register
print 'binary:', binarize(register)

# 6. Measure the state of the location register and output the classical
#       value. (pseudo-randomly select among the indices of T)

# TODO: Merge the remote branch main-method to remote branch master.
# TODO: Pull updates from remote branch master to local branch master.
