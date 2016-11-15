###########################################################
# Author: Jeffrey A. Aborot                               #
# Institution: Algorithms & Complexity Laboratory         #
#              Department of Computer Science             #
#              University of the Philippines Diliman      #
###########################################################

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

def qRAM_text(text_index):
    return text[text_index]

def qRAM_pattern(symbol):
    try:
        return pattern.index(symbol)
    except ValueError:
        return -1 # return index -1 if symbol is not in P

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
print 'register state:', register

# 2. Put into a superposition the index register.
for i in range(N):
    register[i][0] = i
print 'register state:', register

# 3. Put the symbol register into the state |Beta(T[i])> where Beta(.) is
#   	a function which maps a symbol in Sigma to an log Sigma-bit
#       binary number.
for i in range(N):
    register[i][1] = qRAM_text(i)
print 'register state:', register

# 4. Put the location register into the state |Psi(T[i])> where Psi(.) is
#       a function which maps a symbol in Sigma to its location of
#       first occurrence in P.
for i in range(N):
    register[i][2] = qRAM_pattern(register[i][1])
print 'register state:', register

# 5. Put the location register into the state |i-Psi(T[i])>.
for i in range(N):
    register[i][2] = i - register[i][2]
print 'register state:', register

# 6. Measure the state of the location register and output the classical
#       value. (pseudo-randomly select among the indices of T)