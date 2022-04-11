import numpy as np
import math

#errorCap is for testing only in scenarios when the number of errors is <= number of errors that can be corrected
#otherwise just put it to be 0, where it will allow any number of errors
def generateErrors(word, n, p, errorCap):
    if errorCap == 0: errorCap = n
    rng = np.random.default_rng(0)
    #errors = rng.integers(low = 0, high = 2, size = n)
    #errors = np.random.random_integers(low = 0, high = 2, size = n)
    #errors = np.random.binomial(1, p, n)
    errors = rng.binomial(1, p, n)
    print(word)
    for i in range(0, len(word)):
        word[i] = (word[i] != errors[i]) #does xor
        if errors[i] == 1:
            errorCap = errorCap - 1
        if (errorCap == 0):
            break
    print(word)
    return word

#1 error correcting atm
#talk about theory of generating the syndrome table in the report
def generateSynTable(n, d, pcMat):
    word = np.zeros(n)
    synTable = np.empty((0, d), int)#is it true that the length of the elements of the syn table is the min dist?
    for i in range(0, len(word)):
        word[i] = 1
        synTable = np.append(synTable, np.array([np.matmul(pcMat, word)%2]), axis = 0)
        word[i] = 0
    return synTable

#1 error correcting atm
#-1 = no error, -2 = unknown error maybe later
def synDecode(synTable, word):
    errorBit = -1
    print(word)
    for i in range(0, len(synTable)):
        if (np.array_equal(word, synTable[i])):
            errorBit = i
    return errorBit

def main():
    n = 7
    k = 4
    d = 3
    word = np.array([1, 0, 0, 0, 1, 1, 0])
    word = generateErrors(word, n, 0.2, 1)
    pcMat = np.array([[1, 1, 0, 1, 1, 0, 0], [1, 0, 1, 1, 0, 1, 0], [0, 1, 1, 1, 0, 0, 1]])
    synWord = np.matmul(pcMat, word)%2
    #synTable = np.array([[1,1,0],[1,0,1],[0,1,1],[1,1,1],[1,0,0],[0,1,0],[0,0,1]])
    synTable = generateSynTable(n, d, pcMat)
    errorBit = synDecode(synTable, synWord)
    print(errorBit)


if __name__ == '__main__':
    main()