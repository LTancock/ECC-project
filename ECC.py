import numpy as np
import math
import matplotlib

#errorCap is for testing only in scenarios when the number of errors is <= number of errors that can be corrected
#otherwise just put it to be 0, where it will allow any number of errors
def generateErrors(word, n, p, errorCap):
    #print(word)
    if errorCap == 0: errorCap = n
    #FIX, using rng makes them all output the same answer
    rng = np.random.default_rng(0)
    word2 = np.copy(word)
    #errors = rng.integers(low = 0, high = 2, size = n)
    #errors = np.random.random_integers(low = 0, high = 2, size = n)
    errors = np.random.binomial(1, p, n)
    #errors = rng.binomial(1, p, n)
    #print(word)
    for i in range(0, len(word)):
        word2[i] = (word[i] != errors[i]) #does xor
        if errors[i] == 1:
            errorCap = errorCap - 1
        if (errorCap == 0):
            break
    #print(word)
    return word2

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
#-1 = no error, -2 = unknown error
#when changed to multi-error detecting, make array of error bits so no error is an empty array
def synDecode(synTable, word):
    errorBit = -2
    if np.array_equal(word, [0, 0, 0]): errorBit = -1
    #print(word)
    for i in range(0, len(synTable)):
        if (np.array_equal(word, synTable[i])):
            errorBit = i
    return errorBit

def genToPc(genMat, n, k):
    #takes last n-k columns, transpose and set as rows, then put n-k identity matrix
    pcMat = np.transpose(genMat[:,range(n-k+1, n)])
    pcMat = np.hstack((pcMat, np.identity(n-k)))
    return pcMat

def pcToGen(pcMat, n, k):
    #takes last n-k columns, transpose and set as rows, then put n-k identity matrix
    genMat = np.transpose(pcMat[:,range(0, k)])
    genMat = np.hstack((np.identity(k), genMat))
    return genMat

#compute synTable before since this is run multiple times
def test(word, pcMat, synTable, n, d, p):
    #print(word)
    wordE = generateErrors(word, n, p, 0)
    #print(word)
    synWord = np.matmul(pcMat, wordE)%2
    #print(synWord)
    errorBit = synDecode(synTable, synWord)
    #print(errorBit)
    #print(wordE)
    if errorBit >= 0:
        wordE[errorBit] = (wordE[errorBit] != 1)
    #print(word)
    #print(wordE)
    return ([np.array_equal(wordE, word), errorBit])

def fullTest(genMat, n, k, d, p, tests):
    #run test tests times and return data for plotting
    word = genMat[0]
    pcMat = genToPc(genMat, n, k)
    synTable = generateSynTable(n,d,pcMat)
    data = np.empty(tests)
    for i in range(0, tests):
        data[i] = test(word, pcMat, synTable, n, d, p)[0]
    return data

def main():
    n = 7
    k = 4
    d = 3
    p = 0.3
    #word = np.array([1, 0, 0, 0, 1, 1, 0])
    #word = generateErrors(word, n, p, 1)
    genMat = np.array([[1, 0, 0, 0, 1, 1, 0], [0, 1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1, 1], [0, 0, 0, 1, 1, 1, 1]])
    data = fullTest(genMat, n, k, d, p, 10)
    print(data)
    #pcMat = np.array([[1, 1, 0, 1, 1, 0, 0], [1, 0, 1, 1, 0, 1, 0], [0, 1, 1, 1, 0, 0, 1]])
    #print(pcMat)
    #pcMat = genToPc(genMat, n, k)
    #print(pcMat)
    #synTable = generateSynTable(n, d, pcMat)
    #print(test(word, pcMat, synTable, n, d, p))
    #test2 = pcToGen(pcMat, n, k)
    #synWord = np.matmul(pcMat, word)%2
    #synTable = np.array([[1,1,0],[1,0,1],[0,1,1],[1,1,1],[1,0,0],[0,1,0],[0,0,1]])
    #errorBit = synDecode(synTable, synWord)
    #print(errorBit)


if __name__ == '__main__':
    main()