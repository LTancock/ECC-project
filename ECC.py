import numpy as np
import math
from matplotlib import pyplot as plt

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

    
    #print(word)
    #print(wordE)
    if errorBit >= 0:
        wordE[errorBit] = (wordE[errorBit] != 1)
    
    return ([np.array_equal(wordE, word), errorBit])

def move(msg, dir):
    newMsg = np.copy(msg)
    #print("move first msg",newMsg)
    #print("move dir", dir)
    if (dir):
        i = len(newMsg) - 1
        newMsg[i] = 0 if newMsg[i] == 1 else 1
        #carry happened
        while(newMsg[i] == 0):
            if i > 0:
                newMsg[i-1] = 0 if newMsg[i-1] == 1 else 1
            else:
                break
            i = i - 1
    else:
        i = len(newMsg) - 1
        newMsg[i] = 0 if newMsg[i] == 1 else 1
        
        #carry happened
        while(newMsg[i] == 1):
            if i > 0:
                newMsg[i-1] = 0 if newMsg[i-1] == 1 else 1
            else:
                break
            i = i - 1
    #print("move second message", newMsg)
    return newMsg
                        

#models random walk, with [0,0,0,0] -> [1,1,1,1] when going down and vice-versa
def nextMsg(msg, maxStep):
    #picks 0 = down or 1 = up
    dir = np.random.randint(0,2)
    newMsg = move(msg, dir)
    return newMsg

def fullTest(genMat, n, k, d, p, tests):
    #run test and return data for plotting
    msg = [0,0,0,0]
    word = np.matmul(msg, genMat)%2
    pcMat = genToPc(genMat, n, k)
    synTable = generateSynTable(n,d,pcMat)
    data = np.empty(tests)
    numErrors = tests
    for i in range(0, tests):
        data[i] = test(word, pcMat, synTable, n, d, p)[0]
        msg = nextMsg(msg, 1)
        word = np.matmul(msg, genMat)%2
        numErrors -= data[i]
    print(numErrors) 
    return data

#word a and b should have same length
def hammingDist(worda, wordb):
    dist = 0
    for i in range(0, len(worda)):
        if worda[i] != wordb[i]: dist += 1
    return dist

#experiment with looking at before error correction and after
#deprecated
def markovTest(word, wordE, oldMsg, genMat, pcMat, synTable, n, d, p):
    #print("oldMsg", oldMsg)
    downMsg = move(oldMsg, 0)
    #print("oldMsg", oldMsg)
    #print("d", downMsg)
    upMsg = move(oldMsg, 1)
    #print("oldMsg", oldMsg)
    #print("u", upMsg)
    downWord = np.matmul(downMsg, genMat)%2
    upWord = np.matmul(upMsg, genMat)%2
    wordE = generateErrors(word, n, p, 0)
    synWord = np.matmul(pcMat, wordE)%2
    errorBit = synDecode(synTable, synWord)
    if errorBit >= 0:
        wordE[errorBit] = (wordE[errorBit] != 1)
    corWord = downWord
    if np.array_equal(wordE, upWord): corWord = upWord
    downDist = hammingDist(wordE, downWord)
    upDist = hammingDist(wordE, upWord)
    #return ([np.array_equal(wordE, word), np.array_equal(corWord, word), wordE])
    return ([np.array_equal(wordE, word), np.array_equal(corWord, word), np.array_equal(upWord, wordE) if upDist < downDist else np.array_equal(downWord, wordE), wordE])

#generate error outside and pass in
#separate oldMsg to collateTest
def MBCTest(word, oldMsg, genMat, n, p):
    downMsg = move(oldMsg, 0)
    upMsg = move(oldMsg, 1)
    downWord = np.matmul(downMsg, genMat)%2
    upWord = np.matmul(upMsg, genMat)%2
    wordE = generateErrors(word, n, p, 0)
    corWord = downWord
    downDist = hammingDist(wordE, downWord)
    upDist = hammingDist(wordE, upWord)
    if (upDist < downDist):
        corWord = upWord
    return ([np.array_equal(word, corWord), corWord])

#if MBC result used k times in a row trust hamming for next word
def collateTest(word, oldMsg, genMat, pcMat, synTable, n, p):
    downMsg = move(oldMsg, 0)
    upMsg = move(oldMsg, 1)
    downWord = np.matmul(downMsg, genMat)%2
    upWord = np.matmul(upMsg, genMat)%2
    wordE = generateErrors(word, n, p, 0)
    synWord = np.matmul(pcMat, wordE)%2
    errorBit = synDecode(synTable, synWord)
    downDist = hammingDist(wordE, downWord)
    upDist = hammingDist(wordE, upWord)
    corWord = downWord
    if upDist < downDist: corWord = upWord
    if errorBit >= 0:
        wordE[errorBit] = (wordE[errorBit] != 1)
    ham = False
    if np.array_equal(wordE, upWord) or np.array_equal(wordE, downWord):
        ham = True
    #print(word)
    #print(wordE)
    #print(corWord)
    #print(ham)
    #print(" ")
    return ([np.array_equal(word, wordE), np.array_equal(word, wordE) if ham else np.array_equal(word, corWord), wordE if ham else corWord, ham])
    #return hamming = np.array_equal(word, wordE), coll = if wordE == upWord or downWord then np.array_equal(word, wordE)
            #else np.array_equal(word, corWord), oldMsg = if wordE == upWord or downWord then wordE else corWord
    
def hammingTest(word, pcMat, synTable, n, p):
    wordE = generateErrors(word, n, p, 0)
    synWord = np.matmul(pcMat, wordE)%2
    errorBit = synDecode(synTable, synWord)
    if errorBit >= 0:
        wordE[errorBit] = (wordE[errorBit] != 1)
    return ([np.array_equal(word, wordE), wordE])

#run markov test and return data for plotting
#currently using correct previous msg rather than decoded (not necessarily correct) previous message
def fullMarkovTest(genMat, n, k, d, p, tests, MBCLim):
    msg = [0,0,0,0]
    word = np.matmul(msg, genMat)%2
    #wordE = generateErrors(word, n, p, 0)
    pcMat = genToPc(genMat, n, k)
    synTable = generateSynTable(n,d,pcMat)
    data = np.empty(tests)
    #assistData = np.empty(tests)
    mData = np.empty(tests)
    cData = np.empty(tests)
    numSucc = 0
    numMSucc = 0
    numCSucc = 0
    oldMsg = msg
    oldMsgC = msg
    oldWord = word
    oldWordC = word
    #limit for how many MBCs will be trusted in a row before it will trust Hamming to reset bad string
    MBC = 0
    for i in range(0, tests):
        #print(msg)
        #print(word)
        #data[i] = test(word, pcMat, synTable, n, d, p)[0]
        if i == 0:
            temp = hammingTest(word, pcMat, synTable, n, p)
            data[i] = mData[i] = cData[i] = temp[0]
            oldWord = oldWordC = temp[1]
            #data[i] = mData[i] = test(word, pcMat, synTable, n, d, p)[0]
            #mData[i] = test(word, pcMat, synTable, n, d, p)[0]
        else: 
            #Passing in error word from outside makes it bad???
            mData[i], oldWord = MBCTest(word, oldMsg, genMat, n, p)
            if (MBC < MBCLim):
                data[i], cData[i], oldWordC, ham = collateTest(word, oldMsgC, genMat, pcMat, synTable, n, p)
                if ham:
                    MBC = 0
                else:
                    MBC += 1
            else:
                temp = hammingTest(word, pcMat, synTable, n, p)
                data[i] = cData[i] = temp[0]
                oldWordC = temp[1]
                MBC = 0
        #oldWord = mData[2]

        oldMsg = oldWord[0:4]
        oldMsgC = oldWordC[0:4]
        #print(oldWord)
        #print(oldMsg)
        msg = nextMsg(msg, 1)
        word = np.matmul(msg, genMat)%2
        numSucc += data[i]
        numMSucc += mData[i]
        numCSucc += cData[i]
    print("Hamming s", numSucc)
    print("MBC s", numMSucc)
    print("Combined s", numCSucc)
    return ([data, mData, cData, numSucc, numMSucc, numCSucc])

#success prob vs prob of error, line graph with each of the 3
#success prob of combined vs MBCLim
def plot_succ_v_err(data, mData, cData, p_arr, tests):
    #print(numStep)
    #for i in range(0, len(data)):
    #    data[i] = data[i]/tests
    #    mData[i] = mData[i]/tests
    #    cData[i] = cData[i]/tests
    #x = np.empty(numStep)
    #for i in range(0, numStep):
    #    x[i] = start + (step * i)
    plt.title("Success probability vs probability of error for each bit")
    plt.plot(p_arr, data/tests, c = "r", marker = "o", label = "Hamming code")
    plt.plot(p_arr, mData/tests, c = "b", marker = "o", label = "MBC algorithm")
    plt.plot(p_arr, cData/tests, c = "g", marker = "o", label = "Combined algorithm")
    plt.xlabel("Probability of error for each bit")
    plt.ylabel("Probability for the algorithm to match the sent word")
    plt.legend(loc="upper right")
    plt.show()

def plot_noSucc_v_err(data, mData, cData, p_arr):
    
    plt.title("Number of successes vs probability of error for each bit")
    plt.plot(p_arr, data, c = "r", marker = "o", label = "Hamming code")
    plt.plot(p_arr, mData, c = "b", marker = "o", label = "MBC algorithm")
    plt.plot(p_arr, cData, c = "g", marker = "o", label = "Combined algorithm")
    plt.xlabel("Probability of error for each bit")
    plt.ylabel("Number of times the algorithm output was correct")
    plt.legend(loc="lower right")
    plt.show()

#prob 0.05
def plot_succ_v_MBCLim(data, cData, lims, tests):
    plt.title("Success probability vs value of MBCLim")
    plt.plot(lims, data/tests, c = "r", marker = "o", label = "Hamming code")
    plt.plot(lims, cData/tests, c = "g", marker = "o", label = "Combined algorithm")
    plt.xlabel("Number of times MBC is relied on before trusting hamming code")
    plt.ylabel("Probability for the algorithm to match the sent word")
    plt.legend(loc = "lower right")
    plt.show()


#Hamming code and contextual decoding implemented
#Hamming code is better, but we can combine the 2 since we know when we aren't correct (>1 error)
#Therefore, we can theoretically use the hamming code, and if there is more than 2 errors then we can use the contextualised
def main():
    n = 7
    k = 4
    d = 3
    p = 0.1
    MBCLim = 2
    size = 20
    p_arr = np.empty(size)
    MBCLim_arr = np.empty(size)
    for i in range (0, len(p_arr)):
        p_arr[i] = i * 0.05
        MBCLim_arr[i] = i
    tests = 100
    #print(MBCLim_arr[0], MBCLim_arr[2])
    #word = np.array([1, 0, 0, 0, 1, 1, 0])
    #word = generateErrors(word, n, p, 1)
    genMat = np.array([[1, 0, 0, 0, 1, 1, 0], [0, 1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1, 1], [0, 0, 0, 1, 1, 1, 1]])
    data = np.empty(size)
    mData = np.empty(size)
    cData = np.empty(size)
    
    for i in range(0, len(data)):
        rawdata, rawmData, rawcData, data[i], mData[i], cData[i] = fullMarkovTest(genMat, n, k, d, p_arr[i], tests, MBCLim)
    plot_succ_v_err(data, mData, cData, p_arr, tests)
    #plot_noSucc_v_err(data, mData, cData, p_arr)

    #for i in range(0, len(data)):
    #    rawdata, rawmData, rawcData, data[i], mData[i], cData[i] = fullMarkovTest(genMat, n, k, d, p, tests, MBCLim_arr[i])
    #plot_succ_v_MBCLim(data, cData, MBCLim_arr, tests)

    #print(data)
    #print(mData)
    #print(cData)
    #counter = 0
    #for i in range(0, tests-1):
    #    if(data[i] == 0):
    #        if(mData[i]==1): counter += 1
    #print(counter)


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