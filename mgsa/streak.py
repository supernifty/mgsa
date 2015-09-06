import sys
sys.setrecursionlimit(20000)

def probOfStreak(numCoins, minHeads, headProb, saved=None):
    if saved == None:
        saved = {}

    ID = (numCoins, minHeads, headProb)

    if ID in saved: 
        return saved[ID]
    else:
        if minHeads > numCoins or numCoins <= 0:
            result = 0
        else:
            result = headProb**minHeads
            for firstTail in xrange(1, minHeads+1):
                pr = probOfStreak(numCoins-firstTail, minHeads, headProb, saved)
                result += (headProb**(firstTail-1))*(1-headProb)*pr
        saved[ID] = result

        return result

print probOfStreak(100, 19, 0.90)
