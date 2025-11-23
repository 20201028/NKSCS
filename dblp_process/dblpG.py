import xml.etree.ElementTree as ET
import collections
import re
import heapq
# from nltk.stem.porter import *
# from nltk.stem import WordNetLemmatizer
graph_path = "D:\\Desktop\\CS\\dataProcess\\dblp.txt"
# txt_path = "D:\\Desktop\\CS\\dataProcess\\dblp-txt.txt"
# txt_path = "D:\\Desktop\\CS\\dataProcess\\sub.txt"
txt_path = "D:\\Desktop\\CS\\dataProcess\\dblp-2024-06-02-txt.txt"
stem_path = "D:\\Desktop\\CS\\dataProcess\\stemmer.lowercase.txt"
stop_path = "D:\\Desktop\\CS\\dataProcess\\stopword.txt"
att_path = "D:\\Desktop\\CS\\dataProcess\\dblpAtt.txt"
graphInf_path = "D:\\Desktop\\CS\\dataProcess\\dblpAInf.txt"

kwFreq = 0.2
kw2Id = {}
id2Kw = {}
# read every author's paper and cor-author   
userMap = {}
user2Title = {}
user2CorAuthor = {}
try:
    with open(txt_path, 'r') as f:
        count = 0
        title = None
        for line in f:
            count += 1
            line = line.strip()
                
            if count % 2 == 1:
                title = line
                if title.endswith('.'):
                    title = title[:-1]
            else:
                authors = line.split("\t")
                list_ = []
                    
                for author in authors:
                    if author not in userMap:
                        userId = len(userMap)
                        userMap[author] = userId
                        user2Title[userId]=set()
                        user2CorAuthor[userId]=set()
                        list_.append(userId)
                    else:
                        list_.append(userMap[author])
                    
                for i in range(len(list_)):
                    userId = list_[i]
                    user2Title[userId].add(title)
                    for j in list_:
                        if j != userId:
                            user2CorAuthor[userId].add(j)
                    # user2CorAuthor[userId].update(list_)
    
except Exception as e:
    print(e)
# read stemmer and stopword
def stemmer(word,stem):
    if word in stem:
        word = stem[word]
    return word
stem = {}
stop = set()
try:
    with open(stem_path, 'r') as f:
        for line in f:
            s = line.strip().split("/")
            for i in range(1, len(s)):
                stem[s[i]] = s[0]
    
except Exception as e:
    print(e)

try:
    with open(stop_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                stop.add(line)
    
except Exception as e:
    print(e)

# 

count = 0
vAttSum = 0
# Step 1: Count the frequency
# request each author's paper at least 3
for author, titles in user2Title.items():
    if len(titles) < 3: 
        user2Title[author].clear()
        continue
    kwMap = {}
    for title in titles:
        s = re.split(r'[^a-zA-Z]', title.lower())
        s = [w for w in s if w!='based' and w!='using']
        # s = title.split(' ')
        # stemmer = PorterStemmer()
        stemWord = [stemmer(w,stem) for w in s if w]
        stemSet = (ss for ss in stemWord)
        for word in stemSet:          
            if word not in stop:
                if word in kwMap:
                   kwMap[word] += 1
                else:
                    kwMap[word] = 1
    # Step 2: Extract the most frequent keywords
    # top = int(len(kwMap) * kwFreq)
    top = 20
    sorted_items = sorted(kwMap.items(), key=lambda x: x[1], reverse=True)
    
    top_k_strings = []
    if len(kwMap) > top:
        top_k_strings = sorted_items[:top]
    else:
        top_k_strings = sorted_items
    newTitle = set()
    vAttSum += len(top_k_strings)
    for string, c in top_k_strings:
        if string not in kw2Id:
            kw2Id[string] = str(count)
            id2Kw[str(count)] = string
            count+=1
        newTitle.add(kw2Id[string])

    # newTitle = set()
    # fre = 0
    # if sorted_items:
    #     lastFre = sorted_items[0][1]
    #     fre = lastFre*1
    # for string, c in sorted_items:
    #     if c < fre or c<2: break
    #     if string not in kw2Id:
    #         kw2Id[string] = str(count)
    #         id2Kw[str(count)] = string
    #         count+=1
    #     newTitle.add(kw2Id[string])



    user2Title[author] = newTitle

try:
    with open(graphInf_path, 'w') as stdout:
        stdout.write('node: ' + str(len(userMap)) + '\tatt:' + str(count) + "\tvAttSum: " + str(vAttSum) + "\n")
        for user, id in userMap.items():
            line = str(id) + "\t" + str(user)
            
            for word in user2Title[id]:
                line += "\t" + word + '\t' + id2Kw[word]
            
            stdout.write(line + "\n")
except Exception as e:
    print(e)



userNum = 0
edgeNum = 0
exuser=set()
try:
    with open(att_path, 'w') as stdout:
        for user, kw in user2Title.items():
            
            line = str(user) + '\t'   
            
            if kw:
                line += ",".join(kw)  
                userNum += 1
            else:
                # line += str(-1) 
                exuser.add(user)
                continue
            # for i in range(len(kw)-1):
            #     line += kw[i]  + ","
            # line+=kw[-1]
            stdout.write(line  + "\n")
except Exception as e:
    print(e)

try:
    with open(graph_path, 'w') as stdout:
        for user, nei in user2CorAuthor.items(): 
            if(user in exuser): continue   
            neis = nei.difference(exuser)
            for neighbor in neis:
                if(neighbor<user): continue
                edgeNum += 1
                line = str(user) + "\t" + str(neighbor)
                stdout.write(line  + "\n")
        stdout.write('#vNum: '+str(userNum)  + "\n")
        stdout.write('#eNum: '+str(edgeNum)  + "\n")
except Exception as e:
    print(e)