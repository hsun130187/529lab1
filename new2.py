from __future__ import print_function
from __future__ import unicode_literals

import sys
from scipy import stats
import math


method=2

def file_to_data(file_to_read):
    """This function is uesed to read the data from the training set and make them into a dictionary.
        parameter: file_to_read. It is a file Object that stores the data.
         return type: dictionary
         returns: A dictionary with the data from the data file.
    """
   
    data = []
    for line in file_to_read.readlines():
        parts = line.split(',')
        if len(parts) == 2:
            idx, attrs = parts
            data.append({'id': idx,
                         'attrs': attrs.strip('\r\n')})
            
        elif len(parts) == 3:
            index, attrs, clas = parts
            data.append({'attrs': attrs,'class': clas.strip('\r\n')})
            
        else:
            print('not a correct format.')
            return
        
    return data

def get_most_frequent_class(data):
    """This function is used to identify and return the most frequent class in the data set.
        parameter:data, its type is dictionary
        return type: string,'ie' or 'EI' or 'N' 
    """
    clas = ''
    pValues = data_class_p_values(data)
    max_class = max(pValues)
    if pValues[0] is max_class:
        return 'EI'
    elif pValues[1] is max_class:
        return 'IE'
    else:
        return  'N'
    

'''def dna_count_class(dna_data):
    """counts the number of each class and returns them in this order - ei, ie
    and n.

    :type dna_data: dict
    :param dna_data: Set of parsed dna data.

    :rtype: tuple
    :returns: A tuple of the counts for (EI, IE, N)
    """
    ei_count = 0
    ie_count = 0
    n_count = 0
    for dna in dna_data:
        if dna['class'] == 'IE':
            ie_count += 1
        elif dna['class'] == 'EI':
            ei_count += 1
        else:
            n_count += 1
    return (ei_count, ie_count, n_count)'''

def criti_Value(data,method):
    if method==0:
        pValues = data_class_p_values(data)
        total = 0
        for p in pValues:
            if p != 0:
                total -= p * math.log(p, 2)
        return total
       
    if method==1:
        pValues = data_class_p_values(data)
        gini_value = 0
        if p_values:
            for p in pValues:
                gini_value += p*p
        return 1 - gini_value
        
    else:
        p_values = data_class_p_values(data)
        p_used=max(p_values)
        return 1-p_used
    
def gain_Value(data, values,attrIndex,method):
    size = len(data)
    sum_total = 0
 
    for value in values:
    
        subset = split_data(data, value, attrIndex)
        if subset:
            sum_total += (len(subset)/size) * criti_Value(subset,method)
    return criti_Value(data,method) - sum_total   


def split_data(data, value, attrIndex):
    """Gets a subset of the data where attr has the given value.

    :type dna_data: dict
    :param dna_data: Set of parsed dna data.

    :type values: list
    :param values: List of dna values

    :type attr: int
    :param attr: Attribute in the dna data

    :rtype: float
    :returns: Calculated gain based on the gini value.
    """
    subset = []
    for dna in data:
        if dna['attrs'][attrIndex] == value:
            # If the value of the attribute is what we are testing, add this dna to the subset
            subset.append(dna)
    return subset

def data_class_p_values(data):
    """Calculates probabilty of each of the 3 classes for a set of strands.

    :type dna_data: dict
    :param dna_data: Set of parsed dna data.

    :rtype: tuple
    :returns: A tuple of the probability for (EI, IE, N)
    """
    class_count=data_class_count(data)
    total = len(data)
    
    return (class_count[0]/total, class_count[1]/total, class_count[2]/total)
       
        
def issame_class(data):
    """Returns the class of the data if all the data shares the same class.

    :type dna_data: dict
    :param dna_data: Set of parsed dna data.

    :rtype: bool
    :returns: True if the dna data is all the same class, False otherwise
    """
    classes = [example['class'] for example in data]
   # if len(classes)==1:
        #return True
    #else:
        #return False
    return classes.count(classes[0]) == len(classes)


def allZero(labels):
    if labels:
        return False
    else:
        return True

def chooseBestFeatureToSplit(dataSet,labels):
    gain={}
    for ele in labels:
        gain[ele]=gain_Value(dataSet, values, ele,method)
    return max(gain, key=gain.get)    
    
def data_class_count(data):
    """counts the number of each class and returns them in this order - ei, ie
    and n.

    :type dna_data: dict
    :param dna_data: Set of parsed dna data.

    :rtype: tuple
    :returns: A tuple of the counts for (EI, IE, N)
    """
    ei_no = 0
    ie_no = 0
    n_no = 0
    for dna in data:
        if dna['class'] == 'IE':
            ie_no += 1
        elif dna['class'] == 'EI':
            ei_no += 1
        else:
            n_no += 1
    return (ei_no, ie_no, n_no)

def isBenificial(dataSet,bestFeat):
    p_values = data_class_p_values(dataSet)
    expected_count = []
    real_count = []
    featValues = [example['attrs'][bestFeat] for example in dataSet]#here need to be modified
    uniqueVals = set(featValues)
    for value in uniqueVals:
            childDataSet=[]
            childDataSet=split_data(dataSet,value,bestFeat)
            class_count = data_class_count(childDataSet)
            real_count.append(class_count[0])
            real_count.append(class_count[1])
            real_count.append(class_count[2])
            child_total = sum(class_count)
            expected_count.append(p_values[0]*child_total)
            expected_count.append(p_values[1]*child_total)
            expected_count.append(p_values[2]*child_total)
    #standards=stats.chi2.ppf(0.99,6)
    standards=0
    #print(standards)
    calValue=chi2Stac(real_count,expected_count)
    if calValue>standards:
        return True
    else:
        return False
   
            
def chi2Stac(realVal,expectedVal):
    result=0
    for i in range(len(expectedVal)):
        if expectedVal[i] != 0:
            result += ((realVal[i] - expectedVal[i])**2) / expectedVal[i]
    return result
            
def createTree(dataSet,labels):
    
    if allZero(labels)==True:
        #print('labels all zero')
        return get_most_frequent_class(dataSet)
    
    if issame_class(dataSet):
        #print('same classs')
        return get_most_frequent_class(dataSet)
    
    #print("time")
    bestFeat = chooseBestFeatureToSplit(dataSet,labels)
    #print(bestFeat)
    if not isBenificial(dataSet,bestFeat):
        #print("not benificial")
        #labels[bestFeat]=0
        return get_most_frequent_class(dataSet)
    else:
        #print("not benificial")
        myTree = {bestFeat:{}}
        labels.remove(bestFeat)
        #labels[bestFeat]=0
       # print(labels)
        #featValues = [example['attrs'][bestFeat] for example in dataSet]#here need to be modified
        #uniqueVals = set(featValues)
        for value in values:
            subset=split_data(dataSet,value,bestFeat)
            if not subset:
                myTree[bestFeat][value] = get_most_frequent_class(dataSet)
            else:
                subLabels=labels[:]
                myTree[bestFeat][value] = createTree(split_data(dataSet,value,bestFeat),subLabels)#not sure about this 
            
        return myTree
    
def classify(inputTree,featLabels,testVec):
    firstStr = next(iter(inputTree))
    secondDict = inputTree[firstStr]
    #featIndex = featLabels.index(firstStr)
    featIndex = firstStr
    #print('feat index is',featIndex)
    classLabel='nothing'
    for key in secondDict.keys():
        #print('testvex[2]is',testVec[featIndex])
        if testVec[featIndex] == key:
            if type(secondDict[key]).__name__=='dict':
                classLabel = classify(secondDict[key],featLabels,testVec)
                #print('is dic')
            else:
                classLabel = secondDict[key]
                #print("isvalue")
    return classLabel

def save_result_File(classification, result_file):
    """Saves the classification from the ID3 algorithm to a file.

    :type classification: list
    :param classification: The classification output from the ID3 algorithm for the testing data.

    :type classification_file: File Object
    :param classification_file: File to write the classification to.
    """
    print('id', 'class', file=result_file, sep=',')
    for item in classification:
        idx, cls = item.values()
        print(idx, cls, file=result_file, sep=',')
    return

def classify_data(inputTree,featLabels,data):
        """Find the classification for the given testing data.

        :type data: list
        :param data: Parsed testing data to classify.

        :rtype: list
        :returns: A classification for each data item in the testing data set.
        """
        classification = []
        for dna_data in data:
            idx, attrs = dna_data.values()
           # print(attrs)
            cls = classify(inputTree,featLabels,attrs)
            classification.append({'id': idx,
                                   'class': cls})
        return classification
    
     
file=open('training.csv','r')
data1=file_to_data(file)


file3=open('testing.csv','r')
data3=file_to_data(file3)

#useGini=False


values={'A','T','G','C'}
labels=[]
featLabels=[]
for i in range(60):
    labels.append(i)
    featLabels.append(i)

mytree=createTree(data1,labels)

resultfile=open('result2-ClassificationE-No.csv','w')

resultClass=classify_data(mytree,featLabels,data3)
save_result_File(resultClass, resultfile)

print('stop')

file.close()
file3.close()
resultfile.close()
