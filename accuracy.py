import sys

if __name__ == "__main__":
    f1 = open(sys.argv[1])
    gt = f1.readlines()
    f1.close()
    f1 = open(sys.argv[2])
    t = f1.readlines()
    f1.close()

    def intersection(lst1, lst2):
        lst3 = [value for value in lst1 if value in lst2]
        return lst3

    k = int(sys.argv[3])
    precision = 0.0
    recall = 0.0
    j = 0
    for x,y in zip(gt,t):
        z1 = [int(i) for i in x.split()]
        z2 = [int(i) for i in y.split()]
        l = len(intersection(z1,z2))
        j += 1
        precision += (l * 1.0) / k
        recall += (l * 1.0) / len(z1)

    print("Precision for", k, " predictions:", precision/j)
    print("Recall for", k, " predictions:", recall/j)
