import sys
def dbn_output(name,dbn,rna,dbn_file):
    f = open(dbn_file,"w")
    f.write("#Name: "+name + "\n")
    f.write("#Length: "+str(len(rna)) +"\n")
    f.write("#PageNumber: NA\n")
    f.write(rna + "\n")
    f.write(dbn)
    f.close()

if __name__ == "__main__":
    dbn_file = sys.argv[1]
    dbn = "((..))"
    rna = "AAGATT"
    dbn_output("test",dbn,rna,dbn_file)