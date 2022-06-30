"""
This script will implement the Fox Algorithm on the wheel graph
W(7), and find all the covering spaces for 3-fold coverings. It
will then print out a data file of each cover with it's basic
topological information - the fundamental group.

Currently, it just returns the order of the group, which is a loss
of information since both Z and Z*Z have order "infinity". Should 
COUNT GENERATORS!

this file can be loaded into the interactive notebook with
sage.repl.load.load('FoxAlg-W7.sage',globals())

Can run straight from a termain with
sage FoxAlg-W7.sage 

NOTE: This is written for SAGEMATH! The .py is just to trick emacs
into giving us Python formatting ;-) Then you have to do something like
ln -s $PWD/FoxAlg-W7.py $PWD/FoxAlg-W7.sage
to make the sage invocation actually work... 

-----------------------------------------------------
Since this algorithm depends sensitively on which graph you actually
are using, a bunch of information needs to be specified up front. Namely,
the graph, the relators that give the correct knot group, and the order 
of the covers.

"""

##### Required data
# Wheel graph W6
fname = "W7-3.dat"
K=DiGraph({0:[1,2,3,4,5,6], 1:[2], 2: [3], 3:[4], 4:[5], 5:[6], 6:[1]})
Y=[[-12,-1,7],[-11,-6,12],[-10,-5,11],[-9,-4,10],[9,-8,-3],[-2,8,-7],[5,4,3,2,1,6]]
g=3 #order of the cover

##### Definitions first

# Rembember to cast input to this fellow as str!
def IsoImageTietze(im):
    im_list=im.split('*')
    #print(im_list)
    Tiet=[]
    for gen in im_list:
        
        #print(gen.partition("^")[0])
        # grab the number of the generator
        num=(gen.partition("^")[0]).partition("x")[2]
        #print(num)
        # grab the power of the generator
        pow=gen.partition("^")[2]
        if pow == "":
            #print("0")
            Tiet.append(int(num)+1)
        elif int(pow)>0:
            #print("+")
            for i in range(abs(int(pow))):
                Tiet.append(int(num)+1)
            
        elif int(pow)<0:
            #print("-")
            for i in range(abs(int(pow))):
                Tiet.append(-(int(num)+1))
    #print(Tiet)
    return(Tiet)

def word_hom(Y,PermList,g,E,i):
    YC=[0 for i in range(len(Y))]
    #print(YC)
    #print(Per)
    for j in range(len(Y)):
        # split into positive and negative cases. I think we need to remake each permutation
        # separately, since positive goes all before, and negative includes the permutation in question!
        Per=Permutations(3).identity()
        if Y[j]>0:
            #print("positive!")
            # this loop seems to work, but I think we have to skip the identity and start on the second 
            # term in the word.
            k=0
            while k<j:
                #print("------")
                #print(k)
                if Y[k]>0:
                    # Careful - permutations act FROM the right!
                    Per=Per*PerList[abs(Y[k])-1]
                else:
                    Per=Per*(PerList[abs(Y[k])-1].inverse())
                
                k+=1
                
            #print(Per.cycle_string())
            YC[j]=Y[j]+E*(Per(i+1)-1)
        else:
            #print("negative!")
            k=0
            while k<=j:
                #print("--------")
                #print(k)
                if Y[k]>0:
                    Per=Per*PerList[abs(Y[k])-1]
                else:
                    Per=Per*(PerList[abs(Y[k])-1].inverse())
                
                k+=1
                
            #print(Per.cycle_string())
            YC[j]=-1*(abs(Y[j])+E*(Per(i+1)-1))     
            
     
    return YC

# This function is intended as a "modular counter" to keep track of
# the covers - which are parameterized by selections of permutation
# elements for a particular graph.
# It takes a list of size n, corresponding to the number of free edges,
# a number p corresponding to the number of permutatoin elements (g!),
# a list i.e. L=[1,5,5] corresponding to the current count, and an instruction
# i for incrementing (and i should be 1 for normal operation)
# it returns the next L, like [2,1,1]...
#
# Careful with this function, I don't have many error checks
# around how far it counts...
def CountBase(n,p,L,i):
    if L[n-i] == p:
        L[n-i] = 1
        CountBase(n,p,L,i+1)
    elif L[n-i] < p:
        L[n-i]+=1
    return L

#
#########################################################
# setting the initial sizes of things.
E=K.size()
H=FreeGroup(E) #Initial free group

# knot group of complement, and the simplification isomorphism:
G=H/[H(Y[i]) for i in K.vertices()]
I=G.simplification_isomorphism()

# Now assign permutation elements to a maximal number of edges without
# contradicting group relations.
# First count how many "free" elements we have, which are just the elements
# mapped to themselves under the isomorphism.
P = Permutations(g)
free=0
for i in range(E):
    if str(I(G([i+1]))) == str(G([i+1])):
        #num = randrange(1,factorial(g))
        #PerList[i]=P[num]
        #fname += str(num)
        #outstring += str(PerList[i])
        free+=1
print(free)

# Now assign the free permutations, starting with all identity.
fout = open(fname,"w")
outstring = ""
L = [1 for i in range(free)]
L[free-1]=0 # count from zero!
#print(L)
# huge loop here over each cover!
NCOVER=0
per = Permutations(g)
P=per.list()
while NCOVER < factorial(g)**free:
    CountBase(free,factorial(g),L,1)
    outstring=str(L)+"\t"
    print("==========",L,"==============")
    PerList = [0 for i in range(E)]
    j=0
    #print(P[j],L)
    for i in range(E):
        if str(I(G([i+1]))) == str(G([i+1])):
            #num = randrange(1,factorial(g))
            #print("blah")
            PerList[i]=P[L[j]-1]
            #print(j,P[L[j]-1])
            #fname += str(num)
            #outstring += str(PerList[i])
            j+=1
        #print(i)
    print(PerList)
    #fname = fname + ".dat" # set the output file name.
    #outstring += "\n"

    print(j)
    # now loop through the edges that weren't set and set them,
    # using the isomorphism.
    for i in range(E):
        print(i)
        if str(I(G([i+1]))) != str(G([i+1])):
            print("new perm...",j)
            Per=per.identity()
            T=IsoImageTietze(str(I(G([i+1]))))
            for k in range(len(T)):
                if T[k]>0:
                    print(T,k,PerList[T[k]-1])
                    Per=Per*PerList[T[k]-1]
                else:
                    Per=Per*(PerList[abs(T[k])-1].inverse())
            PerList[i]=Per
    print(PerList)

    # add information to the output string.
    for i in PerList:
        outstring += str(i)
        #outstring += "\n"
    outstring+="\t"
    # Now create the relators in the cover using the word homomorphism:
    YC = []
    for i in range(g):
        for j in range(len(Y)):
            YC.append(word_hom(Y[j],PerList,g,E,i))

    # Create the underlying free group in the cover:
    H1=FreeGroup(E*g)

    # Calculate the fundamental group of the knot complement in the cover times
    # a free group of order g.
    G1=H1/[H1(YC[i]) for i in range(g*len(Y))]

    # now remove the free group by trivializing a word satisfying the Schreier
    # condition...
    SC=[[] for _ in range(g)] # list of words under the homomorphism.
    for i in range(g):
        S=[1 for k in range(i)] # change 2->1 after testing!
        SC[i]=word_hom(S,PerList,g,E,0)

    print(SC)

    # apply these relations to get the fundamental group of the knot complement.
    Rel=YC+SC
    G1=H1/[H1(Rel[i]) for i in range(g+g*len(Y))]
    print(G1)

    # Get the simplified version of this (not used at this time, but
    # possibly useful for analysis)
    G1.simplification_isomorphism()
    G1.order()
    outstring += str(G1.simplified()) +"\t" + str(G.order()) + "\t"

    # Now trivialize closed loops of the graph. First split the permutations
    # into Tietze words...
    VV=[]
    for i in range(len(PerList)):
        T=[]
        for s in PerList[i].cycle_string(singletons=True):
            if s=="(":
                V=[]
            elif s==",": # skip it!
                continue
            elif s==")": #end of permutation
                T.append(V)
            else:
                V.append(int(s))

        VV.append(T)
    print(VV)
    # Now send them to the cover....
    GG=[]
    for k in range(len(PerList)):
        GEN=[]
        for i in range(len(VV[k])): # each permutation indicates a new relation...
            VVV=[] #this is getting out of hand...
            for j in range(len(VV[k][i])):
                VVV=VVV+word_hom([k+1],PerList,g,E,VV[k][i][j]-1)
            GEN.append(VVV)
        GG=GG+GEN
    print(GG)

    # All the relators...
    Rel=YC+SC+GG;Rel

    # Finally, find the actual fundamental group of the cover!
    G1=H1/[H1(Rel[i]) for i in range(g+g*len(Y)+len(GG))]

    # Followed by basic analysis and recognition...
    # the structure description often fails, and it doesn't look like
    # any simple error handing can deal with it. So just skip!
    #G1.structure_description()
    G1.order()
    G1.simplification_isomorphism()
    print("Fundamental Group of the Cover:")
    G1.simplified()
    print(G1.simplified())
    outstring += str(G1.simplified()) + "\t" + str(G1.order()) + "\t" + str((G1.simplified()).ngens()) + "\t"
    if G1.simplified().relations() == ():
        outstring += "0" + "\n"
    else:
        outstring += str(G1.simplified().relations()) + "\n"
    #outstring += str(G1.simplified()) + "\t" + str(G1.order()) + "\t" + str(G1.structure_description()) + "\n"

    # write to the output file...
    #fout = open(fname,"w")
    fout.write(outstring)
    print("finished cover:", NCOVER+1)
    NCOVER+=1
fout.close()
