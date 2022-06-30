"""
This script will implement the Fox Algorithm for a set of wheel
graphs selected by Monte Carlo and determine the order
of the fundamental groups of the 3-fold coverings.

This variation will actually do a Monte Carlo over a single
graph at a time...(rand generator couldn't handle
min=max!)

The calls to structure_description that kills the code can
be safely removed here!

-----------------------------------------------------
Since this algorithm depends sensitively on which graph you actually
are using, some information needs to be specified up front. Namely,
the graph, the relators that give the correct knot group, and the order 
of the covers.

"""

##### MC-selected graph.
fname = "MC-W12-216.dat"
#NMin = 4
#NMax = 4
N = 12 # uncomment later N to fix....
TOT = 216 # number of covers to calculate
g=3 #order of the cover

##### Definitions first

# Rembember to cast input to this fellow as str!
def IsoImageTietze(im):
    im_list=im.split('*')
    #print(im_list)
    Tiet=[]
    for gen in im_list:
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

def word_hom(Y,PerList,g,E,i):
    YC=[0 for i in range(len(Y))]
    #print(YC)
    #print(Per)
    for j in range(len(Y)):
        # split into positive and negative cases. I think we need to remake each permutation
        # separately, since positive goes all before, and negative includes the permutation in question!
        Per=Permutations(g).identity()
        if Y[j]>0:
            #print("positive!")
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


# Abelianization function to find the homology group.
def Abelize(G):
    H = FreeGroup(G.ngens())
    GEN=[]
    for g1 in G.gens():
        for g2 in G.gens():
            GEN.append(g1*g2*g1^(-1)*g2^(-1))
    GG=H/GEN
    #return(GEN)
    return(GG)


#
#########################################################
fout = open(fname,"w")

# Going to start a giant loop here...
# Produce a dictionary of lists to graph our directed graph.
# and the relators for our particular presentation of the graph.
COUNT = 0
while COUNT < TOT:
    outstring=""
    #N = randrange(NMin,NMax)
    Kdict = {}
    keys = range(N)
    Y=[]
    for i in keys:
        Kdict[i]=[]
        Yrel=[]
        if i == 0:
            for j in range(1,N):
                Kdict[i].append(j)
                Yrel.append(j)
        elif i < N-1:
            Kdict[i].append(i+1)
            if i == 1:
                Yrel=[-i,-2*(N-1),N]
            else:
                Yrel=[-i,-(i+N-2),(i+N-1)]
        else:
            Kdict[i].append(1)
            Yrel=[-i,-(i+N-2),(i+N-1)]
        Y.append(Yrel)
    K=DiGraph(Kdict)
    print("Cover: ",COUNT,Kdict,Y)

    # setting the initial sizes of things.
    E=K.size()
    H=FreeGroup(E) #Initial free group
    
    # knot group of complement, and the simplification isomorphism:
    G=H/[H(Y[i]) for i in K.vertices()]
    I=G.simplification_isomorphism()
    print(I)

    # Now assign permutation elements randomly to the free edges.
    P=Permutations(g)
    PerList=[0 for i in range(E)]
    for i in range(E):
        if str(I(G([i+1]))) == str(G([i+1])):
            PerList[i]=P[randrange(1,factorial(g))]

    #print(PerList)
    # now loop through the edges that weren't set and set them,
    # using the isomorphism.
    for i in range(E):
        #print(i)
        if str(I(G([i+1]))) != str(G([i+1])):
            #print("new perm...",i)
            Per=P.identity()
            T=IsoImageTietze(str(I(G([i+1]))))
            #print(T)
            for k in range(len(T)):
                if T[k]>0:
                    Per=Per*PerList[T[k]-1]
                else:
                    Per=Per*(PerList[abs(T[k])-1].inverse())
                PerList[i]=Per
    print("List of Permutations:",PerList)

    # add information to the output string.
    outstring += str(Kdict) + "\t" + str(PerList) + "\t"

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
        S=[1 for k in range(i)] 
        SC[i]=word_hom(S,PerList,g,E,0)
    #print(SC)

    # apply these relations to get the fundamental group of the knot complement.
    Rel=YC+SC
    G1=H1/[H1(Rel[i]) for i in range(g+g*len(Y))]
    print(G1)

    # Get the simplified version of this (not used at this time, but
    # possibly useful for analysis)
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
    #print(VV)
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
    #print(GG)

    # All the relators...
    Rel=YC+SC+GG;Rel

    # Finally, find the actual fundamental group of the cover!
    G1=H1/[H1(Rel[i]) for i in range(g+g*len(Y)+len(GG))]

    # Followed by basic analysis and recognition...
    # the structure description often fails, and it doesn't look like
    # any simple error handing can deal with it. So just skip!
    #G1.structure_description()
    #G1.order()
    #G1.simplification_isomorphism()
    print("Fundamental Group of the Cover:", G1.simplified())
    outstring += str(G1.simplified()) + "\t" + str(G1.order()) + "\t" + str((G1.simplified()).ngens()) + "\t"
    # print the number of relations...
    if G1.simplified().relations() == ():
        outstring += "0" + "\t"
    else:
        outstring += str(G1.simplified().relations()) + "\t"
    # print the structure description, if possible:
    # Note that GAP doesn't deal well with free products,
    # So we are going to label those in this script.
    # we are also going to keep track of the number of exceptional
    # cases that GAP can't deal with and are also not clearly
    # Free Products.
    # Also going to wrap this in a timer tag, since sometimes
    # sage just hangs trying to solve the word problem...

    """
    alarm(5)
    try:
        grpstr = str(G1.structure_description())
        if grpstr == "Z":
            grpstr = "F1"
    except ValueError:
        if G1.simplified().relations() == ():
            grpstr = "F" + str(G1.simplified().ngens())
        else:
            grpstr = "--"
            pass
    except AlarmInterrupt:
        grpstr = "--"
        pass
    except:
        print("Some Unhandled exception!")
    cancel_alarm()
    """
    # instead, just make it blank!
    grpstr = "--"

    # write to the output file...
    outstring += grpstr + "\n"
    fout.write(outstring)
    COUNT+=1
    print("finished cover:", COUNT)

fout.close()
