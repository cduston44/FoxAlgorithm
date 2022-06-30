"""
This script will try to replicate the multiple exception call
errors I am getting without all the excess baggage associated with
the Fox Algorithm.
"""
g=3 # Covers, but this is just a multiplicative factor here
E=5 # edges in the graph
V=4 # Sides in the graph (relations)
W=5 # Max word size for the relations

H=FreeGroup(E*g)

N=1000
ExCases = 0 # To keep track of exceptional cases....
for i in range(N):
    
    # Randomly create some relations
    Y = []
    for j in range(V*g):
        YY = []
        for k in range(randrange(1,W)):
            num = randrange(1,g*E)
            sign = randrange(1,2) # positive or negative
            if sign == 1:
                YY.append(num)
            elif sign == 2:
                YY.append(-num)
        #print(YY)
        Y.append(YY)
    print(i,Y)

    # Now calculate the structure of the group
    Y1=[H(Y[i]) for i in range(g*V)]
    G=H/Y1
    #G=H/[H(Y[i]) for i in range(g*V)]
    print(Y1,G)

    alarm(5)
    try:
        grpstr = str(G.structure_description())
    except ValueError:
        if G.simplified().relations() == ():
            grpstr = "F" + str(G.simplified().ngens())
        else:
            grpstr = "--"
            ExCases += 1
            pass
    except AlarmInterrupt:
        grpstr = "--"
        ExCases += 1
        pass
    except:
        print("Some Unhandled exception!")
        pass
    cancel_alarm()

    print("Group Type is:",grpstr)
print("Exceptional cases:",ExCases)
