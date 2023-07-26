def P2k(k,t,m): #here k is half of the weight
    Ptot=0
    for j in range (0,k):
        #this is the combinatorial formula given in Chiriac and Jorza's paper "The trace of T2 takes no repeated values"
        Ptot= (-1)^j*binomial(2*k-2-j,j)*m^j*t^(2*k-2-2*j)+Ptot 
    return(Ptot)
def TrTm(k,m): #Using Zagier's appendix to Lang and weight 2k, only valid for k>2 or weight 4 and above (at weight 2 the trace is 0 anyway)
    Term1=-1/2*P2k(k,0,m)*pari(4*m).qfbhclassno() #The term where t=0
    Term2=0
    for t in range (1, floor(2*sqrt(m))+1):
        Term1=Term1-P2k(k,t,m)*pari(4*m-t^2).qfbhclassno() #qfbhclassno() is the Hurwitz class number
    divm=m.divisors()
    for i in range (0,len(divm)):
        d=divm[i]
        Term2=Term2-1/2*min(d,m/d)^(2*k-1)
    Tracem=Term1+Term2
    return(Tracem)

def a2CoefTm(k,m): #This is sum(a(i)a(j)), or the second coefficient of the Hecke Polynomial for Tm with weight 2k
    A=TrTm(k,m)^2
    B=TrTm(k,m^2)
    dk=floor(k/6)
    if k%6==1:
        dk=dk-1
    if k==1:
        dk=0
    Eigens=1/2*(A-B)-1/2*m^(2*k-1)*dk
    return(Eigens) #This works for any prime

for i in range (2,70):#used to check
    if a2CoefTm(i,2)== a2CoefTm(i,3):
        print(i)
