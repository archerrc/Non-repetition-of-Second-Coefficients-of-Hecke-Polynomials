#note that in all these computations N is odd because N needs to be coprime to 2 and 4 for these formulas to be correct
def omega(n): #this is the number of prime divisors of n
    omegan=len(list(factor(n)))
    return(omegan)

def mu2(t,N): #this computes the value of mu for T2 and level N
    R.<x>=PolynomialRing(Integers(N))
    mutwopoly=x^2-t*x+2 #Finds all solutions in Z/NZ
    mutwo=mutwopoly.roots(multiplicities=False)
    mu2total=0
    for i in range (0,len(mutwo)):
        if gcd(mutwo[i],N)==1: #Finds which solutions are elements of (Z/NZ)*
            mu2total=mu2total+1
    return(mu2total)

def mu4(t,N):#this computes the values of mu relavant to TrT4 and level N
    R.<x>=PolynomialRing(Integers(N))
    mufourpoly=x^2-t*x+4 #Finds all solutions
    mufour=mufourpoly.roots(multiplicities=False)
    mu4total=0
    for i in range (0,len(mufour)):
        if gcd(mufour[i],N)==1: #Finds which solutions are elements of (Z/NZ)*
            mu4total=mu4total+1
    return(mu4total)
    
def A2T2(N,k): #here k is half of the weight, this is the A2 component given in the book "Traces of Hecke Operators" for T2
    if k>1:
        return(0)
    if k==1:
        Atwo=2
        if N%2==1:
            Atwo=3
    return(Atwo)

def A2T4(N,k): #here k is half of the weight, this is the A2 component given in the book "Traces of Hecke Operators" for T4
    if k>1:
        return(0)
    if k==1:
        AT4=4
        if N%2==1:
            AT4=7
    return(AT4)

def P2k(k,t,m): #here k is half of the weight
    Ptot=0
    for j in range (0,k):
        #this is the combinatorial formula given in Chiriac and Jorza's paper "The trace of T2 takes no repeated values"
        Ptot= (-1)^j*binomial(2*k-2-j,j)*m^j*t^(2*k-2-2*j)+Ptot 
    return(Ptot)

def DimFormula(N,k): #dimension formula for modular forms of level N and weight 2*k (rewritten to run a bit faster)
    Index=Gamma0(N).index()
    epsilon2=Gamma0(N).nu2()
    epsilon3=Gamma0(N).nu3()
    numCusps=Gamma0(N).ncusps()
    genus=1+Index/12-epsilon2/4-epsilon3/3-numCusps/2
    if k==1:
        return(genus)
    dim=(2*k-1)*(genus-1)+floor(k/2)*epsilon2+floor(2*k/3)*epsilon3+(k-1)*numCusps
    return(dim)
#the next two functions can be used to independently calculate the trace of the second or fourth Hecke operator for any (odd) level and (even integer) weight
def TrT2(k,N): #k is half of the weight and N is the level, note that mu(t,N)=mu(-t,N) for any t
    Term2=-1/4*P2k(k,2,2)*(2*mu2(2,N))-1/2*P2k(k,1,2)*(2*mu2(1,N))-1/2*P2k(k,0,2)*mu2(0,N)
    Term3=-2^(omega(N))
    Trace2=A2T2(N,k)+Term2+Term3
    return(Trace2)

def TrT4(k,N):
    Term2=((2*k-1)/12)*Gamma0(N).index()*4^(k-1)
    Term3=-1/2*P2k(k,3,4)*(2*mu4(3,N))-1/2*P2k(k,2,4)*(4/3*(2*mu4(2,N)))-P2k(k,1,4)*(2*mu4(1,N))-1/2*P2k(k,0,4)*3/2*mu4(0,N)
    Term4=-2^(omega(N))-2^(2*k-2)*Gamma0(N).ncusps() #The number of cusps in level N showing up here is interesting
    if N%9==0: #considers when gcd(t,N/t)=3
        if N%27!=0:
            Term4= Term4 -2^(omega(N))
        else:
            Term4= Term4 -2^(omega(N)+1)
    Trace4=A2T4(N,k)+Term2+Term3+Term4
    return(Trace4)

def EigenSum(k,N): #This is sum(a(i)a(j)) where the a(i) are eigenvalues of the trace of the second Hecke operator.
    A=TrT2(k,N)^2
    B=TrT4(k,N)
    dk=DimFormula(N,k)
    Eigens=1/2*(A-B)-2^(2*k-2)*dk
    return(Eigens)

def EigenSumCheck(N): #Check to level 4000 complete
    kmax=2 #these are the highest k-values that need to be checked for a given level, note that these values are one higher than the table, because the table is computing the difference between a2(k) and a2(k+1)
    if 0<N<178: #I know there's a better implementation for this, but this works
        kmax=718
    elif 177<N<316:
        kmax=400
    elif 315<N<714:
        kmax=300
    elif 713<N<2871:
        kmax=200
    elif 2870<N<5119:
        kmax=100
    elif 5118<N<11564:
        kmax=75
    elif 11563<N<18104:
        kmax=50
    elif 18103<N<32235:
        kmax=40
    elif 32234<N<46425:
        kmax=30
    elif 46424<N<50370:
        kmax=25
    elif 50369<N<54837:
        kmax=24
    elif 54836<N<59922:
        kmax=23
    elif 59921<N<65745:
        kmax=22
    elif 65744<N<72453:
        kmax=21
    elif 72452<N<80237:
        kmax=20
    elif 80236<N<89338:
        kmax=19
    elif 89337<N<100068:
        kmax=18
    elif 100067<N<112840:
        kmax=17
    elif 112839<N<128205:
        kmax=16
    elif 128204<N<146910:
        kmax=15
    elif 146909<N<169998:
        kmax=14
    elif 169997<N<198947:
        kmax=13
    elif 198946<N<235922:
        kmax=12
    elif 235921<N<284197:
        kmax=11
    elif 284196<N<348923:
        kmax=10
    elif 2<N<2:
        kmax=9
    elif 2<N<2:
        kmax=8
    elif 2<N<2:
        kmax=7
    elif 2<N<2:
        kmax=6
    elif 2<N<2:
        kmax=5
    elif 2<N<2:
        kmax=4
    elif 2<N<2:
        kmax=3
    a2List=[]
    Output=[]
    check=0
    if N>1185712:
        Output.append(('Check','unnecessary'))
        return(Output)
    GammaN=Gamma0(N)
    omegaN=omega(N)#necessary calculations that depend on the level and not the weight
    mu2t1=mu2(1,N)
    mu2t2=mu2(2,N)
    mu2t0=mu2(0,N)
    mu4t1=mu4(1,N)
    mu4t2=mu4(2,N)
    mu4t3=mu4(3,N)
    mu4t0=mu4(0,N)
    Index=GammaN.index()
    epsilon2=GammaN.nu2()
    epsilon3=GammaN.nu3()
    numCusps=GammaN.ncusps()
    genus=1+Index/12-epsilon2/4-epsilon3/3-numCusps/2
    TraceT2=1 #initializing variables
    TraceT4=1
    EigenSum=1
    dk=1
    for k in range (1,kmax+1): #here each term is computed directly so that computations dependant on the level are not needlessly repeated
        TraceT2=-1/4*P2k(k,2,2)*(2*mu2t2)-1/2*P2k(k,1,2)*(2*mu2t1)-1/2*P2k(k,0,2)*mu2t0-2^(omegaN)+A2T2(N,k) #t2 for this weight
        TraceT4=((2*k-1)/12)*Index*4^(k-1)-1/2*P2k(k,3,4)*(2*mu4t3)-1/2*P2k(k,2,4)*(4/3*(2*mu4t2))-P2k(k,1,4)*(2*mu4t1)-1/2*P2k(k,0,4)*3/2*mu4t0-2^(omegaN)-2^(2*k-2)*numCusps+A2T4(N,k)
        if N%9==0: #fixes a bug that was introduced when making the computation faster
            if N%27!=0:
                TraceT4= TraceT4 -2^(omega(N))
            else:
                TraceT4= TraceT4 -2^(omega(N)+1)
        if k==1:
            dk=genus
        else:
            dk=(2*k-1)*(genus-1)+floor(k/2)*epsilon2+floor(2*k/3)*epsilon3+(k-1)*numCusps
        EigenSum=1/2*(TraceT2^2-TraceT4)-2^(2*k-2)*dk
        Output.append(EigenSum)
    for j in range (0, len(Output)):
        check=Output[j]
        for l in range (0, len(Output)):
            if l==j:
                l=j
            elif check==Output[l]:
                a2List.append([check, 2*(j+1)])
    return(a2List)