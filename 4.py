from random import randint
import math

def EllipticCurveGeneration(z, k):
    # Генерация кривой с помощью параметризации
    while True:
        # характеристика поля
        p = 4*z**2 + 1
        if p in Primes():
            # число точек
            r = 4*z**2 - 2*z + 1
            if r in Primes():
                # след
                T = 1 + 2*z
                break
            # число точек
            else:
                r = 4*z**2 + 2*z + 1
                if r in Primes():
                    T = 1 - 2*z
                    break
        z = z + 1
    # дискриминант мнимого квадратичного порядка    
    D = 4*p - T**2
    print ("След T = {}".format(T))
    print ("Характеристика поля p = {}".format(p))
    print ("Число точек r = {}".format(r))
    print ("Дискриминант D = {}".format(D))
    print ("Переменная z = {}".format(z))
    
    # вычсление полинома Гильберта для D 
    H = hilbert_class_polynomial(-D)
    # генерация конечного поля характеристики p
    K = GF(p)
    # корни полинома Гильберта
    rg = H.roots(ring = K, multiplicities = False)
    if (rg == []):
        return 0, 0, 0
    j = rg[0]
    print ("j-инвариант j = {}".format(j))
    
    # расчет коэффициэнтов кривой в конечном поле K
    A = K(3*j / (1728 - j))
    B = K(2*j / (1728 - j))
    # генерация кривой над конечным полем K
    C = EllipticCurve(K, [A,B])
    # проверяем количество точек эллиптической кривой
    if gcd(C.order(), r) == 1:
        return 0, 0, 0
    print ("{}".format(C))
    
    # генерация расширения конечного поля
    EK.<t> = GF(p**k)
    # генерация кривой над расширенным конечным полем
    EC = EllipticCurve(EK, [A, B])
    if gcd(EC.order(), r**2) == 1:
        return 0, 0, 0
    print ("{}".format(EC))
    return C, EC, r    

# выбираем точку P порядка r на кривой E(F_p)
def GetP(C, EC, r):
    P = C.random_point()
    while (True):
        if (P.order() == r):
            break
        P = C.random_point()
    P = EC(P[0], P[1], P[2])
    return P

# выбираем точку PK порядка r на кривой E(F_q), q = p^k
def GetPK(EC, r):
    PK = EC.random_point()
    while (True):
        PK *= (EC.order() // r ** 2)
        if PK.order() == r:
            break
        PK = EC.random_point()
    return PK

# вычисляем группу кручения
def TorsionSubgroup(P, PK, EZ, r):
    tor = []
    PK_1 = []
    # добавляем бесконечно-удаленную точку
    PK_1.append(EZ) 
    # элементы получаем по формуле 𝑎𝑃+𝑏𝑃𝐾,0≤𝑎,𝑏≤𝑟−1
    for j in range (1, r):
            PK_j = j * PK
            PK_1.append(PK_j)
    tor.append(PK_1)
    for i in range (1, r):
        P_sum = []
        P_i = i * P
        P_sum.append(P_i)
        for j in range (1, r):
            PK_j = j * PK
            P_sum.append(P_i + PK_j)
        tor.append(P_sum)
    # подсчет количества элементов в группе кручения
    l = 0
    for i in range (0, r):
        l +=  len(tor[i])
    print("Количество элементов в группе кручения: {}".format(l))
    return tor

# проверка свойств билинейности и знакопеременности для спаривание Вейля
def WeilPairing(tor, r):
    print("\n_______Спаривание Вейля_______")
    S1 = tor[0][1]
    S2 = tor[0][2]
    T1 = tor[1][0]
    T2 = tor[2][0]
    
    # билинейность e(S1 + S2, T) == e(S1, T) e(S2, T); e(S, T1 + T2) == e(S, T1) e(S, T2); 
    S = S1 + S2
    eSST = S.weil_pairing(T1, r)
    eST = S1.weil_pairing(T1, r)
    eTS = S2.weil_pairing(T1, r)
    ee = eST * eTS
    
    T = T1 + T2
    eSST1 = T.weil_pairing(S1, r)
    eST1 = T1.weil_pairing(S1, r)
    eTS1 = T2.weil_pairing(S1, r)
    ee1 = eST1 * eTS1
        
    if(eSST == ee and eSST1 == ee1):
        print("\nБилинейность по Вейлю: e(S1 + S2, T) == e(S1, T) e(S2, T); e(S, T1 + T2) == e(S, T1) e(S, T2);")
    print("e(S1 + S2, T) = {}".format(eSST),"\ne(S1, T) = {}".format(eST),"\ne(S2, T) = {}".format(eTS),"\ne(S1, T) e(S2, T) = {}".format(ee))
    print("\ne(S, T1 + T2) = {}".format(eSST1),"\ne(S, T1) = {}".format(eST1),"\ne(S, T2) = {}".format(eTS1),"\ne(S, T1) e(S, T2) = {}".format(ee1))
   
    # знакопеременность e(S, S) = 1
    sm = 0
    lenT = 0
    for i in range (0, r):
        for j in range (0, r):
            S = tor[i][j]
            eSS = S.weil_pairing(S, r)
            lenT += 1
            if(eSS == 1):
                sm += 1
    if(sm == lenT):
        print("\nЗнакопеременность по Вейлю: e(S, S) = 1 для всех S ∈ E[r]")
    #print("sm = ", sm)

# проверка свойств билинейности и знакопеременности для спаривание Тейта
def TatePairing(tor, r, k):
    print("\n\n_______Спаривание Тейта_______")
    S1 = tor[0][1]
    S2 = tor[0][2]
    T1 = tor[1][0]
    T2 = tor[2][0]
    
    # билинейность e(S1 + S2, T) == e(S1, T) e(S2, T); e(S, T1 + T2) == e(S, T1) e(S, T2); 
    S = S1 + S2
    eSST = S.tate_pairing(T1, r, k)
    eST = S1.tate_pairing(T1, r, k)
    eTS = S2.tate_pairing(T1, r, k)
    ee = eST * eTS
    
    T = T1 + T2
    eSST1 = T.tate_pairing(S1, r, k)
    eST1 = T1.tate_pairing(S1, r, k)
    eTS1 = T2.tate_pairing(S1, r, k)
    ee1 = eST1 * eTS1
        
    if(eSST == ee and eSST1 == ee1):
        print("\nБилинейность по Тейту: e(S1 + S2, T) == e(S1, T) e(S2, T); e(S, T1 + T2) == e(S, T1) e(S, T2);")
    print("e(S1 + S2, T) = {}".format(eSST),"\ne(S1, T) = {}".format(eST),"\ne(S2, T) = {}".format(eTS),"\ne(S1, T) e(S2, T) = {}".format(ee))
    print("\ne(S, T1 + T2) = {}".format(eSST1),"\ne(S, T1) = {}".format(eST1),"\ne(S, T2) = {}".format(eTS1),"\ne(S, T1) e(S, T2) = {}".format(ee1))
   
    # знакопеременность e(S, S) = 1
    sm = 0
    lenT = 0
    for i in range (0, r):
        for j in range (0, r):
            S = tor[i][j]
            eSS = S.tate_pairing(S, r, k)
            lenT += 1
            if(eSS == 1):
                sm += 1
    if(sm == lenT):
        print("\nЗнакопеременность по Тейту: e(S, S) = 1 для всех S ∈ E[r]\n")
    else:
        print("\nЗнакопеременность по Тейту не выполняется: sm = ", sm)
        print("")

def simplehash(mes, i):
    s = 0
    for y in mes:
        s += ord(y) + 1# орд возвращает код символа
    return s % i  

def IntoPoint(EC, mess):
    p = int(EK.base_field().characteristic()) # р- модуль, А, В
    A = int(EK.a4())
    B = int(EK.a6())
    K = GF(p)
    for i in range(2, p):
        # хэш-функция от сообщения
        h = simplehash(mess, i)
        # символ Лежандра ((h**3 + A*h + B)/p)
        t = kronecker_symbol(h**3 + A*h + B, p)
    # y^2 = b mod(p)? если да то (t/p)=1
        if t == 1:
            var('y')
            UR=h**3+A*h+B==y**2
            print("UR={}".format(UR))
            y_r = UR.roots(ring=K,multiplicities=False)
            y_res =y_r[randint(0,p)%len(y_r)]
    # выбираем произвольно один из корней
            try:
                P = EK([h, y_res])
                return P
            except:
                continue
    P0 = EC([0,1,0])
    return P0 

# генерируем ключи
def KeyGeneration(tor, r):
    # секретный ключ
    x = randint(1, r - 1)
    # случайная  точка из группы кручения
    i, j = randint(0, r - 1), randint(0, r - 1)
    Q = tor[i][j] 
    # открытый ключ
    V = x * Q 
    return x, Q, V

# вложение сообщения в x-координату точки
def MessToPoint(EC, mess):
    p = int(EC.base_field().characteristic())
    A = int(EC.a4())
    B = int(EC.a6())
    K = GF(p)
    #print ('p = {}'.format(p))
    #print ('A = {}'.format(A))
    #print ('B = {}'.format(B))
    for i in range(2, p):
        # хэш-функция от сообщения
        h = simplehash(mess, i)
        # символ Лежандра (h**3 + A*h + B/p)
        t = kronecker_symbol(h**3 + A * h + B, p)
        if t == 1:
            # решаем уравнение
            var('y')
            UR = h**3 + A * h + B == y**2
            y_r = UR.roots(ring = K, multiplicities = False)
            # выбираем произвольно один из корней
            res = y_r[randint(0,p) % len(y_r)]
            try:
                P = EC([h, res])
                # print ('P = {}'.format(P))
                return P
            except:
                continue
    P0 = EC([0,1,0])
    return P0 

# формирование подписи
def Signing(x, EC, mess):
    M = MessToPoint(EC, mess)
    Y = x * M  
    return Y[0]

def Check (EC, mess, Q, V, S, r):
    print ('Сообщение: {}'.format(mess))
    M = MessToPoint(EC, mess)
    print ('M = {}'.format(M))
    print ('S = {}'.format(S))

    e1 = Q.weil_pairing(S, r)
    e2 = V.weil_pairing(M, r)
    print("e1 = {}".format(e1),"\ne2 = {}".format(e2))
    
    if (e1 == e2) or (e1 == e2^(-1)):
        print ('Корректная подпись\n')
    else:
        print ('Некорректная подпись\n')

# проверка подписи
def Verification(mess, EC, tor, r, Xs, Q, V):
    p = EC.base_field().characteristic()
    A = EC.a4()
    B = EC.a6()
    KT = GF(p)
    
    Ys = KT(Xs**3 + A*Xs + B)
    Ys = KT(sqrt(Ys))
    try:       
        S = EC([Xs, Ys])
    except:
        Ys = KT(-Ys)
        S = EC([Xs, Ys])      
    Check (EC, mess, Q, V, S, r)
    Check (EC, "error", Q, V, S, r)
    
# короткая подпись на билинейных отображениях        
def ShortSignatureWP(EC, tor, r, mess):
    x, Q, V = KeyGeneration(tor, r)
    S = Signing(x, EC, mess)
    Verification(mess, EC, tor, r, S, Q, V)
          
def main():    
    # переменная для параметризации
    z = 2
    # степень вложения
    k = 6
    while True:
        C, EC, r = EllipticCurveGeneration(z, k)
        if r == 0:
            print ('Error')
        break
    P = GetP(C, EC, r)
    print("P = {}".format(P))
    PK = GetPK(EC, r)
    print("PK = {}".format(PK))
    tor = TorsionSubgroup(P, PK, EC(0), r)
    WeilPairing(tor, r)
    TatePairing(tor, r, k)
    mess = "Hello"
    ShortSignatureWP(EC, tor, r, mess)
    
if __name__ == "__main__":
    main()    