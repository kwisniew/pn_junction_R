######komentarz############
#zak�adamy: stan r�wnowagi, brak przy�o�onego napi�cia, generacja = rekombinacja
#modu� ma 4um^2, a warstw zubozona niech ma 0.2um
#3 rownania: poissona, ci�g�o�ci dla elektron�w, ci�g�o�ci dla dziur

####ROZMIAR URZ�DZENIA w [cm]####
dlugosc_x=0.3*10^-3
dlugosc_y=0.3*10^-3
width = 0.3*10^-3

#########Domieszkowanie i poczatkowe koncentracje####################

#dziel� na 200 punkt�w. B�d� mia� 3*200x3*200 = 9*40000 punkt�w.
#Z czego po iksie z 200 punkt�w 10 przypada na warstw� zubo�on� (mam nadziej�, �e wystarczy).  
#Ciekawe ile to si� b�dzie rozwi�zywa�o?
#LICZBA PUNK�W MUSI BY� PARZYSTA!!! (tak p�niej indeksuj�)
l_pkt <- 50
h=width/(l_pkt-1)

domieszkowanie_D <- 10^15
#koncentracja elektorn�w po stronie n
n_D <- domieszkowanie_D
domieszkowanie_A <- 10^15
#koncentracja dziur po stronie p
p_A <- domieszkowanie_A

skalowanie_npC <- max(c(domieszkowanie_D, domieszkowanie_A))

#domieszkowanie zdefiniowane dla ca�ego materia�u 2D
domieszkowanie <- matrix(rep(domieszkowanie_A/skalowanie_npC,(l_pkt))^2,ncol = l_pkt, nrow = l_pkt)
#domieszkowanie na typ N b�dzie dawa�o dodatni wk�ad do �adunku, a po stornie typu P ujemny
#{zasada zachowania �adunku}
domieszkowanie[,seq(1,(l_pkt/2),1)] <- domieszkowanie_D/skalowanie_npC
domieszkowanie[,seq((l_pkt/2+1),l_pkt,1)] <- -domieszkowanie_A/skalowanie_npC



#liczenie pocz�tkowych koncentracji: dziur po stonie n i el. po stonir p
NC <- 2.2*10^18
NV <- 1.8*10^18
Eg <- 1.15
k  <- 1.38*10^(-23)
T  <- 300
q  <- 1.6*10^(-19)
ni <- sqrt(NC*NV)*exp(-Eg*1.6*10^(-19)/(2*k*T))
#p_D czytaj: koncentracja dziur po stonie N (donorowej)
p_D<- ni^2/n_D/skalowanie_npC
n_A<- ni^2/p_A/skalowanie_npC

n_D=domieszkowanie_D/skalowanie_npC
p_A=domieszkowanie_A/skalowanie_npC

##########warunki poczatkowe################
#Ustalam warunki poczatkowe: potencja� zero wsz�dzie opr�cz warstwu zubo�onej, przeskalowane warto�ci g�sto�ci no�nik�w
#czyli: w warstwie zubo�onej n,p=0; n=Nd po stornie n, p=Na po stornie p; n po stornie p to ni^2/p i analogicznie 
# p=ni^2/p

#b�dziemy mieli 3N niewiadomych. Oznacza to, �e musimy stworzy� macierz 3Nx3N aby w pe�ni rozwi�za� liniowy uk�ad r�wna�
#zlinearyzujemy go tak jak Gummel i Scharfettel. Stworzymy te� siatk� zawieraj�c� 200 punkt�w po "x" i 200p. po "y"
#w ka�dym punkcie b�d� 3 interesuj�ce nas zmienne n,p,fi mo�emy zapisa� jako: n_i, p_i, fi_i. B�dziemy wi�c mieli
#200x200 punkt�w - oznaczmy to jako N, do ka�dego  niewiadome, czyli razem 3N niewiadomych. Tyle samo musimy mie� uk�ad
#r�wna�. Dla ka�dego punktu b�d� wi�c 3 r�wnania. Komponuj�c te r�wnania dla ka�dego punktu otrzymamy macie� 3Nx3N.
#ustalamy na pocz�tku warto�ci pocz�tkowe naszych 3N niewiadomych, aby algorytm m�g� wystartowa�.
#Aby by�o wygodnie szuka� najbli�szych s�siad�w b�d� to tak naprawd� 3 po��czone macierze:
#macierz warunk�w poczatkowych dla fi, dla n i dla p

warunki_poczatkowe <- matrix(rep(0,(3*l_pkt)*l_pkt), nrow = 3*l_pkt, ncol = l_pkt)
#wierszami: 1:200 fi, 201:400 n, 401:600 p
warunki_poczatkowe[          1:l_pkt    ,         ]=0
warunki_poczatkowe[(  l_pkt+1):(2*l_pkt),          1:(l_pkt/2)]=n
warunki_poczatkowe[(  l_pkt+1):(2*l_pkt),(l_pkt/2+1):(l_pkt  )]=n_A
warunki_poczatkowe[(2*l_pkt+1):(3*l_pkt),          1:(l_pkt/2)]=p_D
warunki_poczatkowe[(2*l_pkt+1):(3*l_pkt),(l_pkt/2+1):l_pkt]=p

#warstwa zubo�ona w przybli�eniu pe�nego zubo�enia
#staramy si� aby w~=0.4x10^-5cm = 04um, do powiedzmy 0.8um
#- czyli jak d�ugo�� jest 4/8um to 10/20% to b�dzie warstwa zubo�ona
#b�d� to punkty �rodkowe. je�li np. po "x" punkt�w jest 200,
#to warstwa zubo�ona zajmie 20 �rodkowych, 10 po stornie N i 10 po stornie P
#gdy� domieszkowanie na dzie� 08.11.2016 przyjmuje za takie same po obu stronach
eps=13.6*8.9*10^-14
fi_bi <- k*T*log(domieszkowanie_A*domieszkowanie_D/ni^2)/q
w = sqrt(2*eps*fi_bi*(1/domieszkowanie_D + 1/domieszkowanie_A)/q)
#d�ugo�� warstwy zubo�onej po stronie N i P.  
xP <- w/(1+domieszkowanie_A/domieszkowanie_D)
xN <- -(w-xP)
#na podstawie Chenming Hu, s.131
pot_war_zub <- function(x)
{ 
  ifelse(x < xN, fi_bi, ifelse((x<=0 & x>=xN),
                          fi_bi-q*domieszkowanie_D*(x-xN)^2/(2*eps),ifelse((x>0 & x<=xP),
                                                                      q*domieszkowanie_A*(xP-x)^2/(2*eps), 0)))
}
h=width/(l_pkt-1)
x <- c(-((l_pkt/2):1)+(1/2),1:(l_pkt/2)-(1/2))*h
Vt  <- k*T/q
skalowanie_fi <- Vt

#NIE SKALOWA� X PRZED U�YCIEM PONI�SZEJ FUNKCJI pot_war_zub !!!!!!!!!!!!!!
fi<-pot_war_zub(x)/Vt

#skalowanie warto�ci x.
h  <-  h/width
xP <- xP/width
xN <- xN/width
x  <-  x/width
w  <-  w/width

tmp <- matrix(rep(fi, l_pkt), nrow = l_pkt, ncol = l_pkt, byrow = TRUE)
#ustalamy nasze warto�ci pocz�tkowe dla warstwy zubo�onej -potencja�
warunki_poczatkowe[1:l_pkt,]=tmp
#te same warto�ci pocz�tkowe b�d� nam potrzebna do rowi�zania r�wniania poissona - powiedzmy, �e to te r�wniania 
#s� opisywane przez wiersze 1:200, ale te same b�d� nam potrzeba do r�wna� ci�g�o�ci elektron�w i dziur
#(patrz poni�ej)

#liniowy wzrost koncentracji w warstwie zubozonej - w razie potrzeby
{
#for(i in (x[x<=xP & x>=xN]-xN))
#{
#  delta_n = n_D-n_A
#  delta_p =-p_D+p_A
#  warunki_poczatkowe[(  l_pkt+1):(2*l_pkt),which((x-xN)==i)]=n_D-(i/w)*delta_n
#  warunki_poczatkowe[(2*l_pkt+1):(3*l_pkt),which((x-xN)==i)]=    (i/w)*delta_p
#}
}  

warunki_poczatkowe[(  l_pkt+1):(2*l_pkt),x<=xP & x>=xN]=0
warunki_poczatkowe[(2*l_pkt+1):(3*l_pkt),x<=xP & x>=xN]=0

#############Najblizszy sasiad##############
najblizszy_sasiad_wewnatrz <- function(i)
{
  c(i-1,i,i+1,i-l_pkt, i+l_pkt)
}
najblizszy_sasiad_kontakty <- function(i)
{
  i
}
najblizszy_sasiad_boki     <- function(i)
{
  c(i-1,i,i+1)
}

#################bernoulli i jakobian#########################
#Teraz trzeba przej�� do okre�lenia jakobianu - zasadniczej cz�ci metody Newtona.
#Wiemy, �e dzyskretyzowane r�wnianie zawiera w sobie funkcj� bernouliego x/(e^x-1)
#gdzie x to dwie zmienne: fi1-fi2, gdzie fi1/2 to potencja�y w s�siaduj�cych punktach
#powstaje jednak problem co si� dzieje w punkcie 0, w Numerical Models for Differential Problems
#dodaj�, �e B(0)=1 aby funkcja by�a ci�g�a. aby zrobi� ci�g�� pochod� k�adk� dB(0)/dx = 1.5
#zdefiniujmy j�:

#wersja analityczna bernoulliego
{
#bernoulli <- function(x,y)
#{
#  ifelse((x-y)!=0, (x-y)/(exp(x-y)-1),1)
#  
#}
}

#wersja numeryczna por. Selberherr str 169
bernoulli <- function(x,y)
{
  z1=-37.42995
  z2=-8.72095e-06
  z3=6.580134e-08
  z4=37.42995
  z5=745
  
  z<-x-y
  if(z <= z1){z <- -z}
  else if(z>z1 & z<z2){z <- z/(exp(z)-1)}
  else if(z>=z2 & z<=z3){z <- 1-z/2}
  else if(z>z3 & z<z4){z <- z*exp(-z)/(1-exp(-z))}
  else if(z>=z4 & z<z5){z <- z*exp(-z)}
  else if(z>=z5) {z<- 0}
  else print("Nie potrafi�!!")
  z
}

diff_bernoulli <- function(x,y,d)
{
  if(d=="x")
  {
    #print("pochodna po x")
    ifelse((x-y)!=0, -exp(y)*(exp(x)*(x-y-1)+exp(y))/(exp(x)-exp(y))^2, -1.5)
  }
  else
  {
    if(d=="y")
    {
      #print("pochodna po y")
      ifelse((x-y)!=0, exp(y)*(exp(x)*(x-y-1)+exp(y))/(exp(x)-exp(y))^2, -1.5)
    }
    else
    {
      print("zjeba�e�!!!")
    }
  }
}
#dla punkt�w wewn�trznych mamy mamdy r�wnania zgodne z dyskretyzacj� Scharfettera Gummela
#... a w metodzie Newtona liczy si� jakobian, czyli liczymy pochodne cz�stkowe

#ruchliwo�ci (w cm^2V^-1s^-1) i inne potrzebne sta�e i przeskalowania (s.143 Selberherr)
u_n <- 100
u_p <- 25
Dn  <- u_n*Vt
Dp  <- u_p*Vt
skalowanie_D <- max(Dn, Dp)
Dn <- Dn/skalowanie_D
Dp <- Dp/skalowanie_D
lambda_kwadrat <- Vt*eps/(width^2*q*skalowanie_npC)

#przepisanie niewiadomych do 1 wiersza (200^2 niewiadomych fi, 200^2 n i 200^2 p)
#transponujemy bo R zaczyna przepisywa� po kolumnach
#!!! Zmienna wyniki b�dzie ITEROWANA
wyniki <- t(warunki_poczatkowe)
dim(wyniki) <- c(3*l_pkt^2,1)
#sprawdzenie, czy si� dobrze uda�o niewiadome ustawi�
#plot(wyniki[(0*200+1):(1*200)], type = 'l')
#plot(wyniki[(200*200+1):(201*200)], type = 'l')
#plot(wyniki[(400*200+1):(401*200)], type = 'l')


jakobian <- matrix(rep(0,(3*l_pkt^2)^2),nrow = 3*l_pkt^2, ncol = 3*l_pkt^2)
#wspolczynniki dla rowniania poissona - pierwsze l_pkt^2 wierszy macierzy jakobianu.

#indeksy istotne dla i-tego wiersza w rusznych miejscach struktury


#okre�lenie kt�re indeksy odpowiadaj� za kolejne miejsca w strukturze (kontakty,boki{sztuczne},wn�trze)
lewa_strona <- seq(1    ,(3*l_pkt^2),l_pkt)
prawa_strona<- seq(l_pkt,(3*l_pkt^2),l_pkt)
indeksy_kontaktow <- c(lewa_strona, prawa_strona)
indeksy_kontaktow <- sort(indeksy_kontaktow)
indeksy_bokow     <- c(seq(2,l_pkt-1,1), seq(l_pkt*(l_pkt-1)+2,l_pkt^2-1))
indeksy_bokow     <- c(indeksy_bokow, indeksy_bokow+l_pkt^2, indeksy_bokow+2*l_pkt^2)
indeksy_brzegow   <- sort(c(indeksy_bokow, indeksy_kontaktow))
indeksy_wszystkie <- 1:(3*l_pkt^2)
indeksy_wewnatrz  <- indeksy_wszystkie[-indeksy_brzegow]


#rozprawiamy si� z poissonem 
jakobian_poisson_wewnatrz <- function(i)
{
  c( lambda_kwadrat/h^2, -4*lambda_kwadrat/h^2, lambda_kwadrat/h^2, lambda_kwadrat/h^2, lambda_kwadrat/h^2, -1, 1)
}
jakobian_poisson_boki     <- function(i)
{
  c( lambda_kwadrat/h^2, -2*lambda_kwadrat/h^2, lambda_kwadrat/h^2, -1, 1)
}
jakobian_poisson_kontakty <- function(i)
{
  c(0,-1,1)
}

#rozprawiamy si� z r�wnaniami ci�g�o�ci
#!!! FUNCKJE ZACZYNAJ�CE SI� OD "jakobian_" U�YWAJ� GLOBALNEGO WEKTORA WYNIKI!!! 
jakobian_continuity_wewnatrz <- function(i, jaki_nosnik)
{
  if(jaki_nosnik=="n")
  {
    #print("elektrony")
    fi<- wyniki[najblizszy_sasiad_wewnatrz(i)]
    q <- wyniki[najblizszy_sasiad_wewnatrz(i)+l_pkt^2]
    D <- Dn
  }
  else
  {
    if(jaki_nosnik=="p")
    {
      #print("dziury")
      #zamieniami znak w r�wnaniu przy wszystkich fi --> patrz strona 159 Selbererr
      fi<- -wyniki[najblizszy_sasiad_wewnatrz(i)]
      q <-  wyniki[najblizszy_sasiad_wewnatrz(i)+2*l_pkt^2]
      D <-  Dp 
    }
    else
    {
      print("ZJEBA�E�!!!")
    }
  }
  
  
  jn1 <- (D*diff_bernoulli(fi[1],fi[2],"x")*q[1]-D*diff_bernoulli(fi[2],fi[1],"y")*q[2])/h^2
  jn2 <-
  (
     D*diff_bernoulli(fi[1],fi[2],"y")*q[1]+D*diff_bernoulli(fi[3],fi[2],"y")*q[3]+D*diff_bernoulli(fi[4],fi[2],"y")*q[4]+D*diff_bernoulli(fi[5],fi[2],"y")*q[5]+
   (-D*diff_bernoulli(fi[2],fi[1],"x")     -D*diff_bernoulli(fi[2],fi[3],"x")     -D*diff_bernoulli(fi[2],fi[4],"x")     -D*diff_bernoulli(fi[2],fi[5],"x"))*q[2]
  )/h^2
  jn3 <- (D*diff_bernoulli(fi[3],fi[2],"x")*q[3]-D*diff_bernoulli(fi[2],fi[3],"y")*q[2])/h^2
  jn4 <- (D*diff_bernoulli(fi[4],fi[2],"x")*q[4]-D*diff_bernoulli(fi[2],fi[4],"y")*q[2])/h^2
  jn5 <- (D*diff_bernoulli(fi[5],fi[2],"x")*q[5]-D*diff_bernoulli(fi[2],fi[5],"y")*q[2])/h^2
  
  n1  <- D*bernoulli(fi[1],fi[2])/h^2
  n2  <- (-D*bernoulli(fi[2],fi[1])-D*bernoulli(fi[2],fi[3])-D*bernoulli(fi[2],fi[4])-D*bernoulli(fi[2], fi[5]))/h^2
  n3  <- D*bernoulli(fi[3],fi[2])/h^2
  n4  <- D*bernoulli(fi[4],fi[2])/h^2
  n5  <- D*bernoulli(fi[5],fi[2])/h^2
  
  c(jn1,jn2,jn3,jn4,jn5,n1,n2,n3,n4,n5)

}
jakobian_continuity_wewnatrz_elektrony <- function(i)
{

  fi<-  wyniki[najblizszy_sasiad_wewnatrz(i)]
  q <-  wyniki[najblizszy_sasiad_wewnatrz(i)+l_pkt^2]
  D <-  Dn 

  jn1 <- (D*diff_bernoulli(fi[1],fi[2],"x")*q[1]-D*diff_bernoulli(fi[2],fi[1],"y")*q[2])/h^2
  jn2 <-
    (
      D*diff_bernoulli(fi[1],fi[2],"y")*q[1]+D*diff_bernoulli(fi[3],fi[2],"y")*q[3]+D*diff_bernoulli(fi[4],fi[2],"y")*q[4]+D*diff_bernoulli(fi[5],fi[2],"y")*q[5]+
        (-D*diff_bernoulli(fi[2],fi[1],"x")     -D*diff_bernoulli(fi[2],fi[3],"x")     -D*diff_bernoulli(fi[2],fi[4],"x")     -D*diff_bernoulli(fi[2],fi[5],"x"))*q[2]
    )/h^2
  jn3 <- (D*diff_bernoulli(fi[3],fi[2],"x")*q[3]-D*diff_bernoulli(fi[2],fi[3],"y")*q[2])/h^2
  jn4 <- (D*diff_bernoulli(fi[4],fi[2],"x")*q[4]-D*diff_bernoulli(fi[2],fi[4],"y")*q[2])/h^2
  jn5 <- (D*diff_bernoulli(fi[5],fi[2],"x")*q[5]-D*diff_bernoulli(fi[2],fi[5],"y")*q[2])/h^2
  
  n1  <- D*bernoulli(fi[1],fi[2])/h^2
  n2  <- (-D*bernoulli(fi[2],fi[1])-D*bernoulli(fi[2],fi[3])-D*bernoulli(fi[2],fi[4])-D*bernoulli(fi[2], fi[5]))/h^2
  n3  <- D*bernoulli(fi[3],fi[2])/h^2
  n4  <- D*bernoulli(fi[4],fi[2])/h^2
  n5  <- D*bernoulli(fi[5],fi[2])/h^2
  
  c(jn1,jn2,jn3,jn4,jn5,n1,n2,n3,n4,n5)
  
}
jakobian_continuity_wewnatrz_dziury <- function(i)
{
  
  fi<-  wyniki[najblizszy_sasiad_wewnatrz(i)]
  q <-  wyniki[najblizszy_sasiad_wewnatrz(i)+2*l_pkt^2]
  D <-  Dp 
  
  jn1 <- (D*diff_bernoulli(fi[2],fi[1],"y")*q[1]-D*diff_bernoulli(fi[1],fi[2],"x")*q[2])/h^2
  jn2 <-
    (
      D*diff_bernoulli(fi[2],fi[1],"x")*q[1]+D*diff_bernoulli(fi[2],fi[3],"x")*q[3]+D*diff_bernoulli(fi[2],fi[4],"x")*q[4]+D*diff_bernoulli(fi[2],fi[5],"x")*q[5]+
    (-D*diff_bernoulli(fi[1],fi[2],"y")     -D*diff_bernoulli(fi[3],fi[2],"y")     -D*diff_bernoulli(fi[4],fi[2],"y")     -D*diff_bernoulli(fi[5],fi[2],"y"))*q[2]
    )/h^2
  jn3 <- (D*diff_bernoulli(fi[2],fi[3],"y")*q[3]-D*diff_bernoulli(fi[3],fi[2],"x")*q[2])/h^2
  jn4 <- (D*diff_bernoulli(fi[2],fi[4],"y")*q[4]-D*diff_bernoulli(fi[4],fi[2],"x")*q[2])/h^2
  jn5 <- (D*diff_bernoulli(fi[2],fi[5],"y")*q[5]-D*diff_bernoulli(fi[5],fi[2],"x")*q[2])/h^2
  
  n1  <- D*bernoulli(fi[2],fi[1])/h^2
  n2  <- (-D*bernoulli(fi[1],fi[2])-D*bernoulli(fi[3],fi[2])-D*bernoulli(fi[4],fi[2])-D*bernoulli(fi[5], fi[2]))/h^2
  n3  <- D*bernoulli(fi[2],fi[3])/h^2
  n4  <- D*bernoulli(fi[2],fi[4])/h^2
  n5  <- D*bernoulli(fi[2],fi[5])/h^2
  
  c(jn1,jn2,jn3,jn4,jn5,n1,n2,n3,n4,n5)
}

jakobian_continuity_boki     <- function(i, jaki_nosnik)
{
  if(jaki_nosnik=="n")
  {
    #print("elektrony")
    fi<- wyniki[najblizszy_sasiad_boki(i)]
    q <- wyniki[najblizszy_sasiad_boki(i)+l_pkt^2]
    D <- Dn
  }
  else
  {
    if(jaki_nosnik=="p")
    {
      #print("dziury")
      fi<- -wyniki[najblizszy_sasiad_boki(i)]
      q <- wyniki[najblizszy_sasiad_boki(i)+2*l_pkt^2]
      D <- Dp 
    }
    else
    {
      print("ZJEBA�E�!!!")
    }
  }
  
  jq1 <- (D*diff_bernoulli(fi[1],fi[2],"x")*q[1]-D*diff_bernoulli(fi[2],fi[1],"y")*q[2])/h^2
  jq2 <-
    (
       D*diff_bernoulli(fi[1],fi[2],"y")*q[1]+D*diff_bernoulli(fi[3],fi[2],"y")*q[3]+
     (-D*diff_bernoulli(fi[2],fi[1],"x")     -D*diff_bernoulli(fi[2],fi[3],"x"))*q[2]
    )/h^2
  jq3 <- (D*diff_bernoulli(fi[3],fi[2],"x")*q[3]-D*diff_bernoulli(fi[2],fi[3],"y")*q[2])/h^2

  
  q1  <- D*bernoulli(fi[1],fi[2])/h^2
  q2  <- (-D*bernoulli(fi[2],fi[1])-D*bernoulli(fi[2],fi[3]))/h^2
  q3  <- D*bernoulli(fi[3],fi[2])/h^2
  
  c(jq1,jq2,jq3,q1,q2,q3)
}
jakobian_continuity_boki_elektrony     <- function(i)
{

  fi<- wyniki[najblizszy_sasiad_boki(i)]
  q <- wyniki[najblizszy_sasiad_boki(i)+l_pkt^2]
  D <- Dn

  
  jq1 <- (D*diff_bernoulli(fi[1],fi[2],"x")*q[1]-D*diff_bernoulli(fi[2],fi[1],"y")*q[2])/h^2
  jq2 <-
    (
      D*diff_bernoulli(fi[1],fi[2],"y")*q[1]+D*diff_bernoulli(fi[3],fi[2],"y")*q[3]+
        (-D*diff_bernoulli(fi[2],fi[1],"x")     -D*diff_bernoulli(fi[2],fi[3],"x"))*q[2]
    )/h^2
  jq3 <- (D*diff_bernoulli(fi[3],fi[2],"x")*q[3]-D*diff_bernoulli(fi[2],fi[3],"y")*q[2])/h^2
  
  
  q1  <- D*bernoulli(fi[1],fi[2])/h^2
  q2  <- (-D*bernoulli(fi[2],fi[1])-D*bernoulli(fi[2],fi[3]))/h^2
  q3  <- D*bernoulli(fi[3],fi[2])/h^2
  
  c(jq1,jq2,jq3,q1,q2,q3)
}
jakobian_continuity_boki_dziury     <- function(i)
{
  
  fi<- wyniki[najblizszy_sasiad_boki(i)]
  q <- wyniki[najblizszy_sasiad_boki(i)+2*l_pkt^2]
  D <- Dp
  
  
  jq1 <- (D*diff_bernoulli(fi[2],fi[1],"y")*q[1]-D*diff_bernoulli(fi[1],fi[2],"x")*q[2])/h^2
  jq2 <-
    (
      D*diff_bernoulli(fi[2],fi[1],"x")*q[1]+D*diff_bernoulli(fi[2],fi[3],"x")*q[3]+
        (-D*diff_bernoulli(fi[1],fi[2],"y")     -D*diff_bernoulli(fi[3],fi[2],"y"))*q[2]
    )/h^2
  jq3 <- (D*diff_bernoulli(fi[2],fi[3],"y")*q[3]-D*diff_bernoulli(fi[3],fi[2],"x")*q[2])/h^2
  
  
  q1  <- D*bernoulli(fi[2],fi[1])/h^2
  q2  <- (-D*bernoulli(fi[1],fi[2])-D*bernoulli(fi[3],fi[2]))/h^2
  q3  <- D*bernoulli(fi[2],fi[3])/h^2
  
  c(jq1,jq2,jq3,q1,q2,q3)
}

jakobian_continuity_kontakty <- function(i, jaki_nosnik)
{
  if(jaki_nosnik=="n")
  {
    #print("elektrony")
  }
  else
  {
    if(jaki_nosnik=="p")
    {
      #print("dziury")
    }
    else
    {
      print("ZJEBA�E�!!!")
    }
  }
  c(0,0)
}
jakobian_continuity_kontakty_elektrony <- function(i, jaki_nosnik)
{
  c(0,0)
}
jakobian_continuity_kontakty_dziury <- function(i, jaki_nosnik)
{
  c(0,0)
}

##################rozwi�zania dla warunkow poczatkowych###################

#Obliczamy wyniki r�wna� poissona dla WATUNKOW POCZATKOWYCH
#!!! W zmiennej wartsc_funkcji b�d� przechowywane warto�ci funkcji do ka�dego nast�pnego kroku w metodzie Newtona
#b�dzie to wi�c ZMIENNA ITEROWANA
wartosci_funkcji <- rep(0,3*l_pkt^2)

#POISSON WYNIKI
poisson_wewnatrz <- function(i)
{
  fi <- wyniki[najblizszy_sasiad_wewnatrz(i)]
  ni <- wyniki[i+l_pkt^2]
  pi <- wyniki[i+2*l_pkt^2]
  ((fi[1]+fi[3]+fi[4]+fi[5]-4*fi[2])*lambda_kwadrat)/h^2 - ni + pi + domieszkowanie[floor(((i-1)/l_pkt)+1),((i-1)%%l_pkt+1)]
}
poisson_boki     <- function(i)
{
  fi <- wyniki[najblizszy_sasiad_boki(i)]
  ni <- wyniki[i+l_pkt^2]
  pi <- wyniki[i+2*l_pkt^2]
  ((fi[1]+fi[3]-2*fi[2])*lambda_kwadrat)/h^2 - ni + pi + domieszkowanie[floor(((i-1)/l_pkt)+1),((i-1)%%l_pkt+1)]
}
poisson_kontakty <- function(i)
{
  ni <- wyniki[i+l_pkt^2]
  pi <- wyniki[i+2*l_pkt^2]
  -ni + pi + domieszkowanie[floor(((i-1)/l_pkt)+1),((i-1)%%l_pkt+1)]
}

#obliczamy wynik r�wna� ci�g�o�ci dla warunk�w pocz�tkowych
continuity_wewnatrz <- function(i,jaki_nosnik)
{
  #patrz�c na przypisanie g�sto�ci �adunk�w q to i nie mo�e przekracza� l_pkt^2
  if(i>l_pkt^2){print("Za du�y indeks!! Dozowolony maksymalnie l_pkt^2")}
  
  if(jaki_nosnik=="n")
  {
    #print("elektrony")
    fi<- wyniki[najblizszy_sasiad_wewnatrz(i)]
    q <- wyniki[najblizszy_sasiad_wewnatrz(i)+l_pkt^2]
    D <- Dn
  }
  else
  {
    if(jaki_nosnik=="p")
    {
      #print("dziury")
      #zamieniami znak w r�wnaniu przy wszystkich fi --> patrz strona 159 Selbererr
      fi<- -wyniki[najblizszy_sasiad_wewnatrz(i)]
      q <-  wyniki[najblizszy_sasiad_wewnatrz(i)+2*l_pkt^2]
      D <-  Dp 
    }
    else
    {
      print("ZJEBA�E�!!!")
    }
  }
  (D*q[1]*bernoulli(fi[1],fi[2])-D*q[2]*bernoulli(fi[2],fi[1])+
   D*q[3]*bernoulli(fi[3],fi[2])-D*q[2]*bernoulli(fi[2],fi[3])+
   D*q[4]*bernoulli(fi[4],fi[2])-D*q[2]*bernoulli(fi[2],fi[4])+
   D*q[5]*bernoulli(fi[5],fi[2])-D*q[2]*bernoulli(fi[2],fi[5]))/h^2;
  
}
continuity_wewnatrz_elektrony <- function(i,jaki_nosnik)
{
  if(i>l_pkt^2){print("Za du�y indeks!! Dozowolony maksymalnie l_pkt^2")}

  fi<- wyniki[najblizszy_sasiad_wewnatrz(i)]
  q <- wyniki[najblizszy_sasiad_wewnatrz(i)+l_pkt^2]
  D <- Dn

  (D*q[1]*bernoulli(fi[1],fi[2])-D*q[2]*bernoulli(fi[2],fi[1])+
    D*q[3]*bernoulli(fi[3],fi[2])-D*q[2]*bernoulli(fi[2],fi[3])+
    D*q[4]*bernoulli(fi[4],fi[2])-D*q[2]*bernoulli(fi[2],fi[4])+
    D*q[5]*bernoulli(fi[5],fi[2])-D*q[2]*bernoulli(fi[2],fi[5]))/h^2;
  
}
continuity_wewnatrz_dziury <- function(i,jaki_nosnik)
{
  #patrz�c na przypisanie g�sto�ci �adunk�w q to i nie mo�e przekracza� l_pkt^2
  if(i>l_pkt^2){print("Za du�y indeks!! Dozowolony maksymalnie l_pkt^2")}
  
  fi<-  wyniki[najblizszy_sasiad_wewnatrz(i)]
  q <-  wyniki[najblizszy_sasiad_wewnatrz(i)+2*l_pkt^2]
  D <-  Dp 

   (D*q[1]*bernoulli(fi[2],fi[1])-D*q[2]*bernoulli(fi[1],fi[2])+
    D*q[3]*bernoulli(fi[2],fi[3])-D*q[2]*bernoulli(fi[3],fi[2])+
    D*q[4]*bernoulli(fi[2],fi[4])-D*q[2]*bernoulli(fi[4],fi[2])+
    D*q[5]*bernoulli(fi[2],fi[5])-D*q[2]*bernoulli(fi[5],fi[2]))/h^2;
  
}
continuity_boki <- function(i,jaki_nosnik)
{
  #patrz�c na przypisanie g�sto�ci �adunk�w "q" to "i" nie mo�e przekracza� l_pkt^2
  if(i>l_pkt^2){print("Za du�y indeks!! Dozowolony maksymalnie l_pkt^2")}
  
  if(jaki_nosnik=="n")
  {
    #print("elektrony")
    fi<- wyniki[najblizszy_sasiad_boki(i)]
    q <- wyniki[najblizszy_sasiad_boki(i)+l_pkt^2]
    D <- Dn
  }
  else
  {
    if(jaki_nosnik=="p")
    {
      #print("dziury")
      #zamieniami znak w r�wnaniu przy wszystkich fi --> patrz strona 159 Selbererr
      fi<- -wyniki[najblizszy_sasiad_boki(i)]
      q <-  wyniki[najblizszy_sasiad_boki(i)+2*l_pkt^2]
      D <-  Dp 
    }
    else
    {
      print("ZJEBA�E�!!!")
    }
  }
  (D*q[1]*bernoulli(fi[1],fi[2])-D*q[2]*bernoulli(fi[2],fi[1])+
   D*q[3]*bernoulli(fi[3],fi[2])-D*q[2]*bernoulli(fi[2],fi[3]))/h^2;
  
}
continuity_boki_elektrony <- function(i,jaki_nosnik)
{
  #patrz�c na przypisanie g�sto�ci �adunk�w "q" to "i" nie mo�e przekracza� l_pkt^2
  if(i>l_pkt^2){print("Za du�y indeks!! Dozowolony maksymalnie l_pkt^2")}
  
  fi<- wyniki[najblizszy_sasiad_boki(i)]
  q <- wyniki[najblizszy_sasiad_boki(i)+l_pkt^2]
  D <- Dn

  (D*q[1]*bernoulli(fi[1],fi[2])-D*q[2]*bernoulli(fi[2],fi[1])+
    D*q[3]*bernoulli(fi[3],fi[2])-D*q[2]*bernoulli(fi[2],fi[3]))/h^2;
  
}
continuity_boki_dziury <- function(i,jaki_nosnik)
{
  #patrz�c na przypisanie g�sto�ci �adunk�w "q" to "i" nie mo�e przekracza� l_pkt^2
  if(i>l_pkt^2){print("Za du�y indeks!! Dozowolony maksymalnie l_pkt^2")}
  
  fi<- wyniki[najblizszy_sasiad_boki(i)]
  q <- wyniki[najblizszy_sasiad_boki(i)+2*l_pkt^2]
  D <- Dp
  
  (D*q[1]*bernoulli(fi[2],fi[1])-D*q[2]*bernoulli(fi[1],fi[2])+
   D*q[3]*bernoulli(fi[2],fi[3])-D*q[2]*bernoulli(fi[3],fi[2]))/h^2;
  
}
continuity_kontakty <- function(i,jaki_nosnik)
{
  0
}
continuity_kontakty_elektrony <- function(i,jaki_nosnik)
{
  0
}
continuity_kontakty_dziury <- function(i,jaki_nosnik)
{
  0
}

###############funkcje do obliczania warto�ci r�wna� (poissona i ci�g�o�ci) i jakobianu##############
#!!!!!!!korzystaj� ze zmiennej globalnej "wyniki"
oblicz_wartosc_funkcji  <- function()
{
  #POTENCJA�
  for(i in indeksy_kontaktow[indeksy_kontaktow<=l_pkt^2])
  {
    wartosci_funkcji[i]=poisson_kontakty(i)
  }
  for(i in indeksy_bokow[    indeksy_bokow    <=l_pkt^2])
  {
    wartosci_funkcji[i]=poisson_boki(i)
  }
  for(i in indeksy_wewnatrz[ indeksy_wewnatrz <=l_pkt^2])
  {
    wartosci_funkcji[i]=poisson_wewnatrz(i)
  }
  #ELEKTRONY
  for(i in indeksy_kontaktow[indeksy_kontaktow>l_pkt^2 & indeksy_kontaktow<=2*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_kontakty(i-l_pkt^2, "n")
  }
  for(i in indeksy_bokow[    indeksy_bokow    >l_pkt^2 & indeksy_bokow    <=2*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_boki(i-l_pkt^2, "n")
  }
  for(i in indeksy_wewnatrz[ indeksy_wewnatrz  >l_pkt^2 & indeksy_wewnatrz <=2*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_wewnatrz(i-l_pkt^2, "n")
  }
  #DZIURY
  for(i in indeksy_kontaktow[indeksy_kontaktow>2*l_pkt^2 & indeksy_kontaktow<=3*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_kontakty(i-2*l_pkt^2, "p")
  }
  for(i in indeksy_bokow[    indeksy_bokow    >2*l_pkt^2 & indeksy_bokow    <=3*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_boki(i-2*l_pkt^2, "p")
  }
  for(i in indeksy_wewnatrz[ indeksy_wewnatrz  >2*l_pkt^2 & indeksy_wewnatrz <=3*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_wewnatrz(i-2*l_pkt^2, "p")
  }
  wartosci_funkcji
}
oblicz_wartosc_funkcji2 <- function()
{
  #POTENCJA�
  for(i in indeksy_kontaktow[indeksy_kontaktow<=l_pkt^2])
  {
    wartosci_funkcji[i]=poisson_kontakty(i)
  }
  for(i in indeksy_bokow[    indeksy_bokow    <=l_pkt^2])
  {
    wartosci_funkcji[i]=poisson_boki(i)
  }
  for(i in indeksy_wewnatrz[ indeksy_wewnatrz <=l_pkt^2])
  {
    wartosci_funkcji[i]=poisson_wewnatrz(i)
  }
  #ELEKTRONY
  for(i in indeksy_kontaktow[indeksy_kontaktow>l_pkt^2 & indeksy_kontaktow<=2*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_kontakty_elektrony(i-l_pkt^2)
  }
  for(i in indeksy_bokow[    indeksy_bokow    >l_pkt^2 & indeksy_bokow    <=2*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_boki_elektrony(i-l_pkt^2)
  }
  for(i in indeksy_wewnatrz[ indeksy_wewnatrz  >l_pkt^2 & indeksy_wewnatrz <=2*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_wewnatrz_elektrony(i-l_pkt^2)
  }
  #DZIURY
  for(i in indeksy_kontaktow[indeksy_kontaktow>2*l_pkt^2 & indeksy_kontaktow<=3*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_kontakty_dziury(i-2*l_pkt^2)
  }
  for(i in indeksy_bokow[    indeksy_bokow    >2*l_pkt^2 & indeksy_bokow    <=3*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_boki_dziury(i-2*l_pkt^2)
  }
  for(i in indeksy_wewnatrz[ indeksy_wewnatrz  >2*l_pkt^2 & indeksy_wewnatrz <=3*l_pkt^2])
  {
    wartosci_funkcji[i]=continuity_wewnatrz_dziury(i-2*l_pkt^2)
  }
  wartosci_funkcji
}
oblicz_jakobian         <- function()
{
  for(i in indeksy_kontaktow[indeksy_kontaktow<=l_pkt^2])
  {
    jakobian[i,c(najblizszy_sasiad_kontakty(i),i+l_pkt^2,i+2*l_pkt^2)]=jakobian_poisson_kontakty(i)
  }
  for(i in indeksy_bokow[    indeksy_bokow    <=l_pkt^2])
  {
    jakobian[i,c(najblizszy_sasiad_boki(i),i+l_pkt^2,i+2*l_pkt^2)]=jakobian_poisson_boki(i)
  }
  for(i in indeksy_wewnatrz[ indeksy_wewnatrz <=l_pkt^2])
  {
    jakobian[i,c(najblizszy_sasiad_wewnatrz(i),i+l_pkt^2,i+2*l_pkt^2)]=jakobian_poisson_wewnatrz(i)
  }
  
  #ELEKTRONY
  for(i in indeksy_kontaktow[(indeksy_kontaktow>l_pkt^2 & indeksy_kontaktow<=2*l_pkt^2)])
  {
    jakobian[i,c(najblizszy_sasiad_kontakty(i-l_pkt^2),najblizszy_sasiad_kontakty(i))]=jakobian_continuity_kontakty(i-l_pkt^2, "n")
  }
  for(i in indeksy_bokow[    (indeksy_bokow    >l_pkt^2 & indeksy_bokow    <=2*l_pkt^2)])
  {
    jakobian[i,c(najblizszy_sasiad_boki(    i-l_pkt^2),najblizszy_sasiad_boki(i    ))]=jakobian_continuity_boki(    i-l_pkt^2, "n")
  }
  for(i in indeksy_wewnatrz[ (indeksy_wewnatrz >l_pkt^2 & indeksy_wewnatrz <=2*l_pkt^2)])
  {
    jakobian[i,c(najblizszy_sasiad_wewnatrz(i-l_pkt^2),najblizszy_sasiad_wewnatrz(i))]=jakobian_continuity_wewnatrz(i-l_pkt^2, "n")
  }
  
  #DZIURY
  for(i in indeksy_kontaktow[(indeksy_kontaktow>(2*l_pkt^2) & indeksy_kontaktow<=(3*l_pkt^2))])
  {
    jakobian[i,c(najblizszy_sasiad_kontakty(i-2*l_pkt^2),najblizszy_sasiad_kontakty(i))]=jakobian_continuity_kontakty(i-2*l_pkt^2, "p")
  }
  for(i in indeksy_bokow[     (indeksy_bokow>(   2*l_pkt^2) & indeksy_bokow    <=(3*l_pkt^2))])
  {
    jakobian[i,c(najblizszy_sasiad_boki(    i-2*l_pkt^2),najblizszy_sasiad_boki(i    ))]=jakobian_continuity_boki(    i-2*l_pkt^2, "p")
  }
  for(i in indeksy_wewnatrz[  (indeksy_wewnatrz>(2*l_pkt^2) & indeksy_wewnatrz <=(3*l_pkt^2))])
  {
    jakobian[i,c(najblizszy_sasiad_wewnatrz(i-2*l_pkt^2),najblizszy_sasiad_wewnatrz(i))]=jakobian_continuity_wewnatrz(i-2*l_pkt^2, "p")
  }
  jakobian
}
oblicz_jakobian2        <- function()
{
  for(i in indeksy_kontaktow[indeksy_kontaktow<=l_pkt^2])
  {
    jakobian[i,c(najblizszy_sasiad_kontakty(i),i+l_pkt^2,i+2*l_pkt^2)]=jakobian_poisson_kontakty(i)
  }
  for(i in indeksy_bokow[    indeksy_bokow    <=l_pkt^2])
  {
    jakobian[i,c(najblizszy_sasiad_boki(i),i+l_pkt^2,i+2*l_pkt^2)]=jakobian_poisson_boki(i)
  }
  for(i in indeksy_wewnatrz[ indeksy_wewnatrz <=l_pkt^2])
  {
    jakobian[i,c(najblizszy_sasiad_wewnatrz(i),i+l_pkt^2,i+2*l_pkt^2)]=jakobian_poisson_wewnatrz(i)
  }
  
  #ELEKTRONY
  for(i in indeksy_kontaktow[(indeksy_kontaktow>l_pkt^2 & indeksy_kontaktow<=2*l_pkt^2)])
  {
    jakobian[i,c(najblizszy_sasiad_kontakty(i-l_pkt^2),najblizszy_sasiad_kontakty(i))]=jakobian_continuity_kontakty_elektrony(i-l_pkt^2)
  }
  for(i in indeksy_bokow[    (indeksy_bokow    >l_pkt^2 & indeksy_bokow    <=2*l_pkt^2)])
  {
    jakobian[i,c(najblizszy_sasiad_boki(    i-l_pkt^2),najblizszy_sasiad_boki(i    ))]=jakobian_continuity_boki_elektrony(    i-l_pkt^2)
  }
  for(i in indeksy_wewnatrz[ (indeksy_wewnatrz >l_pkt^2 & indeksy_wewnatrz <=2*l_pkt^2)])
  {
    jakobian[i,c(najblizszy_sasiad_wewnatrz(i-l_pkt^2),najblizszy_sasiad_wewnatrz(i))]=jakobian_continuity_wewnatrz_elektrony(i-l_pkt^2)
  }
  
  #DZIURY
  for(i in indeksy_kontaktow[(indeksy_kontaktow>(2*l_pkt^2) & indeksy_kontaktow<=(3*l_pkt^2))])
  {
    jakobian[i,c(najblizszy_sasiad_kontakty(i-2*l_pkt^2),najblizszy_sasiad_kontakty(i))]=jakobian_continuity_kontakty_dziury(i-2*l_pkt^2)
  }
  for(i in indeksy_bokow[     (indeksy_bokow>(   2*l_pkt^2) & indeksy_bokow    <=(3*l_pkt^2))])
  {
    jakobian[i,c(najblizszy_sasiad_boki(    i-2*l_pkt^2),najblizszy_sasiad_boki(i    ))]=jakobian_continuity_boki_dziury(    i-2*l_pkt^2)
  }
  for(i in indeksy_wewnatrz[  (indeksy_wewnatrz>(2*l_pkt^2) & indeksy_wewnatrz <=(3*l_pkt^2))])
  {
    jakobian[i,c(najblizszy_sasiad_wewnatrz(i-2*l_pkt^2),najblizszy_sasiad_wewnatrz(i))]=jakobian_continuity_wewnatrz_dziury(i-2*l_pkt^2)
  }
  jakobian
}

#### Pierwsza PR�BA#########
wyniki <- t(warunki_poczatkowe)
dim(wyniki) <- c(3*l_pkt^2,1)
#sprawdzam co sie stanie jak ruchliwosci beda takie same
Dp=Dn
jakobian         <- oblicz_jakobian2()
wartosci_funkcji <- oblicz_wartosc_funkcji2()
error <- sum(sapply(wartosci_funkcji[-indeksy_kontaktow],function(x){x^2}))/length(wartosci_funkcji[-indeksy_kontaktow])
tk=0.1
tlumienie <- 1/tk
wyniki[-indeksy_kontaktow] <- solve(tlumienie*jakobian[-indeksy_kontaktow,-indeksy_kontaktow],-wartosci_funkcji[-indeksy_kontaktow], tol=1e-20)+ wyniki[-indeksy_kontaktow]
#poni�szy spos�b jest ponad 4 razy wolniejszy:
{
#wyniki[-indeksy_kontaktow] <- wyniki[-indeksy_kontaktow]-solve(jakobian[-indeksy_kontaktow,-indeksy_kontaktow])%*%wartosci_funkcji[-indeksy_kontaktow]
}
#wyniki_po_pierwszej_probie <- wyniki
#wyniki <- wyniki_po_pierwszej_probie
#########G��WNA P�TLA##############

#Napiszmy p�tl� w kt�rej bedzie wykonywana metoda Newtona. Na razie 10 iteracji
ptm <- proc.time()
for(i in 1:6)
{
 #Potrzebujemy 3 elemtent�w: wynik�w dla n-1, dla tych wynik�w obliczone: warto�ci funkcji n-1 oraz warto�ci jakobianu 
 #wyniki bierzemy z warunk�w pocz�tkowych dla n=1, dla reszty te wyniki s� obliczne z r�wnania z linijki 508
 #potrzebjemy warto�ci funkcji:
  wartosci_funkcji <- oblicz_wartosc_funkcji2()
  print(sum(sapply(wartosci_funkcji[-indeksy_kontaktow],function(x){x^2}))/length(wartosci_funkcji[-indeksy_kontaktow]))
  error <- append(sum(sapply(wartosci_funkcji[-indeksy_kontaktow],function(x){x^2}))/length(wartosci_funkcji[-indeksy_kontaktow]),error)
  #liczymy jakobian
  jakobian         <- oblicz_jakobian2()
  wyniki[-indeksy_kontaktow] <- (
                                 solve(tlumienie*jakobian[-indeksy_kontaktow,-indeksy_kontaktow],-wartosci_funkcji[-indeksy_kontaktow], tol=1e-30)
                                 + wyniki[-indeksy_kontaktow]
                                )
  #wyniki[-indeksy_kontaktow] <- wyniki[-indeksy_kontaktow]-solve(jakobian[-indeksy_kontaktow,-indeksy_kontaktow], tol=1e-20)%*%wartosci_funkcji[-indeksy_kontaktow]
}
proc.time()-ptm
write.csv(wyniki, "wyniki_po_36_iteracjach_tlumienie_10")
####sprawdzam wykresami####
library(rgl)

plot3d(rep((1:l_pkt),l_pkt),matrix(rep(1:l_pkt,l_pkt),ncol=l_pkt,byrow=T),wyniki[1:l_pkt^2], type = 'l')
plot3d(rep((1:l_pkt),l_pkt),matrix(rep(1:l_pkt,l_pkt),ncol=l_pkt,byrow=T),wyniki[(l_pkt^2+1):(2*l_pkt^2)], type = 'l')
plot3d(rep((1:l_pkt),l_pkt),matrix(rep(1:l_pkt,l_pkt),ncol=l_pkt,byrow=T),wyniki[(2*l_pkt^2+1):(3*l_pkt^2)])

powierzchnia_wyniki <- function(x,y)
{
  wyniki[x+(y-1)*l_pkt]
}
powierzchnia_wartosci_funkcji <- function(x,y)
{
  wartosci_funkcji[x+(y-1)*l_pkt]
}
wykres_powierzchni <- function(jakie_rownanie,powierzchnia)
{
  x <- 1:l_pkt
  y <- x
  if(jakie_rownanie=='fi')
  {
    z <- outer(x, y, powierzchnia) 
    persp3d(x,y,z, alpha=0.5)
  }
  else
  {
    if(jakie_rownanie=="n")
    {
      y <- y+l_pkt
      z <- outer(x, y, powierzchnia) 
      persp3d(x,y,z, alpha=0.5)
    }
    else
    {
      if(jakie_rownanie=="p")
      {
        y <- y+2*l_pkt
        z <- outer(x, y, powierzchnia) 
        persp3d(x,y,z, alpha=0.5)
      }
      else 
      {
        print("nie potrafi�!")
      }
    }
  }
}
rgl.set(1)
wykres_powierzchni("fi", powierzchnia_wyniki)
wykres_powierzchni("n" , powierzchnia_wyniki)
wykres_powierzchni("p" , powierzchnia_wyniki)
wykres_powierzchni("fi", powierzchnia_wartosci_funkcji)
wykres_powierzchni("n" , powierzchnia_wartosci_funkcji)
wykres_powierzchni("p" , powierzchnia_wartosci_funkcji)
#macierz w zmiennej "y" oznacza macie� o kolumnach kt�re s� kolejno wype�nione przez 1-dynki, 2-jki, itd...
open3d()
wykres_punktowy <- function(macierz, co_chcesz_wykreslic)
{ #trzeba wpisywa� punkty pojedynczo (x1,y1,z1), (x2,y2,z2) itd. - dlatego poni�szy zapis:
  if(co_chcesz_wykreslic=="fi")
    plot3d(rep((1:l_pkt),l_pkt),matrix(rep(1:l_pkt,l_pkt),ncol=l_pkt,byrow=T),macierz[1:l_pkt,])
  else if(co_chcesz_wykreslic=="n")
    plot3d(rep((1:l_pkt),l_pkt),matrix(rep(1:l_pkt,l_pkt),ncol=l_pkt,byrow=T),macierz[(l_pkt+1):(2*l_pkt),])
  else if(co_chcesz_wykreslic=="p")
    plot3d(rep((1:l_pkt),l_pkt),matrix(rep(1:l_pkt,l_pkt),ncol=l_pkt,byrow=T),macierz[(2*l_pkt+1):(3*l_pkt),])
  else print("Nie potrafi�!")
}
rgl.set(2)
wykres_punktowy(warunki_poczatkowe,"fi")
wykres_punktowy(warunki_poczatkowe,"n" )
wykres_punktowy(warunki_poczatkowe,"p" )

#wykres punktowy wektora
{
#plot3d(rep((1:l_pkt),l_pkt),matrix(rep(1:l_pkt,l_pkt),ncol=l_pkt,byrow=T),wartosci_funkcji[1:l_pkt^2])
#plot3d(rep((1:l_pkt),l_pkt),matrix(rep(1:l_pkt,l_pkt),ncol=l_pkt,byrow=T),wartosci_funkcji[(l_pkt^2+1):(2*l_pkt^2)])
#plot3d(rep((1:l_pkt),l_pkt),matrix(rep(1:l_pkt,l_pkt),ncol=l_pkt,byrow=T),wartosci_funkcji[(2*l_pkt^2+1):(3*l_pkt^2)])
}

########sprawdzam r�no�ci z eRa poni�ej:##########
ptm <- proc.time()
proc.time()-ptm

rgl.dev.list()
rgl.set(8)
rgl.cur()

x <- c(1,-2)
y <- c(5,-3)
A <- matrix(c(2,2,-4,4), ncol = 2)
A
x-solve(A)%*%y

my_which <- function(x){which(x!=0)}
tmp <- apply(jakobian[-indeksy_kontaktow,-indeksy_kontaktow],1, my_which)
for(i in 1:7200){if(length(tmp[[i]])==0) print("wektor zerowy!")}
