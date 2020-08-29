import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import animation
from matplotlib.animation import FuncAnimation
from scipy.signal import find_peaks
from scipy.stats import linregress

"""
a ideia desse codigo é unir as funcionalidades do codigo modos3, com as animações feitas pela amanda com os resultados
que consegui obtir no codigo modosturbo, analizando a onda em uma região especifica do espaço, onde  as oscilações são observadas
"""

M = .5
raio = 2*M
L = 100
N = 5000
horizonte = raio+1e-10
inf_neg =(horizonte)+2*M*np.log((horizonte)/(2*M)-1)
"""
M = .25
raio = 2*M
L = 100
N = 5000
inf_neg = (raio+1e-10)+2*M*np.log((raio+1e-10)/(2*M)-1)
"""
# tmax = 20000
# tmax2 = 2000

# t = np.linspace(0, tmax+1, tmax+1)

rt = np.linspace(inf_neg, L, N+1) #coordenada tartaruga r*
r = np.zeros(N+1)                  #coordenada "normal" r

def biseccao(n,r0,L,tol = 1e-6):
	L+=10#isso ainda faz sentido?
	#L+10 é maximo valor que presume-se que pode ser obtido da coordenada tartaruga
	# rn = randint(int(r0),int(L))
	# print('\n chute inicial',rn)

	ra,rb,contagem =r0,L,0
	rn=.5*(ra+rb)
	erro =1
	while erro > tol and contagem < 10000:

		fa = ra + 2*M*np.log((ra/(2*M) - 1)) - rt[n]
		fb = rb + 2*M*np.log((rb/(2*M) - 1)) - rt[n]
		fn = rn + 2*M*np.log((rn/(2*M) - 1)) - rt[n]
		erro = fn

		if fn*fa<=0:
			rb = rn
			rn = .5*(rn+ra)
		elif fn*fb <0:
			ra = rn
			rn = .5*(rn+rb)
		else:
			print('não encontrada raiz no intervalo entre', ra,' e ',rb)
		contagem +=1
	return rn

r0 = horizonte
for n in range(0,N+1):
	# r[n] = tartaruga(n,r0)
	r[n] = biseccao(n,r0,L)
	r0 = r[n]

#As distâncias vertical e horizontal da nossa malha de pontos
dx = (L-inf_neg) / N
# dt = dx**2 / 2#pq?
dt = dx - 0.0005 
x0 =int(N-N/3)
# x0=800

t = np.arange(0,500,dt)
print('tamanho de t: ', len(t))
tmax= len(t)
#Definfindo o potencial 
l = 0
s=0#perturbaçao escalar
V = np.zeros(N+1)
V = (1 - 2*M/r)*(l*(l+1)/r**2 + 2*(1-s**2)*M/r**3)
Vmax = max(V)#definindo o potencial maximo

#Definindo uma Gaussiana genérica f(x) = a*np.exp(-(x-b)**2/(2*c**2)):
#a=300 V, x vai até 4.5
a = Vmax
b = rt[x0]
c = a
f = np.zeros(N+1)
f = a*np.exp(-.5*((rt-b)/(c))**2)

# Definindo a condição de contorno g(x)
g = np.zeros(N+1)

Y = np.zeros((N+1,len(t)))
# Y = np.zeros((N+1,tmax+1))

c1 = (dt**2)/(dx**2)
c2 = 2 - 2*(dt**2)/(dx**2) 

Y[1:-1,0] = f[1:-1]
Y[1:-1,1] = Y[1:-1,0] + g[1:-1]*dt



''' 
segundo o critério de instabilidade de von newmann, temos que a equação só admite solução para (dt²/dx²)sin²(O/2)+(dt²/4)V(O) < 1,]
isso equivale a impor que (dt²/dx²)+(dt²/4)Vmax <1 ou ainda, Vmax < 4/dt**2)*(1-dt**2/dx**2)

'''

for l in range(1,len(t)-1):
	Y[1:-1,l+1] = -Y[1:-1,l-1] +c1*(Y[2:,l] + Y[:-2,l]) + (c2 - V[1:-1])*Y[1:-1,l]
print(' Vmax', max(V))
print('passo dt', dt)
print('passo dx', dx,'preciso que dt<dx')

"""
#método de ordem superior
c1 = (dt**2)/(dx**2)
c2 = 2 - 2*(dt**2)/(dx**2) 

Y[1:-1,1] = f[1:-1]

Y[1:-1,0] = -g[1:-1]*dt + (c1/2)*(f[2:]+f[:-2]) + f[1:-1]*(c2 - V[1:-1])/2

for l in range(1,4):
	Y[1:-1,l+1] = -Y[1:-1,l-1] +c1*(Y[2:,l] + Y[:-2,l]) + (c2 - V[1:-1])*Y[1:-1,l]

c3 = 12*(dt**2)*V[2:-2] - 30


for l in range(2,tmax-4):
	Y[2:-2,l+2] = 16*Y[2:-2,l+1] + 16*Y[2:-2,l-1] - Y[2:-2,l-2] - c1*(-Y[4:,l] +16*Y[3:-1,l] -30*Y[2:-2,l] +16*Y[1:-3,l] -Y[:-4,l] ) + c3*Y[2:-2,l]
"""
##condicao de Von newmann
VN = (4/dt**2)*(1-dt**2/dx**2)
print('\n condicao de VN:', VN)
if VN>Vmax:print('estabilidade de von newmann satisfeita')
else: print('O método é INSTÁVEL !!!, reveja os parametros atribuidos')

#Y_t = np.absolute(Y[x0,:])
# Y_t = Y[x0,:]

#transposicao da matriz Y, para facilitar nas animações
Y = Y.T

#################################################### Animação ############################################################
'''
#primeira animacao: potencial e onda
fig1, ax1 = plt.subplots(figsize=(5, 3))
ax1.set(xlim=(inf_neg, L), ylim=(-1.5*a, 1.5*a))
ax1.plot(rt,V,'lightcoral')

line = ax1.plot(rt, Y[0, :], color='darkcyan', lw=1)[0]
point1, = ax1.plot([rt[x0]], [Y[0,x0]],marker='o', color='orange', markersize=3)
time_text1 = ax1.text(0.5, 0.95, '', transform=ax1.transAxes)
time_text2 = ax1.text(0.5, 0.88, '', transform=ax1.transAxes)

def animate(i):
    line.set_ydata(Y[i, :])
    line.set_xdata(rt)
    time_text1.set_text('frame = %.0f ' % i )
    time_text2.set_text('instante = %.0f s ' % t[i] )
    return line, time_text1,

anim = FuncAnimation(
    fig1, animate, interval=1, frames=tmax)

#segunda animacao: ponto de observacao
def animate1(i):
    point1.set_ydata(Y[i,x0] )
    return point1, 

anim1 = FuncAnimation(
    fig1, animate1, interval=1, frames=tmax)
plt.draw()
plt.show()
'''


###################################### Funcao para selecionar o frame de analise ###########################################
"""
Para determinar o intervalo de analise, uma análise manual da animação é feita, observando onde se analisar, e entre qual intervalo
de frames.

analisando a animacao parece que podemos tentar observar,
em x0=23.8678, frame final =4290, frame inicial = 3090
"""
#pontos onde a análise começa e depois termina

obs = 23.8676
ponto = (obs-inf_neg )/dx
ponto = int(ponto)

print('\n posição de observação:', ponto)
print('\n ponto de observação:', rt[ponto])

# pt_inicial =3090
pt_inicial =3390
#pt_observacao = np.where(rt==rt[x0])[0]
#tamanho = pt_observacao.size

#print('\n tamanho depois do where:',pt_observacao)
#print('\n tamanho do array pt_observacao :',tamanho)

#pt_observacao = np.array(pt_observacao)#pega o unico valor do array/espera-se que só um ponto tenha sido encontrado
#print('tamanho depois da conversao:',pt_observacao)
#pt_observacao = pt_observacao.item()
pt_final = 4290

#Y_t = Y[pt_observacao, 3090:4290]
#print("\n debug: \n")
#print('pt de observacao: ',pt_observacao)
#print('\n Y_t:',Y_t)

Y_t = Y[pt_inicial:pt_final, ponto]
Y_t = np.absolute(Y_t)
t_obs = t[pt_inicial:pt_final]
print('tamanho de Y_t: ', len(Y_t))

print('\n dados de t_reduzido:')
print('t t_reduzido inicia em t[{}] e termina em t[{}]'.format(pt_inicial,pt_final))
print('\n t_reduzido começa em t={} e termina em t={}'.format(t[pt_inicial],t[pt_final]))
##################################################################################################################################

#################################### CALCULO DAS FREQUENCIAS ################################################
#plotas xmax_i, ymax_i [ ]
"""
como analizar os modos?
achar os maximos, duas vezes o delta x entre dois maximos dá o comprimento de onda
Os Xs,Ys podem ser colocados em um fit exponencial e encontrar assim a taxa de deacaimento
Plotando X x lny, o coeficiente angular será -1/C, onde Y = e^(-CX)

#segundo kokotas, Mw = (0,37,-0.09)******atualizar isso segundo konoplaya
"""
picos,_ = find_peaks(Y_t)#encontra os indices em Y_t onde estão os picos
print("picos:",picos)
pontos = [[i,Y_t[i]] for i in picos]

"""
Até aqui o codigo encontrava varios picos onde deveria encontrar somente 1, para corrigir isso, o codigo a seguir
faz a filtragem, mantendo os picos somente quando os seus vizinhos estão a uma distância maior que 100 pontos em x.
Para frequências maiores esse valor talvez precise ser ajustado, para verificar se a filtragem está funcionado corretamente,
basta invocar a função verificar_filtragem().
"""

picost = []#lista vazia que recebe somente os picos unicos
picost = picos

'''
for i in range(0,len(picos)-1):
	if abs(picos[i+1]-picos[i])<100: pass
	else: picosx.append(picos[i])

'''

print("picost:",picost)

picosy = [Y_t[i] for i in picost]#coordena y dos picos

def verificar_filtragem():
	"""
	Essa função possui como único proposito verificar
	se o método desenvolvido para filtrar os picos,
	de forma a haver somente um em cada máximo está funcionando,
	se não estiver, será necessário ajustar o valor da filtragem
	para a frequencia trabalhada"""

	print('picos lista:',picos,'\n tamanho',len(picos))
	print('\n picos após a filtragem: ', picosx)

	plt.subplot(2,2,1)
	plt.title('antes da filtragem')
	plt.plot(t_obs,Y_t)
	plt.plot(picos,[Y_t[i] for i in picos],'ok')

	plt.subplot(2,1,2)
	plt.title('após a filtragem')
	plt.plot(t_obs,Y_t)
	plt.plot(picosx,picosy,'or')

	plt.show()

# verificar_filtragem()

def frequencia_real(xs):
	valores_medios =[]
	for i in range(0,len(xs)-1):
		Deltax = abs(xs[i]-xs[i+1])
		valores_medios.append(Deltax)
	print('valores medio de delta x: ', valores_medios)
	return 1/(2*sum(valores_medios)/len(valores_medios))#retorna o inverso de duas vezes delta x medio

def frequencia_imaginaria(xs,ys):
	lny = np.log(ys)
	inclinacao, intercept, r_value, p_value, std_err = linregress(xs,lny)#regressao linear de x vs ln y
	return inclinacao,lny,intercept

"""
para calcular a frequencia precisamos dos valores de t onde os picos ocorrem, e não somente os índices.
No geral, creio que isso não muda muita coisa pois os termos da matriz t estão igualmente espaçados.
"""
tempos_picos = [t_obs[i] for i in picost]

Wr = frequencia_real(tempos_picos)
Wi,lny,lnA= frequencia_imaginaria(tempos_picos,picosy)

print('frequencia real: ',Wr,'\nfrequencia imaginaria: ', Wi)
print('\nresultado esperado: \n')
# print('Wr: ',2*0.37,'\nWi: ', 2*0.09)

t2=np.linspace(t[pt_inicial],t[pt_final],100)

Y_exp = np.exp(lnA)*np.exp(Wi*t2)

plt.subplot(2,1,1)
plt.yscale("log")
plt.plot(t_obs,Y_t)#curva dos modos normais

plt.plot(tempos_picos,picosy,'ok',MarkerSize=2)#pontos de maximo

plt.plot(t2,Y_exp,'-r')#ajuste exponencial dos modos
plt.title('modos normais e ajuste exponencial')


plt.subplot(2,1,2)
plt.plot(tempos_picos,lny,'ok',MarkerSize=2)
plt.plot(t2,Wi*t2+lnA,'-r')#acho que isso deveria parecer uma reta

plt.xlabel('$t$')
plt.ylabel('$ln Y_t$')
plt.title('reta de regressao')

plt.show()

'''
##################################################################################################################################


# anim.save('v3.mp4') #Está dando erro pra salvar : P

#Configurando cada uma dos gráficos pedidos no problema. 
print(' Vmax', max(V))
print('passo dt', dt)
print('passo dx', dx,'preciso que dt<dx')

##condicao de Von newmann
VN = (4/dt**2)*(1-dt**2/dx**2)
print('condicao de VN:', VN)
print('é preciso ser maior que Vmax')

fig = plt.figure(figsize=(16, 10))

ax1 = fig.add_subplot(211) 

ax1.axhline(0, color='k')
ax1.axvline(0, color='k')

ax1.set_xlabel('$x$') 
ax1.set_ylabel('$f(x)$') 

ax1.plot(rt, V, 'lightcoral',label='Potencial')
ax1.plot(rt, Y[0,:], 'darkcyan',label='Gaussiana')
plt.legend()

ax2 = fig.add_subplot(212)

ax2.axhline(0, color='k')
ax2.axvline(0, color='k')
#ax2.set_yscale('log')

ax2.set_xlabel('$t$')
ax2.set_ylabel('$Y(t)$')

ax2.plot(t[3090:4290], Y_t, 'mediumspringgreen',label='Função de onde em x0')
plt.legend()

plt.tight_layout() 

plt.show() 
'''
