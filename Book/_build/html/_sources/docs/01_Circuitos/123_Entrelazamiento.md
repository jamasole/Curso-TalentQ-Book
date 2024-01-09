---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.15.2
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

+++ {"slideshow": {"slide_type": "slide"}}

# Entrelazamiento en acción



$ \newcommand{\bra}[1]{\langle #1|} $
$ \newcommand{\ket}[1]{|#1\rangle} $
$ \newcommand{\braket}[2]{\langle #1|#2\rangle} $
$ \newcommand{\i}{{\color{blue} i}} $ 
$ \newcommand{\Hil}{{\mathcal H}} $
$ \newcommand{\cg}[1]{{\rm C}#1} $
$ \newcommand{\bn}{{\bf n}} $

+++ {"slideshow": {"slide_type": "skip"}}


- [Desigualdades de Bell](#Bell)

    - [*Teorema CHSH*](#csch)


- [El experimento GHZ](#ghz)


- [Teleportación](#teleportacion)
    - [LOCC](#LOCC)


- [Codificación Superdensa](#superdense)

- [Intercambio de entrelazamiento](#intercambent)

```{code-cell} ipython3
---
slideshow:
  slide_type: '-'
---
import sys
sys.path.append('../')
import macro_tQ as tQ

import numpy as np
import scipy.linalg as la
from IPython.display import display,Markdown,Latex
import matplotlib.pyplot as plt
from qiskit.tools.visualization import array_to_latex
```

+++ {"slideshow": {"slide_type": "slide"}}

##  Desigualdades de Bell

+++ {"slideshow": {"slide_type": "-"}}

Las teorías con *realismo local* asumen que los valores que adquieren las magnitudes que se miden en un experimento *pertenecen* al sistema medido. 

+++ {"slideshow": {"slide_type": "-"}}

Así, en una teoría con realismo local la posición de una partícula es algo bien definido aunque no lo estemos midiendo. 

+++ {"slideshow": {"slide_type": "fragment"}}

La palabra *local* hace referencia a que ningún agente puede propagar su acción a mayor velocidad que la luz. Podría usarse la palabra *causal* en su lugar.

+++ {"slideshow": {"slide_type": "-"}}

Por ejemplo, al hablar del espín de un electrón, no es correcto decir que *la proyección del espín a lo largo del eje $\hat{\bf z}~$ <u>es</u> $~+\hbar/2$*.

+++ {"slideshow": {"slide_type": "-"}}

 Lo correcto es decir que, <u>al medir</u>,  cualquiera de los valores $\pm\hbar/2$ se <u>adquiere</u> o <u>pone de manifiesto</u>  de forma aleatoria.

+++ {"slideshow": {"slide_type": "slide"}}

En 1964  físico nor-irlandés John Bell, que trabajaba en el CERN,   [demostró](https://cds.cern.ch/record/111654/files/vol1p195-200_001.pdf) que **todas** las teorías que respetan el realismo local satisfacen ciertas *desigualdades matemáticas*. 

En particular, <u>todas las magnitudes que evolucionan siguiendo las leyes de la física clásica</u> las satisfacen.

+++ {"slideshow": {"slide_type": "-"}}

Por el contrario, John Bell mostró cómo, en Mecánica Cuántica, el entrelazamiento, incorpora *correlaciones* sutiles entre las medidas que permiten traspasar  dichas desigualdades. 


+++ {"slideshow": {"slide_type": "skip"}}

La discusión, de puramente filosófica a ser objeto de investigación experimental, culminó con el [experimento](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.49.1804) de Alain Aspect y colaboradores en 1982.

Se observó que la Mecánica Cuántica viola las desigualdades de Bell y, por tanto, **no es una teoría con realismo local**. 

+++ {"slideshow": {"slide_type": "fragment"}}

La propuesta de John Bell dio pie a una familia de desigualdades que ponen en evidencia la imposibilidad de obtener ciertas correlaciones en un mundo clásico. 

Vamos a examinar la desigualdad en la forma estudiada por Clauser, Horne, Shimony y Holt [CHSH](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.23.880).  


+++ {"slideshow": {"slide_type": "skip"}}

Posteriormente estudiaremos el experimento de GHZ, el cuál, también pone de manifiesto las correlaciones sutiles que introduce el entrelazamiento de una forma mucho determinista, en lugar de estadística.

+++ {"slideshow": {"slide_type": "slide"}}

## Perfecta anticorrelación

Central en esta discusión es la presencia de entrelazamiento. Vamos a seleccionar el denominado *singlete*  de la Base de Bell.
<br>
<br>
$$
\ket{B_{11}}=\frac{1}{\sqrt{2}}(\ket{01}-\ket{10})
$$

Supongamos que  Alice y Bob poseen sendos medidores de Stern Gerlach apuntando en la dirección $\hat{\bf z} = (0,0,1)$.

+++ {"slideshow": {"slide_type": "fragment"}}

- Si Alice registra +1 el estado colapsa a $\ket{01}$ y, por tanto, Bob solo podrá medir  $\, -1$


-  Si Alice registra -1 el estado colapsa a $\ket{10}$ y, por tanto, Bob solo podrá medir  $\, +1$

+++ {"slideshow": {"slide_type": "fragment"}}

En definitiva, hay una anticorrelación perfecta que se pone de manifiesto en el valor medio del producto de las medidas. 

Si las mediciones de Alice y Bob son $a_i=\pm 1\Rightarrow b_i=\mp 1$ respectivamente, el valor medio
de la variable aleatoria $\{a_i, b_i\}$ es $-1$

$$
\overline{ab} = \frac{1}{N}\sum_{i=1}^N a_i b_i =\frac{1}{N}\sum_{i=1}^N (-1) = -1
$$

+++ {"slideshow": {"slide_type": "slide"}}

Vamos a ver cómo la predicción teórica confirma este hecho.

El estado $\ket{B_{11}}$ ya es autoestado del observable asociado 

$$
Z\otimes Z \ket{B_{11}} = Z\ket{0}Z\ket{1} - Z\ket{1}Z\ket{0} = -\ket{01} + \ket{10} = -\ket{B_{11}}
$$

+++ {"slideshow": {"slide_type": "fragment"}}

Y el valor esperado satura  este autovalor 
<br>

$$
 \langle Z\otimes Z\rangle = \bra{B_{11}}Z\otimes Z\,  \ket{B_{11}}  = -\braket{B_{11}}{B_{11}}= -1\, 
$$

+++ {"slideshow": {"slide_type": "skip"}}

<div class="alert alert-block alert-danger">
<b>Notar</b>
<br>
La anti correlación que hemos hallado es la misma que en un experimento clásico preparado en una bolsa que contiene dos calcetines de dos colores: si Alice saca el blanco, el que saca Bob tiene que ser negro.     
    
</div>


+++ {"slideshow": {"slide_type": "slide"}}

La cosa se pone más divertida cuando los dos polarizadores de Stern Gerlach <u>*no se orientan en la misma dirección*</u>

El observable asociado ahora a la dirección $~\hat{\bf n}~$  será  $~ A = \hat{\bf n}\cdot \boldsymbol{\sigma}$

Aun así, los autovalores de este operador y, por ello, la proyección del espín seguirá siendo $\pm 1$

<div class="alert alert-block alert-info",text-align:center>
<p style="text-align: left ;color: navy;">  
<b>Teorema</b>: 
<br>
el valor medio del producto de las proyecciones de espín a lo largo de sendos ejes $~\hat{\bf m}~$ y $~\hat{\bf n}~$ viene dada por el coseno del ángulo $\theta$  que forman los  ejes de los dos detectores 
<br>
<br>    
$$
\bra{B_{11}}(\hat{\bf m}\cdot \boldsymbol{\sigma}\otimes \hat{\bf n}\cdot \boldsymbol{\sigma})\ket{B_{11}}
 = -\cos \theta = - \hat{\bf m}\cdot \hat{\bf n}
$$
</p>
</div>


+++ {"slideshow": {"slide_type": "skip"}}

<div class="alert alert-block alert-success">
<b>Ejercicio:</b> 
prueba este resultado.
</div>

+++ {"slideshow": {"slide_type": "slide"}}

- Cuando los ejes son paralelos recuperamos la  anti-correlación, en cualquier dirección

$$
-\hat{\bf n}\cdot \hat{\bf n} = - \cos 0 = -1 
$$

- cuando las direcciones de los detectores de Alice y Bob no coinciden sigue habiendo una medida cuántica $a,b = \pm 1$ en ambos. 
<br>
<br>
- Clásicamente, ¿qué esperarías medir?



+++ {"slideshow": {"slide_type": "slide"}}

### La desigualdad de CHSH 

Enn 1970 Clauser, Horne, Shimony y Holt [CHSH](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.23.880) propusieron una figura de mérito capaz de distinguir las correlaciones cuánticas de las clásicas.

+++ {"slideshow": {"slide_type": "fragment"}}

La idea es que Alice y Bob pueden escoger, cada uno,   <u>*dos orientaciones  arbitrarias*</u> para su detector.

Para el detector de Alice las denotamos $\hat{\bf n}_A, \hat{\bf n}'_A$  y  para el de Bob $\, \hat{\bf n}_B, \hat{\bf n}'_B.$ 

Los pasos a seguir son los siguientes:

**1.-** Alice y Bob seleccionan *una orientación* para sus detectores. Hay cuatro parejas posibles 

$$
\begin{array}{c|c} {\rm Alice} & {\rm Bob} \\ \hline \hat{\bf n}_A  &\hat{\bf n}_B \\ \hat{\bf n}_A &\hat{\bf n}'_B \\ \hat{\bf n}'_A & \hat{\bf n}_B \\ \hat{\bf n}'_A &\hat{\bf n}'_B \end{array}
$$

+++ {"slideshow": {"slide_type": "slide"}}


**2.-** Alice y Bob reciben un electrón cada uno de un par entrelazado en el estado $\ket{B_{11}}$


**3.-** Alice y Bob realizan la medida de la proyección del espín a lo largo del eje elegido y anotan el resultado de la medición $(a,b)=(\pm 1, \pm 1)$ 

+++ {"slideshow": {"slide_type": "fragment"}}

**4.-** Repiten el paso anterior un número $i=1,..., N$ grande de veces, y con los datos obtenidos $~(a_i,b_i)=(\pm 1, \pm 1)~$ pueden reconstruir la correlación

$$
 C(\hat{\bf n}_A, \hat{\bf n}_B) = \frac{1}{N}\sum_{i=1}^N a_i b_i ~\in~  [-1,1]
$$

+++ {"slideshow": {"slide_type": "fragment"}}

**5.-** Repiten todo el proceso anterior para las cuatro posibles orientaciones elegidas de forma aleatoria. Con las  $4N$ mediciones construyen la cantidad
<br>
<br>
$$
R = | C(\hat{\bf n}_A, \hat{\bf n}_B) +  C(\hat{\bf n}_A, \hat{\bf n}'_B) +  C(\hat{\bf n}'_A, \hat{\bf n}_B)-  C(\hat{\bf n}'_A, \hat{\bf n}'_B)|
$$


+++ {"slideshow": {"slide_type": "slide"}}

En un mundo clásico, supondríamos que los valores $a_i,b_i$ proceden de *valores predefinidos* para cada sistema individual, sobre el que efectuamos simplemente un promedio estadístico de muchos sistemas. Entonces podemos probar la **desigualdad de CHSH**

<div class="alert alert-block alert-info",text-align:center>
<p style="text-align: left ;color: navy;">  
<b>Teorema</b>: 
La desigualdad de CHSH afirma que 
$$
R \leq 2
$$
</p></div>

<details>
    <summary><p style="text-align:right"> >><i>Prueba</i> </p></summary>
Es fácil ver que se cumple para cada colección $a_i,a'_i,b_i,b'_i \in \pm 1$ la desigualdad
<br>
<br>
$$
a_i (b_i +  b'_i) + a'_i (b_i -  b'_i) = \pm 2
$$
<br>
porque si $b_i +  b'_i=\pm 2$ entonces $b_i - b'_i=0$ y viceversa. Ahora podemos demostrar la desigualdad 
\begin{eqnarray}
R &=& \lim_{N\to \infty}|\frac{1}{N} \sum_{i=1}^N \left( a_i b_i + a_i b'_i + a'_i b_i - a'_i b'_i\right) | \nonumber\\
&=& \lim_{N\to \infty}\frac{1}{N}| \sum_{i=1}^N  \left( a_i (b_i +  b'_i) + a'_i (b_i -  b'_i)\right)| \nonumber\\
&\leq &  \lim_{N\to \infty} \frac{1}{N} \sum_{i=1}^N  |\left( a_i (b_i +  b'_i) + a'_i (b_i -  b'_i)\right)| \ \nonumber\\
&= &  \lim_{N\to \infty} \frac{1}{N} \sum_{i=1}^N  2 \ \nonumber\\
&=& 2
\end{eqnarray}
</details>   

+++ {"slideshow": {"slide_type": "fragment"}}

La Mecánica Cuántica nos proporciona una respuesta teórica para $R$ que sólo depende de los ángulos relativos $\cos\theta_{AB} = \cos(\theta_A-\theta_B)= \hat{\bf n}_A\cdot \hat{\bf n}_B$.
<br>
<br>

$$
R = | \cos\theta_{AB}  + \cos \theta_{A'B} + \cos \theta_{AB'}  - \cos \theta_{A'B'}|\, .
$$
<br>

+++ {"slideshow": {"slide_type": "slide"}}

Ahora sólo hace falta jugar un poco con los detectores.

Por ejemplo, podemos situarlos en el plano $(y,z)$, perpendicular al eje de propagación $x$, de manera que los vectores  $\hat{\bf n}'_A,\hat{\bf n}_{A},\hat{\bf n}_B$ y $\hat{\bf n}_B'$ estén ordenados correlativamente en sentido horario. Finalmente tomaremos
dos ejes coincidentes $\hat{\bf n}_{A}=\hat{\bf n}_B$ paralelos $\Rightarrow \theta_{AB}=0$, y apertura igual para ambos, $\theta_{A'A} = \theta_{BB'} = \varphi$, de modo que $\theta_{A'B'} = 2\varphi$ 

+++ {"slideshow": {"slide_type": "-"}}

<center> 
<img src="./images/CHSH_basis.png" width='35%'/>
</center>

+++ {"slideshow": {"slide_type": "fragment"}}



$$
R = |1 +  2 \cos\varphi  - \cos 2\varphi|\, .
$$

Derivando vemos que esta expresión alcanza su máximo cuando $\sin\varphi = \sin 2\varphi$ lo cual tiene solución $\varphi = \pi/3 = 60^\circ$. Sustituyendo encontramos $R= 2.5>2$. 
 

+++ {"slideshow": {"slide_type": "slide"}}

Vamos a hacer una lista de ángulos relativos entre los dos polarizadores 

```{code-cell} ipython3
---
slideshow:
  slide_type: '-'
---
'angulos de medida'
phi_divs = 100
phi_list=np.linspace(0,2*np.pi,phi_divs) #lista de angulos'

'lista de correlaciones C[i] a calcular'
C=[0,0,0,0]

'lista de valores de R'
R=np.zeros(phi_divs)

'Numero de medidas'
nshots = 2048

from qiskit import QuantumCircuit, Aer, execute
M_simulator = Aer.get_backend('qasm_simulator')

```

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
for j in range(phi_divs): 

    'ángulos de medida'
    phi=phi_list[j] 
    '''A' = 2\phi, A = B = phi, B' = 0'''
    angles_AB=[[phi,phi],[2*phi,phi],[phi,0],[2*phi,0]] # AB , A'B , AB', A'B'

    for i in range(4):        
        'una pareja de ángulos uno para A y otro para B '
        ang_AB=angles_AB[i]

        'hay un circuito para cada proceso de medida'
        qc=QuantumCircuit(2,2)
        'creamos el par de Bell B(11)'
        qc.x(0)
        qc.x(1)
        qc.h(0)
        qc.cx(0,1)

         
        'simulamos la medida en la base de los polarizadores de A y B rotados en torno al eje x'
        qc.rx(-ang_AB[0],0)  # notar el valor negativo del ángulo de rotación 
        qc.rx(-ang_AB[1],1)
        qc.measure([0,1],[0,1])
        
        if j ==0 and i ==0:
            display(qc.draw('mpl'))
        ' medimos '
        counts=execute(qc,backend = M_simulator,shots = nshots).result().get_counts()         
        'construimos el correlador'
        C[i]= 0 
        for bitstring, counts in counts.items():
            C[i] += (-1)**(sum([int(bit) for bit in bitstring])) * counts/nshots 

    'construimos la cantidad R'
    R[j]=np.abs(C[0]+C[1]+C[2]-C[3])
```

+++ {"slideshow": {"slide_type": "slide"}}

Podemos graficarlo y compararlo con el valor analítico

```{code-cell} ipython3
---
slideshow:
  slide_type: '-'
---
' función analítica '
fx= np.abs(1 + 2*np.cos(phi_list) - np.cos(2*phi_list))
plt.plot(phi_list,fx,'k-',linewidth=1)

' resultado de la simulación'
plt.plot(phi_list,R,'.')
plt.axhline(y = 2, color = 'r', linestyle = '-')
plt.axhline(y = -2, color = 'r', linestyle = '-')

' líneas horizotales en +2 y - 2'
plt.axvline(x = np.pi/3, color = 'b', linestyle = '-')
plt.axvline(x = 5*np.pi/3, color = 'b', linestyle = '-')
plt.xlabel(r'$\varphi$')
plt.ylabel(r'$R(\varphi)$')
plt.show()
```

+++ {"slideshow": {"slide_type": "skip"}}

Prueba con otros estados de la base de Bell. Realiza este experimento en un ordenador real.

+++ {"slideshow": {"slide_type": "skip"}}

<div class="alert alert-block alert-success">
<b>Ejercicio  1.2.3.2</b> 
    
Usar bases perpendiculares de  Alice $A\perp A'$,  y Bob $B\perp B'$, y formando un ángulo $\varphi$ entre sí. 
    
Sin pérdida de generalidad puedes tomar $B = X\perp Z = B'$. 
    
Variar $\varphi$ en el intervalo $(0,\pi)$ y hallar el valor de la máxima violación de la desigualdad de Bell.
    
¿Es mayor que en el ejemplo anterior? ¿Podrías  hallar otra violación aún mayor?
</div>

+++ {"slideshow": {"slide_type": "slide"}}

## Experimento GHZ

+++ {"slideshow": {"slide_type": "-"}}

Las desigualdades de Bell-CHCS demuestran que hay una manera de distinguir una teoría con  realismo local de una en la que éste postulado no  exista.

Lo hacen de una forma estadística: es necesario formar ciertas medias de datos sobre estados de 2-cúbits, los  estados de Bell. 

+++ {"slideshow": {"slide_type": "fragment"}}

En 1997 Greenberger, Horne y Zeilinger estudiaron las propiedades de una serie de estados de 3-cúbits:
<br>
<br>
$$
\ket{\rm GHZ} = \frac{1}{\sqrt{2}}( \ket{000} - \ket{111} )\, . 
$$

Es muy fácil ver que, en vez de una violación **probabilística** de una *desigualdad*, el estado GHZ conduce a una  violación **determinista** de una *igualdad*. 

+++ {"slideshow": {"slide_type": "slide"}}

Supongamos que:  
<br>
-  un sistema está  compuesto de  tres subsistemas, $A, B$ y $C$. 
<br>
<br>
-  en cada uno de ellos hay dos *magnitudes* observables, $X$ e $Y$, que al ser medidas adquieren valores binarios $x, y = \pm 1$. 

+++ {"slideshow": {"slide_type": "fragment"}}

Pretendemos saber:

si existe un estado del sistema tal que, con el resultado  de  medidas simultáneas $X$ o $Y$ de cada una de sus partes $A,B$ y $C$ , se obtengan resultados $x$ e $y$ que verifiquen las siguientes ecuaciones:

$$
\begin{array}{ccc}
{\rm medimos} & & \hbox{obtenemos} ~ x,y ~~ \hbox{tales que} \\
XYY ~~~~~ &\to& ~xyy =~~~ 1\nonumber\\
YXY ~~~~~ &\to& ~yxy =~~~ 1 \label{GHZeq}\\
YYX ~~~~~ &\to& ~yyx =~~~ 1\nonumber\\
XXX ~~~~~ &\to& ~xxx = -1\nonumber
\end{array}
$$

+++ {"slideshow": {"slide_type": "slide"}}

Si la teoría satisface los axiomas de realismo local,  los valores de $X$ es $Y$ estarán bien definidos independientemente de la medida que efectuemos. 

Dicho de otra forma, las medidas que efectuamos son compatibles y no afectan a los resultados posibles.

+++ {"slideshow": {"slide_type": "fragment"}}

Entonces podemos utilizar conclusiones que extraigamos de los resultados de un experimento para los de otro (notar que se supone que el estado es el mismo en todos los casos  )

- Las tres primeras afirman que $x=1$ puesto que $y^2=1$. 

- La última afirma que $x=-1$

vemos que hay una contradicción. Concluimos que las medidas efectuadas *no son compatibles entre sí*. 

+++ {"slideshow": {"slide_type": "slide"}}

Como vemos, la paradoja existe en tanto en cuanto atribuyamos a las magnitudes $X$ e $Y$ una noción de realidad  independiente de la medición. 

Sólo de esta manera podemos inferir, de las tres primeras líneas, el resultado  $xxx  =  1$ sin realizar una medición directa de $XXX$  por consistencia.


+++ {"slideshow": {"slide_type": "-"}}

La satisfacción, por tanto, de las cuatro afirmaciones,  pasa por aceptar que las variables $X$ e $Y$, no pueden tener valores predefinidos simultáneamente  en los cuatro experimentos. 


+++ {"slideshow": {"slide_type": "fragment"}}

Cuánticamente es fácil ver que el estado GHZ proporciona una solución al conjunto de ecuaciones.

+++ {"slideshow": {"slide_type": "skip"}}


Supongamos que  $x$ e $y$ son los resultados de medir $X = \sigma_x~$ y $~Y = \sigma_y$
<br>
<br>

$$
X \ket{0} = \ket{1}\, , ~ X \ket{1} = \ket{0}\, ~~ , ~~ Y \ket{0} = i \ket{1}\, , ~ Y \ket{1} = -i \ket{0}\, . ~ 
$$

+++ {"slideshow": {"slide_type": "skip"}}

 Por un lado, 
<br>
<br>
$$xyy\ket{{\rm GHZ}}  = X \otimes Y\otimes Y \big(\frac{1}{\sqrt{2}}(\ket{000}-\ket{111}\big) = \frac{1}{\sqrt{2}} \big(i^2 \ket{111}- (-i)^2
\ket{000}\big) = +  \ket{{\rm GHZ}}\, ,$$

y,
análogamente, obtenemos $xyy=yxy=yyx= +1$.

+++ {"slideshow": {"slide_type": "skip"}}

Por otro, $xxx=-1$ se sigue de


$$  
xxx  \ket{{\rm GHZ}} =   X \otimes X \otimes X \big(\frac{1}{\sqrt{2}}(\ket{000}-\ket{111}\big) =
\frac{1}{\sqrt{2}}(\ket{111}-\ket{000}) =  - \ket{\rm GHZ}\, .
$$


+++ {"slideshow": {"slide_type": "slide"}}

Inicializamos el estado GHZ 

+++ {"slideshow": {"slide_type": "-"}}

<center> <b> GHZ</b>
<img src="./images/GHZ.png" width='35%'/>
</center>

```{code-cell} ipython3
---
run_control:
  marked: false
slideshow:
  slide_type: fragment
---
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, Aer, execute

qc_GHZ = QuantumCircuit(3,3) 

qc_GHZ.x(0)
qc_GHZ.h(0)
qc_GHZ.cx(0,1)
qc_GHZ.cx(0,2)

qc_GHZ.draw(output='mpl')
```

+++ {"slideshow": {"slide_type": "slide"}}

Analicemos el valor medio de alguno de los operadores $YYX, YXY, XYY$ ó $XXX$

```{code-cell} ipython3
---
slideshow:
  slide_type: fragment
---
from qiskit import Aer, execute

M_backend = Aer.get_backend('qasm_simulator')
shots=1000


'elige uno de los términos siguientes'
#multibasis = 'YYX'
#multibasis = 'YXY'
#multibasis = 'XYY'
multibasis = 'XXX'

qc_GHZ.barrier()
mc_add_measure_XYZ(qc_GHZ, multibasis)
qc_GHZ.draw(output='mpl')
```

Ejecuta el código anterior para las 4 posibles cadenas de Pauli

```{code-cell} ipython3
---
slideshow:
  slide_type: fragment
---
cuentas_GHZ = execute(qc_GHZ,M_backend,shots=shots).result().get_counts()
print('<'+str(multibasis)+'> =', mc_val_esp_sigma(cuentas_GHZ)[0])
```

+++ {"slideshow": {"slide_type": "slide"}}

## Teleportación

+++ {"slideshow": {"slide_type": "-"}}

El entrelazamiento conlleva un tipo  nuevo de correlación que acaba constituyendo un recurso importante. 

El protocolo de teleportación fue propuesto por [C.Bennet,G. Brassar, C. Crepeau, R.Josza, A. Peres u W. Wooters](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.70.1895) en 1993.


+++ {"slideshow": {"slide_type": "slide"}}

Supongamos que Alice y Bob tiene dos cúbits que se encuentran en el estado entrelazado $\ket{B_{00}}$


 Alice tiene, un segundo cúbit inicializado en un estado arbitrario $\ket{\phi}$, y se plantea la posibilidad de transferirlo o clonarlo en el laboratorio de Bob.

El circuito de la imagen permite efectuar esa tarea

+++ {"slideshow": {"slide_type": "-"}}

<div>
<img src="images/teleportacion.png" width='60%' style="margin:auto" />
</div>

+++ {"slideshow": {"slide_type": "slide"}}

- El estado inicial es
\begin{eqnarray}
\ket{\phi} \ket{B_{00}} &=&  \ket{\phi}\otimes \frac{1}{\sqrt{2}}\left(\ket{0}\otimes \ket{0} + \ket{1}\otimes\ket{1}\rule{0mm}{4mm}\right)\, \\
&=& \frac{1}{\sqrt{2}} (a \ket{000} + b\ket{100} + a\ket{011} + b\ket{111})
\end{eqnarray}
<br>
<br>
Alice tiene acceso a los dos primeros qubits y Bob al tercero. 

+++ {"slideshow": {"slide_type": "fragment"}}

<br>

- Alice realiza una  <i>medida de Bell</i> a sus dos qubits, intercalando un operador $U_{\rm desent}= (H \otimes I)\cdot U_{\rm CNOT}~$ que rota a la base de Bell
<br>
<br>
\begin{eqnarray}
 (H \otimes I\otimes I)(U_{\rm CNOT}\otimes I)\, \ket{\phi}\ket{B_{00}} & =&  \frac{1}{2}
\left[\rule{0mm}{3mm} \ket{00}(a\ket{0} + b\ket{1}) + \ket{01}(a\ket{1}+ b\ket{0}) \right. \nonumber\\
\rule{0mm}{14mm}
&&\left.  \rule{0mm}{3mm}+ ~\ket{10}(a\ket{0}-b\ket{1}) + \ket{11}(a\ket{1}-b\ket{0})\right] \label{teleport}
\end{eqnarray}

+++ {"slideshow": {"slide_type": "slide"}}


- Alice mide el estado que obra en su poder, y obtiene un 2-bit  clásico, $xy$ de manera equiprobable para las 4 posibilidades.
<br>
De forma correlacionada, el cúbit
 de Bob colapsa a uno de entre 4 estados $\ket{\phi_{xy}},~$ pero no sabe a cuál.

+++ {"slideshow": {"slide_type": "fragment"}}

- Para resolver esta ambigüedad, Alice envía el resultado de su medida $xy$ por un canal clásico a Bob

+++ {"slideshow": {"slide_type": "slide"}}

-  Bob efectúa sobre su cúbit una operación controlada por este 2-bit, $U_{xy} =  X^y Z^x $. 
<br>
<br>
$$
xy ~=~ \left\{ \begin{array}{c} 00 \\ 01 \\ 10 \\ 11 \end{array} \right\}~~ \Longrightarrow  
~~~ X^y Z^x \ket{\phi_{xy}} ~=~ \left\{  \begin{array}{rl}    I :& (a\ket{0} + b\ket{1})  \\    X: & (a\ket{1} + b\ket{0})  \\  Z:& (a\ket{0} - b\ket{1})  
\\  XZ = -iY:&   (a\ket{1} - b\ket{0}) \\
       \end{array} \right. ~~ \longrightarrow  ~~ a \ket{0} + b\ket{1}\,  = \, \ket{\phi}
$$
<br>
Como resultado de esta operación, el cúbit de Bob es finalmente $\ket{\phi}$.

+++ {"slideshow": {"slide_type": "slide"}}

Vamos preparar el estado de Alice

```{code-cell} ipython3
---
slideshow:
  slide_type: '-'
---
# escoge dos ángulos theta y phi en la esfera de Bloch
[theta,phi]=[1.1 , 2.2]

# o hazlo de forma aleatoria
#import random
#[theta,phi] = np.array([np.pi*np.random.rand(),2*np.pi*np.random.rand()])

estado = [np.cos(theta/2),(np.exp(1j*phi))*np.sin(theta/2)]

#visualicémoslo en la esfera de Bloch
from qiskit.visualization import plot_bloch_multivector, array_to_latex  
plot_bloch_multivector(estado)
```

+++ {"slideshow": {"slide_type": "slide"}}


Lo primero que haremos será incializar el circuito, dotado de tres cubits y dos canales clásicos en dos registros. Alice tendrá control sobre los dos primeros $q_0, q_1$ y Bob sobre el tercero $q_2$.

```{code-cell} ipython3
---
run_control:
  marked: false
slideshow:
  slide_type: '-'
---
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, Aer, execute

qrA = QuantumRegister(2, name="qAlice") 
qrB =  QuantumRegister(1, name="qBob") 
crx = ClassicalRegister(1, name="crx") 
cry = ClassicalRegister(1, name="cry") 
```

+++ {"slideshow": {"slide_type": "-"}}

generamos el circuito de teleportación

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
# inicializa el circuito
qc = QuantumCircuit(qrA,qrB, crx, cry)

# Alice y Bob comparten un par de Bell
qc.h(qrA[1])
qc.cx(qrA[1],qrB[0])
qc.barrier()

# inyecta en el circuito el estado a teleportar
qc.initialize(estado,0)

qc.barrier()

# Añade un medidor de estados de Bell en los cúbits de Alice y efectúa las medidas
qc.cx(qrA[0],qrA[1])
qc.h(qrA[0])
qc.barrier()
qc.measure([qrA[0],qrA[1]],[0,1])

# Añade el operador X^y Z^x controlado clásicamente por las salidas x e y 
qc.x(qrB[0]).c_if(cry, 1)  
qc.z(qrB[0]).c_if(crx, 1)  
qc.barrier()

#visualiza el circuito 
qc.draw('mpl')
```

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
# ejecuta varias veces y observa el resultado 
S_simulator = Aer.get_backend('statevector_simulator')
out_vector=execute(qc,S_simulator).result().get_statevector()
plot_bloch_multivector(out_vector)
```

+++ {"slideshow": {"slide_type": "fragment"}}

Observa al ejecutar repetidas veces que los resultados de las medidas intermedias son aleatorios, pero el estado teleportado es siempre el correcto

+++ {"slideshow": {"slide_type": "slide"}}


El siguiente ejercicio ilustra el principio de la *medida diferida*. 


<div class="alert alert-block alert-success">

<b>  Ejercicio </b> 
    
Modifica y ejecuta  el circuito de teleportación de dos formas distintas
<br>   
a)  sustituyendo los controles clásicos por controles cuánticos
<br>
<br>
b) permutando el orden de los controles y los aparatos de medida
<br>
<br>    
Discute la sutileza que distingue estas posibilidades.    
    <i>Lectura recomendada: Nielsen-Chuang sección 4.4.  </i>
</div>


+++ {"slideshow": {"slide_type": "fragment"}}

<div class="alert alert-block alert-success">

<b>  Ejercicio </b> 
    
cambia el estado que comparten Alice y Bob por $\ket{B_{11}}$ y modifica el circuito para que teleporte igualmente.

</div>


+++ {"slideshow": {"slide_type": "slide"}}

<div class="alert alert-block alert-danger">
El protocolo de teleportación parece poner en riesgo conceptos fundamentales. 
<br>

<b> No clonación: </b>
el protocolo tiene como ingrediente esencial la medida y, por tanto, la <i>destrucción</i> del estado de Alice. Como consecuencia, el estado ha sido teleportado pero no clonado. Es decir,   <i>no se viola el Teorema de No Clonación </i>. 
<br>
<br>
    
<b> Causalidad: </b>
el protocolo de teleportación <i> no viola causalidad </i>. 

Es necesario mandar una información clásica (como muy rápido a la velocidad de la luz) para resolver la ambigüedad que le queda a Bob. 

Esta parte es la que hace que la teleportación no sea un proceso instantáneo de acción a distancia.    
    
    
</div>

+++ {"slideshow": {"slide_type": "slide"}}

## Codificación superdensa

+++ {"slideshow": {"slide_type": "-"}}

Una segunda aplicación de las correlaciones que introduce el entrelazamiento permite codificar 2 bits clásicos en términos de un sólo cúbit.

Este protocolo fue propuesto por [C. Bennet y S.J. Wiesner en 1992](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.69.2881)

Supongamos que Alice quiere transmitir  un para de bits $xy\in \{00,01,10,11\}$ a Bob.
Ambos deben estar en posesión de un estado entrelazado que tomaremos $\ket{B_{00}}$

El circuito siguiente muestra el proceso 

+++ {"slideshow": {"slide_type": "slide"}}

<div>
<img src="images/superdense_coding.png" width="60%" style="margin:auto"/>
</div>

+++ {"slideshow": {"slide_type": "slide"}}

- Alice comienza efectuando una operación en su laboratorio sobre el cúbit que maneja, que depende del 2-bit $xy$ que quiere enviar
<br>
$$
U_{xy} =Z^x X^y ~~:~ ~\ket{B_{00}} ~\Rightarrow  ~\ket{B_{xy}} = \left\{ \begin{array}{c}  \frac{1}{\sqrt{2}}(\ket{00} + \ket{11})  = \ket{B_{00}} \\  \frac{1}{\sqrt{2}}(\ket{10} + \ket{01}) =\ket{B_{01}}  \\   \frac{1}{\sqrt{2}}(\ket{00} - \ket{11})= \ket{B_{10}}  \\  \frac{1}{\sqrt{2}}(\ket{10} - \ket{01}) =\ket{B_{11}} \end{array} \right. 
$$
<br>
Una vez efectuada esa operación, el estado que comparten es  $\ket{B_{xy}}$.

+++ {"slideshow": {"slide_type": "skip"}}

- Alice envía a Bob su parte del estado entrelazado a través de un canal cuántico. 
<br>

- Ahora Bob está en posesión del estado $\ket{B_{xy}}$ y sólo necesita hacer una *medida de Bell* sobre el estado que posee para recuperar la etiqueta $xy$

+++ {"slideshow": {"slide_type": "slide"}}

## Intercambio de entrelazamiento

+++ {"slideshow": {"slide_type": "skip"}}

Ya hemos visto cómo el ejemplo más sencillo de entrelazamiento tiene una aplicación muy interesante en la teleportación. Vamos a ver que el entrelazamiento se puede *contagiar* a terceras partes. En inglés se denomina *entanglement swapping*.


+++ {"slideshow": {"slide_type": "-"}}

Consideremos el siguiente circuito

<br>
<br>
<div>
<img src="images/entanglement_swap.png" width="50%" style="margin:auto"/>
</div>
<br>

+++ {"slideshow": {"slide_type": "skip"}}

- Alice (A) entrelaza un cúbit con Charles (C) y otro con Bob (B). 
<br>


+++ {"slideshow": {"slide_type": "skip"}}

- Después, Alice hace una medida de Bell, y comunica el resultado $xy$ a Charles y a Bob respectivamente, 
<br>


+++ {"slideshow": {"slide_type": "skip"}}


- los cuales efectúan las puertas controladas $Z^x$ y $X^y$ respectivamente.

+++ {"slideshow": {"slide_type": "skip"}}

-  El estado a la salida es $\ket{B_{00}}_{CB}$ de Bell compartido por Bob y Charles.

+++ {"slideshow": {"slide_type": "skip"}}

<div class="alert alert-block alert-success">

<b> Ejercicio  </b> 
    
1-  programa el circuito de la figura    

2-  ejecuta varias veces el circuito y muestra que el estado final que comparten Bob y Charles está entrelazado. ¿Es siempre el mismo estado?
    
3- Completa el circuito para que sea  capaz de teleportar un cúbit arbitrario entre Charles y Bob
 
</div>

+++ {"slideshow": {"slide_type": "skip"}}

<a id='LOCC'></a>
## LOCC

Bajo la denominacióm de Operaciones Locales y Comunicación Clásica se denomina un método en Teoría Cuántica de la Información, mediante el cual,

- una *operación local* es efectuada en una parte del sistema. 
<br>

- el resultado es comunicado clásicamente a otra del sistema
<br>

- en esa otra parte se fectúa también una operación local, condicionada al resultado comunicado clásicamente.

+++ {"slideshow": {"slide_type": "skip"}}

Al algoritmo de teleportación es un ejemplo de LOCC en el cuál 

- Alice realiza una medida de Bell
<br>

- Bob actúa localmente con $X$ o $Z$ dependiendo del resultado comunicado clásicamente. 

+++ {"slideshow": {"slide_type": "skip"}}

En el intercambio de entrelazamiento vemos que LOCC se puede utilizar para transformar un estado entrelazado en otro estado entrelazado.

En este caso un estado de Bell entre Alice y Bob pasa a otro estado de Bell entre Charles y Bob
