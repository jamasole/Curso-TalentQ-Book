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

# Números Complejos 

$ \newcommand{\bra}[1]{\langle #1|} $
$ \newcommand{\ket}[1]{|#1\rangle} $
$ \newcommand{\braket}[2]{\langle #1|#2\rangle} $
$ \newcommand{\i}{{ i}} $ 
$ \newcommand{\Hil}{{\mathbb H}} $

 

+++

```{contents}
:local:
:depth: 2
```

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
%run ../macro_tQ.py
import sys

sys.path.append('../')

import macro_tQ as tQ

import numpy as np
import scipy.linalg as la
from IPython.display import display,Markdown,Latex
import matplotlib.pyplot as plt
from qiskit.tools.visualization import array_to_latex
import copy 
```

+++ {"slideshow": {"slide_type": "slide"}}

## Introducción

+++ {"slideshow": {"slide_type": "slide"}}

La mecánica cuántica guarda una relación muy estrecha con los número complejos 

+++ {"slideshow": {"slide_type": "fragment"}}

Recuerda que el *cuadrado* de un número real, $a\in {\mathbb R}$ <u>siempre  es positivo</u>

+++ {"slideshow": {"slide_type": "fragment"}}



$$a^2 >0$$  


+++ {"slideshow": {"slide_type": "fragment"}}

por ejemplo $~~2^2 =  4~~$ pero también $~~(-2)^2 = 4~~$.

+++ {"slideshow": {"slide_type": "fragment"}}

Por eso, la raíz cuadrada de un número real <u>*sólo* existe si dicho número es positivo</u>.  


+++ {"slideshow": {"slide_type": "fragment"}}

Por ejemplo $\sqrt{4}=\pm 2$, mientras que $$\sqrt{-4}= \, ?$$

+++ {"slideshow": {"slide_type": "fragment"}}

**Pregunta**: $~~$ ¿cómo podríamos definir la <u>raíz cuadrada de un número real negativo</u>? 

+++ {"slideshow": {"slide_type": "slide"}}


**Respuesta**: para hacerlo, es necesario <u>ampliar el conjunto de los números reales</u>. 

+++

```{prf:definition} Número $i$
se postula la existencia de un <i>nuevo número</i>, $~i$, que es la solución única de la ecuación $~$

$$ i^2 = -1$$
```

+++ {"slideshow": {"slide_type": "fragment"}}

Equivalentemente podíamos haber requerido que  $\i = \sqrt{-1}$. 

+++ {"slideshow": {"slide_type": "fragment"}}

Con esto podemos ahora encontrar la raíz de cualquier número negativo. Por ejemplo, $-4$

$$
(2i)(2i) = 4 i^2 = - 4
$$

y $2i$ será la raíz buscada.

+++ {"slideshow": {"slide_type": "slide"}}

Con el número $i$ se opera igual que con los números reales


+++ {"slideshow": {"slide_type": "fragment"}}

$$ i + i = 2i $$

+++ {"slideshow": {"slide_type": "fragment"}}

$$ i - i = 0 $$

+++ {"slideshow": {"slide_type": "fragment"}}

$$ i + 2i = 3 i$$

+++ {"slideshow": {"slide_type": "fragment"}}

$$i^3  = i*i*i = i^2 * i = -i$$

+++ {"slideshow": {"slide_type": "fragment"}}

$$\frac{i}{i} = 1$$

+++ {"slideshow": {"slide_type": "slide"}}

Observar que el *inverso multiplicativo*  $1/i$ también es $-i$

\begin{eqnarray}
\frac{i}{i}  &=&   1 \nonumber\\
i (-i) &=& - i^2  =  1
\end{eqnarray}

+++ {"slideshow": {"slide_type": "fragment"}}

por tanto hay una identificación importante

$$ i^{-1} = \frac{1}{i} = -i $$

+++

```{admonition} En resumen
:class: warning

la solución al problema planteado consiste en <i>extender</i> el cuerpo de los números reales ${\mathbb R}$ al de los complejos ${\mathbb C}$ que, ahora, incluyen el número $i$

```

+++ {"slideshow": {"slide_type": "slide"}}

## Formas cartesiana y polar

+++ {"slideshow": {"slide_type": "slide"}}

### Forma Cartesiana





 Un *número complejo*, $z \in {\mathbb C}$, se representa en *forma cartesiana*  mediante <u>dos 
 números reales $x,y\in {\mathbb R}$</u> 
<br>

$$
z = x + \i y   ~~~~ \hbox{donde}~~~ 
\left\{\begin{array}{cl} x &\hbox{es la }   parte~ real\\
y & \hbox{es la } parte~imaginaria
\end{array}
\right.
$$

+++ {"slideshow": {"slide_type": "fragment"}}

Un número complejo se representa  en el **plano complejo**: la parte real  en el *eje horizontal*, y la parte imaginaria  en el *eje vertical*

+++ {"slideshow": {"slide_type": "slide"}}

#### Numeros complejos en python:
 
en *python*, el número  imaginario $\i$ se representa con la letra j.  

Añadiendo +0j convertimos un <i>float</i> en un <i>complex</i>

```{code-cell} ipython3
---
slideshow:
  slide_type: fragment
---
print(isinstance(1+0j,complex))
```

```{code-cell} ipython3
---
slideshow:
  slide_type: fragment
---

```

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
'''Forma cartesiana''' 
z = -3 + 1j

'''Extraemos las partes real e imaginaria'''
x=z.real
y=z.imag
print('z=x+iy=',x + 1j*y)

''' Representación en el plano complejo '''
'''esta función está en archivo: macro_CURSO.py  '''
''' Representación en el plano complejo '''
tQ.plot_2D_plane(left=-int(abs(x))-1,right=int(abs(x))+1,up=int(abs(y))+1,down=-int(abs(y))-1)
tQ.draw_vector(x,y,'b')
```

+++ {"slideshow": {"slide_type": "slide"}}




<div class="alert alert-block alert-info",text-align:center>
<p style="text-align: center;"> <p style="text-align: left ;color: navy;">  
<b> Teorema:</b> (Fórmula de Euler).
 Dado un ángulo $\theta \in (0,2\pi)$ las dos expresiones siguientes son equivalentes
<br>
<br>
$$
\cos\theta + i \sin \theta = e^{i\theta} 
$$

<details>
<summary><p style="text-align:right ; color:grey"> >>Demostración
</p></summary>
<br>
La demostración de la Fórmula de Euler viene de expandir ambos miembros en serie de Taylor en torno a $\theta = 0$ y comprobar que ambas series son iguales
<br> 
\begin{array}{rcl}
e^{i\theta} &=& 1 + \i\theta + \frac{1}{2}(\i\theta)^2 + \frac{1}{3!}(\i\theta)^3+\, ... \\
    &=& 1 -\frac{1}{2}\theta^2 ~+~ ... ~+~ \i \left(\theta - \frac{1}{3!} \theta^3+ \, ...\right) \\
    &=& \cos \theta  + i \sin \theta 
\end{array} 
</details>
</div>

+++ {"slideshow": {"slide_type": "slide"}}

El número $z=x + \i y$ se puede representar en *forma polar*
<br> 

$$
z = \rho e^{i\theta} = \rho (\cos\theta + \i \sin\theta)  
$$ 

+++ {"slideshow": {"slide_type": "fragment"}}

de donde obtenemos las componentes cartesianas 

$$x=\rho\cos\theta ~~,~~y=\rho\sin\theta$$

y por tanto $\rho^2 = x^2 + y^2$

+++ {"slideshow": {"slide_type": "slide"}}

Los números reales $\rho$ y $\theta$ se denominan *módulo* y *fase*.
<br>


+++ {"slideshow": {"slide_type": "fragment"}}

Las fases $\theta$ y $\theta+ 2\pi$ representan el *mismo* número complejo
<br>

$$
z = \rho e^{i\theta} = \rho e^{i(\theta + 2\pi)} 
$$

+++ {"slideshow": {"slide_type": "fragment"}}

Ello se debe a que las funciones $\cos \theta$ y $\sin\theta$ son periódicas

$$
\sin\theta  = \sin(\theta + 2\pi)~~~~~~~\cos\theta  = \cos(\theta +2\pi)
$$

+++ {"slideshow": {"slide_type": "slide"}}

Usar la forma polar es útil en situaciones en las que aparecen productos y potencias del número $i$

\begin{eqnarray*}
i &=& e^{i\pi/2} \\
-1 &=& e^{i\pi} \\ \rule{0mm}{8mm}
i &=& e^{3i\pi/2} =e^{-i\pi/2} \\  \rule{0mm}{8mm}
i^i &=& (e^{i\pi/2})^i = e^{i^2 \pi/2} = e^{-\pi/2}  \\  \rule{0mm}{8mm}
i^{2+i} &=& (e^{i\pi/2})^{(2+i)}= e^{i\pi/2(2+i)} = e^{i\pi} e^{-\pi/2} = -e^{-\pi/2}  \rule{0mm}{8mm}
\end{eqnarray*}

+++ {"slideshow": {"slide_type": "slide"}}

<div class="alert alert-block alert-success">
    <b>Ejercicio:</b>
    
haz una lista de 10 números (fases) $\theta_i, ~i=0,...,9$, equi-espaciadas entre 0 y 2$\pi$ y pinta los números complejos $z_i = r \exp(i \theta_i)$ y sus complejos conjugados, $z^*_i$.
</div>


+++ {"slideshow": {"slide_type": "slide"}}

### Conjugacion compleja
 
Todo número complejo, $z$, lleva *asociado* otro, $z^*$, denominado el *complejo conjugado* que se obtiene cambiando $\i \to -\i$
<br>

$$
z = x+\i y ~~~~\leftrightarrow~~~~ z^* = x - \i y \hspace{2cm}
$$
<br>

Es evidente que $(z^*)^* = z$, es una *involución*.

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
'''Conjugacion compleja'''
zc = z.conjugate()
print('z*=(x+iy)*=',z.real +1j*zc.imag)

''' Representación en el plano complejo '''
tQ.plot_2D_plane(left=-int(abs(x))-1,right=int(abs(x))+1,up=int(abs(y))+1,down=-int(abs(y))-1)
tQ.draw_vector(z.real,z.imag,'b')
tQ.draw_vector(zc.real,zc.imag,vcolor='r')

```

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
'''Forma polar'''
r=2
th=3.5
z = r*np.exp(th*1j)
print(z)
x=z.real
y=z.imag
print('z=r exp(i th)=',np.round(x,3) + np.round(y,3)*1j)


''' Dibujamos en el plano complejo '''
tQ.plot_2D_plane(left=-int(abs(x))-1,right=int(abs(x))+1,up=int(abs(y))+1,down=-int(abs(y))-1)
tQ.draw_vector(z.real,z.imag)
tQ.draw_vector(z.conjugate().real,z.conjugate().imag,vcolor='r')
```

+++ {"slideshow": {"slide_type": "fragment"}}

 En forma polar la *conjugación compleja* se obtiene *cambiando el signo* de la fase
<br>

$$
z = \rho e^{i\theta} ~\rightarrow ~z^* = \rho e^{-i\theta} = \rho \cos\theta - \i \rho \sin\theta \hspace{4cm}
$$

+++ {"slideshow": {"slide_type": "slide"}}

### Conversión entre formas cartesiana y polar

+++ {"slideshow": {"slide_type": "fragment"}}

La conversión de la representación *polar a cartesiana* es muy sencilla
gracias a las fórmula de Euler


<br>

$$
z = r e^{i\theta} = x + i y ~~~\hbox{ con }  ~~~\left\{\begin{array}{l} x=r \cos \theta \\ \rule{0mm}{4mm} y = r\sin \theta
\end{array} \right.
$$

<br>

+++ {"slideshow": {"slide_type": "slide"}}

La conversión inversa, *de cartesiana a polar* es un poco más delicada. Formalmente sería

<br>


$$
z = x + i y  = r e^{i\theta} ~~~ \hbox{ con } ~~~ \left\{\begin{array}{l} r=\sqrt{x^2+y^2} \\  \rule{0mm}{4mm} \theta = \arctan(y/x)
\end{array} \right.
$$

+++ {"slideshow": {"slide_type": "slide"}}

A la hora de la verdad hay que fijar el signo de la función $\arctan(y/x)$. La siguiente función examina esto mirando a los signos de $x$ y de $y$ para saber en qué cuadrante estamos.

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
'Conversión de Cartesianas a Polares'
def cartes2polar(z):

    r = np.abs(z)    
    y = z.imag
    x = z.real
    
    if r==0:
        print('el número 0+i0 no admite representación polar')
        th='indefinido'
    elif x==0 and y>0: 
        th=np.pi/2
    elif x==0 and y<0:
        th=3*np.pi/2
    elif x>0 and y>=0:
        th=np.arctan(y/x)
    elif x<0 and y>=0:
        th=np.arctan(-y/x)+np.pi/2
    elif x<0 and y<0:
        th=np.arctan(y/x)+np.pi
    elif x>0 and y<0:
         th=np.arctan(-y/x)+3*np.pi/2.       
            
    return r,th

#el signo correcto también se puede conseguir usando la funcion np.arctan2(x,y)
```

<div class="alert alert-block alert-success">
    <b>Ejercicio:</b>
calcula  la forma polar del número complejo $z = 4 + 3 i$ a mano y verificalo con la función que acabamos de definir

```{code-cell} ipython3
---
slideshow:
  slide_type: '-'
---
'A la inversa no es necesario definir ninguna funcion, ya que numpy directamente escribe un numero complejo en forma cartesiana'
z = 3*np.exp(1j*0.5)
print(np.round(z,2))
```

+++ {"slideshow": {"slide_type": "slide"}}

## Operaciones básicas


Los numeros complejos ${\mathbb C}$ forman una estructura matemática denominada *cuerpo*. Esto  quiere decir que admiten dos operaciones *internas*: la **suma** y la **multiplicación**. Vamos a estudiarlas por separado


+++ {"slideshow": {"slide_type": "skip"}}

### Suma

En representación *cartesiana* se *suman las partes real e imaginaria por separado*

<br>

$$
(a + \i b) + (c + \i d) = (a+c) + \i (b+d)\hspace{6cm}
$$ 

<br>



La resta es obvia, ya que $a,b,c,d$ pueden ser números negativos. 

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
'''Suma en cartesianas'''

z1 = 1+4j
z2 = 2+2j

'''Suma y resta'''
zs = z1+z2
zd = z1-z2

print('z1+z2=',zs)
print('z1-z2=',zd)
print('************************************')



tQ.plot_2D_plane(left=-2,right=4,up=7,down=-1) #cambiar las dimensiones para que encuadrar la figura
tQ.draw_vector(z1.real,z1.imag,'b')
tQ.draw_vector(z2.real,z2.imag,'b')
tQ.draw_vector(zs.real,zs.imag,vcolor='r')
tQ.draw_vector(zd.real,zd.imag,vcolor='g')
```

+++ {"slideshow": {"slide_type": "slide"}}

En *forma polar*, la suma de dos números complejos no admite ninguna simplificación, y deben transformarse primeramente a forma cartesiana, para sumarse. 

$$
z + w = \rho e^{i\theta} + \sigma e^{i\phi} = (\rho\cos\theta + \sigma\cos\phi) + i(\rho\sin\theta +  \sigma\sin\phi) \hspace{6cm}
$$

```{code-cell} ipython3
---
run_control:
  marked: true
slideshow:
  slide_type: slide
---
'python directamente escribe un numero complejo en forma cartesiana'
z1 = 3*np.exp(1j*0.5)
z2 = 1*np.exp(-1j*0.7)

'''Suma y resta'''
zs = z1+z2
zd = z1-z2

print('z1+z2=',np.round(zs,4))
print('z1-z2=',np.round(zd,4))
print('************************************')

tQ.plot_2D_plane(left=-2,right=4,up=3,down=-1) #cambiar las dimensiones para que encuadrar la figura
tQ.draw_vector(z1.real,z1.imag,'b')
tQ.draw_vector(z2.real,z2.imag,'b')
tQ.draw_vector(zs.real,zs.imag,vcolor='r')
tQ.draw_vector(zd.real,zd.imag,vcolor='g')
```

+++ {"slideshow": {"slide_type": "slide"}}

###  Multiplicación 



En *forma cartesiana* la multiplicación es complicada, debiendo multiplicarse todos los factores entre sí, y
teniendo en cuenta que $\i^2= -1$

$$
(a + \i b) (c + \i d) =ac +  a\i d +\i bc +\i^2 bd = (ab - bd) + \i(ac + bd)\hspace{6cm}
$$ 

Para hallar el producto de dos números complejos $z=r e^{i\theta}$ y $w=s e^{i\phi}$ escritos en forma polar, se multiplican los módulos y se suman las fases

$$
z w = r e^{i\theta} s e^{i\phi} = rs\,   e^{i(\theta + \phi)} 
$$

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
z1 = 3*np.exp(1j*0.5)
z2 = 1*np.exp(-1j*0.7)
'''Producto'''
print('z1*z2 = ', np.round(z1*z2,4))
print('z1**3 = ', np.round(pow(z1,6)))
print('************************************')
```

+++ {"slideshow": {"slide_type": "slide"}}

### Valor absoluto

El cuadrado de un número real $a\in {\mathbb R}$ es otro número real positivo $a^2 >0$. Ello nos permite definir
el valor absoluto $|a| = \sqrt{a^2}$ que es el mismo para $a$ y para $-a$.

+++ {"slideshow": {"slide_type": "fragment"}}

Esto no sucede con un número complejo $z$. $~$ En efecto,  

$$z^2 = x^2 - y^2 +2\i xy$$ 

es complejo.$~$ Sin embargo, el producto de un número por su conjugado es un número *real* y *positivo*

$$
z z^*  = (x + \i y) (x-\i y) = x^2 + y^2 >0
$$

+++ {"slideshow": {"slide_type": "fragment"}}

Ello nos permite definir el *valor absoluto* de un número complejo

$$
|z| = \sqrt{z z^*} = \sqrt{x^2 + y^2}
$$

+++ {"slideshow": {"slide_type": "slide"}}



El *valor absoluto* de una fase es 1

$$
|e^{i\theta}| = \sqrt{ e^{i\theta}   e^{-i\theta}}=\sqrt{ e^{i(\theta-\theta)}}=\sqrt{e^0} = 1
$$

+++ {"slideshow": {"slide_type": "fragment"}}

El *valor absoluto* de un número complejo con el *módulo*  escrito en forma polar

$$
|z| = \sqrt{zz^*} = \sqrt{\rho e^{i\theta} \rho e^{-i\theta}}=\sqrt{\rho^2} = \rho
$$

+++ {"slideshow": {"slide_type": "fragment"}}

<div class="alert alert-block alert-success">
    <b>Ejercicio:</b>
     Verifica el valor absoluto de un producto de números complejos es el producto de sus valores absolutos
</div>

+++ {"slideshow": {"slide_type": "slide"}}

### División


Al igual que la multiplicación, en forma cartesiana, la división **no es simple**. Sea $z = a+ \i b$ y $w=c+\i d$ 
<br>
<br>

$$
\frac{z}{w} = \frac{z}{w}\frac{w^*}{w^*} = \frac{( a+ \i b)(c-\i d)}{|w|^2} = \frac{ac+bd + \i(bc-ad)}{c^2+d^2}  = \frac{ac+bd}{c^2+d^2} +\i\frac{bc-ad}{c^2+d^2}
$$


+++ {"slideshow": {"slide_type": "fragment"}}

En forma polar la división es tan sencilla como la multiplicación. Se toma el  cociente de los módulos y la resta de las fases

$$
\frac{z}{w} = \frac{\rho e^{i\theta}}{\sigma e^{i\phi}} = \frac{\rho}{\sigma}
e^{i(\theta-\phi)}
$$

```{code-cell} ipython3
---
slideshow:
  slide_type: slide
---
'''Valor absoluto'''
print('|z1|=',abs(z1))
print('comprobación |z1|=',np.sqrt(z1*z1.conjugate()).real) 
print('************************************')


'''Division'''
print('z1/z2=',np.round(z1/z2,5))
print('comprobación z1/z2=', np.round(z1*z2.conjugate()/(z2*z2.conjugate()),5))
```

+++ {"slideshow": {"slide_type": "slide"}}

## Ejemplos y propiedades

+++ {"slideshow": {"slide_type": "-"}}

### Sumas nulas

En muchas ocasiones nos encontraremos la siguiente representación del numero cero (complejo) $0 = 0 + \i 0$

$$
\sum_{k=0}^{N-1} e^{2\pi \i k/N} =   e^{2\pi \i\, 0/N} +  e^{2\pi \i\, 1/N}  +~...~ +   e^{2\pi \i\, (N-2)/N}+   e^{2\pi \i\, (N-1)/N} ~=~  ~0
 $$

Para convencerse de que esta identidad es cierta vamos a representar los números complejos y su suma. Puedes cambiar $N$ y también multiplicar por un módulo constante
 

```{code-cell} ipython3
---
code_folding: [9]
slideshow:
  slide_type: slide
---
%run ../macro_CURSO.py

# cambiar el número N
N=9
rho=1

''' Creamos las fases'''
lista_de_fases=np.exp(2*np.pi*1j*np.array(range(N))/N)
#print('lista de fases =', np.round(lista_de_fases,2))


''' Dibujamos los números complejos '''
tQ.plot_2D_plane(fsize=(6,6))
for vec in rho*lista_de_fases:
    draw_vector(x=vec.real,y=vec.imag)

#draw_unit_circle()
plt.gca().add_patch(plt.Circle((0.,0.),1.,color='black',fill=False)) 


''' Calculamos la suma. '''
#print(lista_de_fases)
print(np.round(sum(rho*lista_de_fases),10))
```

Sea  $1\leq j \leq N-1$ un entero por el que multiplicamos todas las fases. El resultado es el mismo

$$
\sum_{k=0}^{N-1} e^{2\pi \i j k/N} =   e^{2\pi \i\, 0/N} +  e^{2\pi \i\, j/N}  +~...~ +   e^{2\pi \i\, j(N-2)/N}+   e^{2\pi \i\, j(N-1)/N} ~=~  ~0
 $$

+++

Sin embargo si $j = 0, N, 2N,... = 0\,\hbox{mod} N$, entonces la suma no se anula y su valor es igual a $N$. 

Tomemos por ejemplo $j=N$ 

$$
\sum_{k=0}^{N-1} e^{2\pi \i (3N) k/N} = \sum_{k=0}^{N-1} e^{2\pi \i  k}  =  \sum_{k=0}^{N-1} 1 =~  ~N
 $$

+++

<div class="alert alert-block alert-success">
<b>Ejercicio:</b>
 Modifica la lista de fases para convencerte de que todos los resultados anteriores son correctos
</div>

+++

Una manera de resumir todos los casos anteriores en una sola expresión involucra la función $\delta$ de Kronecker

$$
\delta_{ij} = \left\{ \begin{array}{rcl} 0 & \hbox{si} & i\neq 0 \\ 1 & \hbox{si} & i = j \end{array} \right.
$$

+++

Con ella podemos enunciar el siguiente resultado 

+++

<div class="alert alert-block alert-info",text-align:center>
<p style="text-align: Left; color: navy"> 
<br>
$$
\frac{1}{N}\sum_{k=0}^{N-1} e^{2\pi \i \, j k/N} =  \delta_{j\, 0{\rm mod} N}
$$
</p>
</div>

+++

que usaremos con profusión al estudiar la transformada de Fourier cuántica.

+++

### Desigualdad triangular

+++

El módulo de la suma de dos números complejos verifica que

$$
| z+w| \leq |z| + |w| 
$$

Donde la igualdad sólo se verifica cuando ambos números complejos son paralelos en el plano complejo. 

```{code-cell} ipython3
'''Comprueba que sólo cuando z1 y z2 son paralelos, se satura la desigualdad triangular'''

'''Suma en cartesianas'''
z1 = 1+2j

ang = 0. #el ángulo entre z1 y z2
z2 = z1*(1.2*np.exp(1j*ang))

'''Suma '''
zs = z1+z2

print('|z1|+|z2|=',abs(z1)+abs(z2))
print('|z1+z2|=',abs(z1+z2))


tQ.plot_2D_plane(left=-2,right=4,up=7,down=-1) #cambiar las dimensiones para que encuadrar la figura
tQ.draw_vector(z1.real,z1.imag,'b')
tQ.draw_vector(z2.real,z2.imag,'b')
tQ.draw_vector(zs.real,zs.imag,vcolor='r')
```
