# Método de Lattice Boltzmann (LBM)  

El **método de Lattice Boltzmann (LBM)** es una herramienta numérica alternativa dentro de la **dinámica de fluidos computacional (CFD)**, basada en una **representación mesoscópica** de los fluidos. En lugar de resolver directamente las ecuaciones de **Navier–Stokes**, este método modela la evolución espacio-temporal de una **función de distribución de partículas** sobre una red discreta (*lattice*). En cada nodo de la red, las partículas se desplazan en direcciones discretas y colisionan localmente, generando de forma emergente el comportamiento macroscópico del fluido.  

El LBM **simula el flujo de fluidos a partir de interacciones microscópicas**, discretizando el espacio de velocidades y el dominio físico, lo que permite un tratamiento **computacionalmente eficiente** en comparación con los métodos CFD tradicionales. Esta eficiencia, junto con su formulación local, lo hace ideal para la **paralelización** y para su ejecución en **GPU** o **supercomputadoras**.  

Entre sus ventajas está su capacidad para manejar **geometrías complejas**, **fronteras irregulares** y **flujos multifásicos**, así como incorporar otros efectos físicos como **reacciones químicas, transferencia de calor, turbulencia o acústica**. Gracias a su estructura matemática, el método puede capturar fenómenos **multiescalares** que van desde lo microscópico hasta lo macroscópico con notable precisión.  

### OpenLB
**OpenLB** es una de las implementaciones más reconocidas de código abierto del método LBM. Está desarrollada en **C++** y diseñada para ejecutarse de manera eficiente en **arquitecturas paralelas**. Su enfoque modular permite a los usuarios definir **condiciones de frontera personalizadas**, **modelos multifísicos** (como flujo térmico o multifásico) y trabajar con **geometrías irregulares**.  

OpenLB es muy utilizado en investigación académica para simular **flujos en medios porosos**, **transferencia de calor**, **microfluidos** y **geotermia**, entre otros. Además, cuenta con una comunidad activa y documentación extensa que facilita su aprendizaje y extensión. Es especialmente apreciado por su **transparencia del código**, su **flexibilidad**, y su uso en entornos **educativos y de investigación**.  

### Palabos  
**Palabos (Parallel Lattice Boltzmann Solver)** es otra implementación **libre y de código abierto** del método LBM. Se distingue por ofrecer una **infraestructura de alto nivel** para simular sistemas fluidodinámicos complejos con facilidad, sin necesidad de modificar el núcleo del código. Al igual que OpenLB, está escrito en **C++** y utiliza técnicas avanzadas de **paralelización** que permiten ejecutar simulaciones en **clusters HPC** y en **GPU**.  

Palabos proporciona un entorno muy versátil para el modelado de **interacciones fluido–estructura**, **flujos multifásicos**, **transportes reactivos**, **dinámica térmica** y problemas con **interfaces complejas**. Además, está diseñado con una arquitectura **modular y orientada a objetos**, lo que facilita agregar modelos físicos personalizados. Gracias a su robustez y capacidad de manejo de datos grandes, Palabos también se emplea en **industria**, especialmente en **aerodinámica**, **biomedicina** y **energías renovables**.  

### Comparación general  
Tanto **OpenLB** como **Palabos** son proyectos abiertos, bien mantenidos y ampliamente usados en la comunidad científica. Las diferencias principales radican en su **estructura interna del código** y en su **nivel de abstracción**:  
- **OpenLB** ofrece un entorno más **transparente y didáctico**, ideal para investigar o desarrollar nuevas formulaciones del LBM.  
- **Palabos**, en cambio, ofrece una **infraestructura más madura y modular**, enfocada en la creación rápida de aplicaciones de simulación complejas y extensas.  
Ambos aprovechan al máximo el **paralelismo computacional** y son adecuados para **flujos complejos y multifísicos**.  

En resumen, el **método de Lattice Boltzmann** es una poderosa alternativa numérica dentro de la CFD moderna, y tanto **OpenLB** como **Palabos** representan dos de sus implementaciones más relevantes, combinando precisión física, flexibilidad numérica y eficiencia computacional para abordar problemas fluidodinámicos de gran complejidad.  
