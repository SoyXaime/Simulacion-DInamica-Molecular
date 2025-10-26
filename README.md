# 🔬 Simulación de Dinámica Molecular con Potencial Lennard-Jones

Este proyecto implementa una simulación de **Dinámica Molecular (DM)** para estudiar la evolución y el equilibrio de un sistema de partículas interactuantes mediante el potencial **Lennard-Jones**, integrando las ecuaciones de movimiento con el método **Velocity-Verlet**.  

Se incluyen herramientas para:
- 🧱 **Generar la red inicial**
- ⚙️ **Ejecutar la simulación** (con corrección de energía)
- 💾 **Guardar datos de velocidades y energías**
- 📊 **Analizar y visualizar resultados** mediante **Jupyter Notebooks**

---

## 🧠 Modelo Físico

El potencial Lennard-Jones viene dado por:


$$U(r) = 4\varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r}\right)^6 \right]$$


Este describe:
- Repulsión fuerte a cortas distancias  
- Atracción débil a distancias mayores  

Las ecuaciones de movimiento se integran con el método **Velocity-Verlet**, que conserva bien la energía y es estándar en simulaciones MD.

---

## 📂 Estructura del Repositorio

Simulacion-DInamica-Molecular/
├─ RED/                    # Código Fortran para generar la red inicial
├─ Equilibrio/     	   # Ajuste y estabilización de energía
├─ Dinamica/               # Simulación principal + guardado de datos en .txt
├─ notebooks/              # Análisis, gráficas y comprobación de equilibrio
└─ resultados/             # Resultados generados (NO versionados)

> Los archivos de salida (`.txt`, `.bin`, etc.) se guardan en `resultados/` para mantener el repo limpio.

---

## ⚙️ Compilación y Ejecución

### 1) Compilar (con `gfortran`)
		gfortran *.f90 -o crea_red

		gfortran *.f90 -o equilibrio

		gfortran -O3 -fopenmp *.f90 -o condiciones_eq_omp 
		
		gfortran -O3 -fopenmp *.f90 -o dinamica 


### 2) Crear carpeta para resultados:
		mkdir -p /resulados

### 3) Comandos de ejecución: (Ejecutar cada programa dentro de su carpeta)
		./crea_red

		./equilibrio

		OMP_NUM_THREADS=3 ./condiciones_eq_omp 
	
		OMP_NUM_THREADS=3 ./dinamica 
	

📊 Análisis de resultados
Ver notebooks




