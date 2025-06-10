# Simulación de Fuga de Gas del Casing al Anular

Este script en Python simula una fuga de gas desde el casing hacia el anular en un pozo petrolero. El modelo tiene en cuenta la geometría del pozo, propiedades de los fluidos, condiciones iniciales, y calcula la evolución de presiones, caudales y propagación de la fuga en el tiempo.

---

## 📂 Estructura

El script está dividido en las siguientes secciones:

1. **Parámetros de entrada**: Se configuran todas las condiciones iniciales, dimensiones del pozo y propiedades físicas.
2. **Funciones auxiliares**: Funciones de cálculo para el factor Z, propiedades pseudo-críticas, fricción, etc.
3. **Cálculos iniciales**: Cálculo de presiones estáticas para casing y anular.
4. **Inicialización de simulación**: Se preparan los datos y validaciones para la simulación dinámica.
5. **Bucle de simulación**: Se calcula paso a paso la evolución temporal del sistema.

---

## ⚙️ Requisitos

Instalar las siguientes librerías de Python (si aún no las tenés):

```bash
pip install numpy pandas matplotlib
