# Solow-Swan dinámico

Simulador interactivo de la **ecuación fundamental del modelo de Solow-Swan**, construido a partir del capítulo 1 de *Apuntes de crecimiento económico* de Xavier Sala-i-Martin.

El objetivo es ver, instante por instante, cómo cada término de la ecuación empuja o frena la acumulación de capital per cápita:

$$\dot k = \underbrace{sAk^{\alpha}}_{\text{inversión bruta}} - \underbrace{\delta k}_{\text{depreciación}} - \underbrace{nk}_{\text{dilución}} - \underbrace{xk}_{\text{progreso técnico}}$$

## Qué incluye

| Archivo | Contenido |
|---|---|
| `index.html` | Simulador interactivo (abre en cualquier navegador o publícalo con GitHub Pages). |
| `R/solow_swan.R` | La misma simulación en R con `ggplot2`, para reproducir las gráficas en RStudio. |

## Qué muestra el simulador

- **Ecuación viva**: el valor de cada término (sAk^α, δk, nk, xk y k̇) en el momento *t* que elijas, con una barra que compara la inversión bruta contra la inversión necesaria para mantener *k* constante.
- **Diagrama de Solow (gráfico 1.1)**: f(k), s·f(k) y (δ+n+x)k, con el estado estacionario k\* y el capital de la regla de oro.
- **Dinámica de transición (gráfico 1.10)**: curva de ahorro CA = sAk^(α−1) y curva de depreciación CD; la distancia entre ambas es γ_k.
- **Descomposición en el tiempo** de todos los términos de la ecuación.
- **Trayectorias** de k, y y c, y **tasas de crecimiento** γ_k, γ_y = αγ_k, junto con la aproximación lineal [1.28].
- **Choques permanentes** en s, A, n, δ o x en el momento T, con escenarios precargados de la sección 1.4 (gráficos 1.11, 1.12 y 1.13), convergencia (1.6–1.7) y ahorro excesivo frente a la regla de oro.
- **Tabla de indicadores** antes y después del choque: k\*, y\*, c\*, k_oro, c_oro, s_oro = α, β\* = (1−α)(δ+n+x) y la vida media ln2/β\*.

## Ecuaciones del capítulo 1 que se usan

| Ecuación | Expresión |
|---|---|
| [1.15] Ecuación fundamental | k̇ = sAk^α − (δ+n)k |
| [1.16] Estado estacionario | k\* = (sA/(δ+n))^(1/(1−α)) |
| [1.17] Consumo de estado estacionario y regla de oro | c* = f(k*) − (δ+n)k*; máximo cuando f′(k_oro) = δ+n ⇒ s_oro = α |
| [1.20′] Tasa de crecimiento | γ_k = sAk^(−(1−α)) − (δ+n) |
| [1.24] Con progreso técnico | k̂̇ = sf(k̂) − (δ+n+x)k̂ |
| [1.27]–[1.28] Convergencia | β\* = (1−α)(δ+n), γ_k ≈ −β\*[log k − log k\*] |

Con los valores del libro (α = 0.30, δ = 0.10, n = 0.01) el simulador reproduce β\* ≈ 7.7% anual y una vida media de unos 9 años.

## Cómo publicarlo con GitHub Pages

1. Crea un repositorio nuevo en GitHub (por ejemplo `solow-swan-dinamico`).
2. Sube `index.html`, `README.md` y la carpeta `R/` a la rama `main`.
3. En el repositorio entra a **Settings → Pages**, en *Source* elige **Deploy from a branch**, rama `main`, carpeta `/ (root)` y guarda.
4. En uno o dos minutos el simulador estará en `https://TU-USUARIO.github.io/solow-swan-dinamico/`.

## Notas técnicas

- La ecuación diferencial se resuelve numéricamente con Runge-Kutta de cuarto orden (paso de 0.05 años).
- Si x > 0, k y y se interpretan por trabajador efectivo; el producto per cápita crece a γ_y + x.
- Dependencias cargadas por CDN: Chart.js 4.4.1 y KaTeX 0.16.9. El script de R usa `ggplot2`, `tidyr` y `dplyr` (R ≥ 4.1 por el operador `|>`).

## Referencia

Sala-i-Martin, X. *Apuntes de crecimiento económico*, capítulo 1: "El modelo neoclásico de crecimiento de Solow-Swan".
